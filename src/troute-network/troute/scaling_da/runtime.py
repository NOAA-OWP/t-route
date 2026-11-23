"""Setup and per-loop step of area-scaled streamflow DA for main_v04."""

from __future__ import annotations

import logging
import math
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from troute.scaling_da.spread import MIN_FLOW_CMS, corrections, finite_corrections, make_plan
from troute.scaling_da.trees import build_trees

LOG = logging.getLogger(__name__)

DEPTH_EXPONENT = 0.6  # wide-channel Manning: h ~ Q ** 0.6
DEPTH_RATIO_CAP = 2.0  # the flow ratio used for depth is clipped to [1/2, 2]
# The spread is computed and written in blocks of this many timesteps, which bounds
# its memory whatever the loop length. Within a block, a segment that receives a
# correction is rewritten and floored at 0; other values stay as routed.
BLOCK_STEPS = 24

# Positions in a routing result tuple, and values per timestep in flowveldepth.
_IDS, _FVD, _DA, _NUDGE = 0, 1, 3, 8
_QVD = 3


def read_site_list(path: Path) -> List[str]:
    """Gage ids, one per line; blank lines and # comments are skipped."""
    lines = Path(path).read_text().splitlines()
    return [s for s in (line.split("#", 1)[0].strip() for line in lines) if s]


def read_theta_csv(path: Path) -> Dict[str, float]:
    """Gage id to theta from a two-column CSV with an optional header row."""
    theta: Dict[str, float] = {}
    for n, line in enumerate(Path(path).read_text().splitlines(), start=1):
        cells = [c.strip() for c in line.split("#", 1)[0].split(",")]
        if cells == [""]:
            continue
        if len(cells) < 2 or not cells[0]:
            raise ValueError(f"{path}:{n}: expected 'gage_id,theta', got {line!r}")
        try:
            value = float(cells[1])
        except ValueError:
            if not theta:
                continue  # header
            raise ValueError(f"{path}:{n}: theta must be a number, got {cells[1]!r}") from None
        if cells[0] in theta:
            raise ValueError(f"{path}:{n}: duplicate gage id {cells[0]!r}")
        if not (math.isfinite(value) and value > 0):
            raise ValueError(f"{path}:{n}: theta must be finite and > 0, got {cells[1]!r}")
        theta[cells[0]] = value
    if not theta:
        raise ValueError(f"{path}: no rows")
    return theta


def _observed_in_file(path: Path, qc_threshold: float) -> Set[str]:
    from troute.nhd_io import _read_timeslice_file

    obs, quality = _read_timeslice_file(path)
    if obs is None or len(obs) == 0:
        return set()
    o = obs.to_numpy(dtype=np.float64)
    q = quality.reindex(index=obs.index, columns=obs.columns).to_numpy(dtype=np.float64)
    # same screening as nhd_io.get_obs_from_timeslices
    ok = np.isfinite(o) & (o > 0) & np.isfinite(q) & (q >= 0) & (q <= 1) & (q >= qc_threshold)
    return {str(s).strip() for s, row in zip(obs.index, ok) if row.any()}


def observed_sites(da_config: Mapping[str, Any], da_sets: Sequence[Mapping[str, Any]], cpu_pool: int) -> Set[str]:
    """Gage ids with at least one usable observation in any TimeSlice file of the run."""
    folder = da_config.get("usgs_timeslices_folder")
    names = {str(f) for ds in da_sets or () for f in ((ds or {}).get("usgs_timeslice_files") or ())}
    files = [Path(folder) / n for n in sorted(names)] if folder else []
    qc = da_config.get("qc_threshold")
    qc = 1.0 if qc is None else float(qc)
    if cpu_pool > 1 and len(files) > 1:
        from joblib import Parallel, delayed

        parts = Parallel(n_jobs=cpu_pool)(delayed(_observed_in_file)(f, qc) for f in files)
    else:
        parts = [_observed_in_file(f, qc) for f in files]
    observed = set().union(*parts)
    if not observed:
        raise ValueError(
            f"scaling DA: the run's {len(files)} USGS TimeSlice file(s) in {folder!r} hold no usable "
            "observations; check usgs_timeslices_folder and the simulation period"
        )
    return observed


def _check_keys(kind: str, keys: Iterable[str], known: Set[str]) -> None:
    """Raise if no key matches the domain; warn about the ones that don't."""
    keys = set(keys)
    unused = keys - known
    if keys and unused == keys:
        raise ValueError(
            f"streamflow_scaling_parameters.theta.{kind}: none of its {len(keys)} key(s) match this "
            f"domain, e.g. {sorted(keys)[:5]} vs {sorted(known)[:5]}; check the id format"
        )
    if unused:
        LOG.warning("scaling DA: theta.%s keys matching nothing in this domain: %s", kind, sorted(unused)[:10])


class ScalingDA:
    """Static setup (sources, holdout, trees) and the per-loop upstream spread."""

    def __init__(self, trees, held_out: Set[str], held_out_segs: Set[int], no_observations: Set[str]):
        self.trees = trees
        self.held_out = held_out
        self.held_out_segs = held_out_segs
        self.no_observations = no_observations
        self.loops: List[Dict[str, Any]] = []
        # every segment the spread reads: tree members and pruned neighbours
        self.segments = np.unique(np.concatenate([trees.seg, trees.pruned_seg]))
        self.plan = make_plan(trees, {int(s): i for i, s in enumerate(self.segments)})
        self._tree_of_gage = {int(trees.seg[r]): int(trees.tree[r]) for r in trees.roots}

    @classmethod
    def build(cls, network, da_config: Mapping[str, Any], da_sets: Sequence[Mapping[str, Any]] = (),
              *, cpu_pool: int = 1, observed: Optional[Set[str]] = None) -> Optional["ScalingDA"]:
        """Set up scaling DA for a built HYFeatures network, or return None if it is off.

        Parameters
        ----------
        network : HYFeaturesNetwork
        da_config : mapping
            data_assimilation_parameters.
        da_sets : sequence
            The run's TimeSlice file sets, one per loop.
        cpu_pool : int
            Workers for reading the TimeSlice files.
        observed : set of str, optional
            Gage ids with observations; skips reading the TimeSlice files.
        """
        stream_da = da_config.get("streamflow_da") or {}
        if not stream_da.get("streamflow_scaling", False):
            return None
        start = time.perf_counter()
        opts = {"max_reach_km": 200.0, "holdout_sites_file": None, "theta": None}
        given = stream_da.get("streamflow_scaling_parameters")
        opts.update(given or {})
        theta_opts = {"default": 0.77, "by_vpu": None, "per_tree_csv": None}
        theta_opts.update(opts["theta"] or {})
        attrs = getattr(network, "segment_attributes", None)
        if attrs is None or "total_dasqkm" not in attrs.columns:
            raise ValueError("area-scaled streamflow DA needs total_dasqkm on the routed hydrofabric layer")

        crosswalk = {int(seg): str(site) for seg, site in ((network.gages or {}).get("gages") or {}).items()}
        held_out: Set[str] = set()
        holdout_file = opts["holdout_sites_file"]
        if holdout_file:
            held_out = set(read_site_list(holdout_file))
            unknown = sorted(held_out - set(crosswalk.values()))
            if unknown:
                raise ValueError(f"holdout_sites_file lists {len(unknown)} id(s) that are not gages of this "
                                 f"domain, e.g. {unknown[:5]}")
        if observed is None:
            observed = observed_sites(da_config, da_sets, cpu_pool or 1)
        observed = {str(s).strip() for s in observed}

        rconn = network.reverse_network
        waterbodies = network.waterbody_connections or {}
        lakes = {int(k) for k in waterbodies} | {int(v) for v in waterbodies.values()}
        diffusive = set()
        for tw, d in (network.diffusive_network_data or {}).items():
            diffusive |= {int(tw)} | {int(s) for s in (d or {}).get("mainstem_segs", ())}

        dropped: Dict[str, str] = {}
        routed: Dict[str, int] = {}
        for seg, site in sorted(crosswalk.items()):
            if seg in diffusive or seg not in rconn:
                dropped[site] = "not routed by MC"
            else:
                routed[site] = seg
        candidates = {s: g for s, g in routed.items() if g not in lakes}
        dropped.update({s: "on a lake" for s in routed if s not in candidates})
        sources = {s: g for s, g in candidates.items() if s in observed and s not in held_out}
        no_obs = {s for s in candidates if s not in observed and s not in held_out}
        # every gage nudging corrects stops the trees, so no flow is corrected twice
        stops = {seg for seg, site in crosswalk.items() if site in observed and site not in held_out}
        stops |= lakes | diffusive

        vpu = attrs["vpuid"].dropna().astype(str) if "vpuid" in attrs.columns else pd.Series(dtype=str)
        vpu_theta: Dict[str, float] = {}
        for key, value in (theta_opts["by_vpu"] or {}).items():
            vpu_theta[str(key)] = float(value)
        csv_file = theta_opts["per_tree_csv"]
        site_theta = read_theta_csv(csv_file) if csv_file else {}
        _check_keys("per_tree_csv", site_theta, set(crosswalk.values()))
        _check_keys("by_vpu", vpu_theta, {vpu.get(g) for g in routed.values()} - {None})
        fallback = float(theta_opts["default"])
        theta = {s: site_theta.get(s, vpu_theta.get(vpu.get(g), fallback)) for s, g in sources.items()}

        df = network.dataframe
        length_km = pd.Series(np.asarray(df["dx"].values, dtype=np.float64) / 1000.0, index=df.index)
        area_km2 = attrs["total_dasqkm"].astype(np.float64)
        trees, no_area = build_trees(
            rconn, sources, area_km2[~area_km2.index.duplicated()], length_km[~length_km.index.duplicated()],
            stops, theta, float(opts["max_reach_km"]),
        )
        dropped.update(no_area)

        held_out_segs = {seg for seg, site in crosswalk.items() if site in held_out}
        rt = cls(trees, held_out, held_out_segs, no_obs)
        if dropped:
            LOG.warning("scaling DA: %d gage(s) get no tree: %s", len(dropped), sorted(dropped.items())[:10])
        LOG.info("scaling DA: %d source gage(s), %d tree(s) with %d segment(s); %d held out, %d without "
                 "observations; setup %.2f s", len(sources), len(trees.sites), trees.seg.size,
                 len(held_out), len(no_obs), time.perf_counter() - start)
        return rt

    def nudging_inputs(self, usgs_df: pd.DataFrame, lastobs_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """usgs_df and lastobs_df without the held-out gages (both indexed by segment)."""
        def drop(df):
            if not self.held_out_segs or df is None or df.empty:
                return df
            keep = ~df.index.isin(list(self.held_out_segs))
            return df if keep.all() else df.loc[keep]

        return drop(usgs_df), drop(lastobs_df)

    def _locate(self, run_results) -> Tuple[np.ndarray, np.ndarray]:
        """Result index and row of every segment the spread reads."""
        ids = [np.asarray(r[_IDS], dtype=np.int64).ravel() for r in run_results]
        where = pd.DataFrame({
            "result": np.repeat(np.arange(len(ids)), [a.size for a in ids]),
            "row": np.concatenate([np.arange(a.size) for a in ids]),
        }, index=np.concatenate(ids))
        where = where[where.index.isin(self.segments)]
        if not where.index.is_unique:
            raise RuntimeError("scaling DA: a tree segment appears in more than one routing result")
        where = where.reindex(self.segments)
        if where["result"].isna().any():
            missing = self.segments[where["result"].isna().to_numpy()]
            raise RuntimeError(f"scaling DA: {missing.size} tree segment(s) have no routing result, "
                               f"e.g. {missing[:5].tolist()}")
        return where["result"].to_numpy(dtype=np.int64), where["row"].to_numpy(dtype=np.int64)

    def apply(self, run_results, *, dt: float, final: bool) -> bool:
        """Spread this loop's gage corrections upstream, rewriting run_results in place.

        Returns True when final is set and a flow changed, so the caller takes
        q0 again and the corrected state reaches the next run.
        """
        start = time.perf_counter()
        plan = self.plan
        if plan.roots.size == 0 or not run_results:
            return False
        result_of, row_of = self._locate(run_results)
        groups = [(k, np.flatnonzero(result_of == k)) for k in np.unique(result_of)]
        nts = np.asarray(run_results[groups[0][0]][_FVD]).shape[1] // _QVD
        q = np.empty((self.segments.size, nts), dtype=np.float32)
        for k, idx in groups:
            q[idx] = np.asarray(run_results[k][_FVD])[row_of[idx], 0::_QVD]

        # nudge column 0 is the initial state; columns 1..nts match the published steps
        dq_gage = np.zeros((plan.roots.size, nts), dtype=np.float64)
        for r in run_results:
            gage_segs = np.asarray(r[_DA][0]).astype(np.int64).ravel() if len(r[_DA]) else ()
            for i, seg in enumerate(gage_segs):
                t = self._tree_of_gage.get(int(seg))
                if t is not None:
                    dq_gage[t] = np.asarray(r[_NUDGE])[i, 1:]
        # the confluence split reads the gage's flow without the nudge, as NGWPC does
        gage_rows = plan.rows[plan.roots]
        gage_flow = q[gage_rows].astype(np.float64) - dq_gage
        q[gage_rows] = gage_flow.astype(q.dtype)

        old = q[plan.rows]
        dq_spread = finite_corrections(dq_gage)
        changed = np.zeros(plan.rows.size, dtype=bool)
        volume = 0.0
        for t0 in range(0, nts, BLOCK_STEPS):
            t1 = min(nts, t0 + BLOCK_STEPS)
            inc = corrections(plan, q[:, t0:t1], dq_spread[:, t0:t1], gage_flow[:, t0:t1])
            live = np.flatnonzero(np.any(inc != 0.0, axis=1))
            if live.size == 0:
                continue
            before = old[live, t0:t1].astype(np.float64)
            stored = np.maximum(before + inc[live], 0.0).astype(q.dtype)
            q[plan.rows[live], t0:t1] = stored  # later blocks read only their own columns
            applied = stored.astype(np.float64) - before
            changed[live] |= np.any(applied != 0.0, axis=1)
            volume += float(applied.sum())
        nodes = np.flatnonzero(changed)

        is_changed = np.zeros(self.segments.size, dtype=bool)
        is_changed[plan.rows[nodes]] = True
        node_at_row = np.full(self.segments.size, -1, dtype=np.int64)
        node_at_row[plan.rows] = np.arange(plan.rows.size)
        for k, idx in groups:
            idx = idx[is_changed[idx]]
            if idx.size == 0:
                continue
            fvd, rows = run_results[k][_FVD], row_of[idx]
            q_new, q_old = q[idx], old[node_at_row[idx]].astype(np.float64)
            depth = fvd[rows, 2::_QVD].astype(np.float64)
            live = (q_old > MIN_FLOW_CMS) & (q_new > MIN_FLOW_CMS) & (depth > 0.0)
            if live.any():
                ratio = np.clip(np.where(live, q_new / np.where(live, q_old, 1.0), 1.0),
                                1.0 / DEPTH_RATIO_CAP, DEPTH_RATIO_CAP)
                fvd[rows, 2::_QVD] = (depth * np.where(live, ratio ** DEPTH_EXPONENT, 1.0)).astype(fvd.dtype)
            fvd[rows, 0::_QVD] = q_new.astype(fvd.dtype)

        stats = dict(loop=len(self.loops), final=bool(final),
                     trees=int(np.count_nonzero(np.any(dq_gage != 0.0, axis=1))),
                     segments=int(nodes.size), volume_m3=volume * float(dt),
                     seconds=time.perf_counter() - start)
        self.loops.append(stats)
        LOG.info("scaling DA loop %d%s: %d tree(s) with a gage correction, %d segment(s) changed, "
                 "injected volume %.4g m3, %.2f s", stats["loop"], " (final, enters the state)" if final else "",
                 stats["trees"], stats["segments"], stats["volume_m3"], stats["seconds"])
        return bool(final and nodes.size > 0)
