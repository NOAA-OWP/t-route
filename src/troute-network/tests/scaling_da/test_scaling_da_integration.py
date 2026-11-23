"""Area-scaled streamflow DA through main_v04 on a small HFv4.0 subset.

Skipped unless TROUTE_HFV4_FIXTURE points at the subset gpkg. The gpkg is not
committed; it is cut from an HFv4.0 file with the hfv4_subset tool in
t-route-utils (upstream of fp-16820 and fp-8326 in VPU01). Forcing and USGS
TimeSlice files are synthetic and written to a temporary directory: lateral
inflow proportional to each nexus's area with a storm pulse, and observations
equal to 1.3 times the open-loop flow at each gage.
"""

from __future__ import annotations

import os
import sqlite3
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import pytest

FIXTURE = os.environ.get("TROUTE_HFV4_FIXTURE", "")

pytestmark = pytest.mark.skipif(
    not (FIXTURE and Path(FIXTURE).is_file()), reason="set TROUTE_HFV4_FIXTURE to an HFv4.0 subset gpkg"
)

T0 = datetime(2019, 10, 28, 0, 0)
DT = 300
NTS = 144  # 12 h
LOOP_H = 6  # two loops: the first keeps the uncorrected state, the second hands it off
FACTOR = 1.3

COLUMNS = {
    "flowpath": {"key": "flowpath_id", "downstream": "flowpath_toid"},
    "flowline": {"key": "flowline_id", "downstream": "flowline_toid"},
}
CHANNEL = {
    "dx": "lengthm", "n": "n", "ncc": "nCC", "s0": "So", "bw": "BtmWdth", "waterbody": "lake_id",
    "tw": "TopWdth", "twcc": "TopWdthCC", "musk": "MusK", "musx": "MusX", "cs": "ChSlp", "alt": "alt",
}


# ----------------------------------------------------------------- inputs


def write_forcing(gpkg: str, path: Path) -> None:
    """Lateral inflow at every nexus, hourly over the run plus 1 h."""
    import xarray as xr

    con = sqlite3.connect(gpkg)
    rows = con.execute("select flowpath_toid, sum(areasqkm) from flowpaths group by flowpath_toid").fetchall()
    con.close()
    ids = np.array([r[0].encode() for r in rows], dtype="S")
    area = np.array([float(r[1] or 0.0) for r in rows])
    hours = np.arange(NTS * DT // 3600 + 2)
    # 0.002 m3/s/km2 base plus a 0.02 m3/s/km2 storm pulse peaking at 4 h
    rate = 0.002 + 0.02 * np.exp(-0.5 * ((hours - 4.0) / 1.5) ** 2)
    runoff = (rate[:, None] * area[None, :]).astype(np.float32)
    times = pd.to_datetime([T0 + timedelta(hours=int(h)) for h in hours])
    ds = xr.Dataset({"runoff_rate": (("time", "feature_id"), runoff)}, coords={"time": times, "feature_id": ids})
    ds["runoff_rate"].attrs["units"] = "m3 s-1"
    ds.to_netcdf(
        path,
        encoding={
            "time": {"units": "minutes since 1970-01-01 00:00:00", "dtype": "int64"},
            "runoff_rate": {"_FillValue": -9999.0},
        },
    )


def char_array(strings: Sequence[str], width: int) -> np.ndarray:
    """(n, width) single bytes for a netCDF char variable, right-justified like TimeSlice ids."""
    raw = b"".join(s.encode("ascii").rjust(width) for s in strings)
    return np.frombuffer(raw, dtype="S1").reshape(len(strings), width)


def write_timeslices(folder: Path, obs: Dict[str, np.ndarray]) -> None:
    """One 15-min TimeSlice file per observation time; obs[site][k] is at T0 + (k+1)*DT.

    A file at T0 with zero discharge (screened out by the QC) anchors the
    observation time axis at T0; without it the first loop's observations are
    applied one cadence early.
    """
    import netCDF4

    os.makedirs(folder, exist_ok=True)
    step = 900 // DT
    for k in range(-1, NTS, step):
        t = T0 + timedelta(seconds=(k + 1) * DT)
        if k < 0:
            sites = sorted(obs)
            values = np.zeros(len(sites), dtype=np.float32)
        else:
            sites = [s for s in sorted(obs) if np.isfinite(obs[s][k])]
            values = np.array([obs[s][k] for s in sites], dtype=np.float32)
        name = folder / f"{t:%Y-%m-%d_%H:%M:%S}.15min.usgsTimeSlice.ncdf"
        with netCDF4.Dataset(name, "w", format="NETCDF4") as ds:
            for dim, size in (("stationIdInd", len(sites)), ("stationIdStrLen", 15), ("timeStrLen", 19)):
                ds.createDimension(dim, size)
            sid = ds.createVariable("stationId", "S1", ("stationIdInd", "stationIdStrLen"))
            tim = ds.createVariable("time", "S1", ("stationIdInd", "timeStrLen"))
            q = ds.createVariable("discharge", "f4", ("stationIdInd",))
            qual = ds.createVariable("discharge_quality", "i2", ("stationIdInd",))
            if sites:
                sid[:] = char_array(sites, 15)
                tim[:] = char_array([f"{t:%Y-%m-%d_%H:%M:%S}"] * len(sites), 19)
                q[:] = values
                qual[:] = 100
            ds.sliceCenterTimeUTC = f"{t:%Y-%m-%d_%H:%M:%S}"
            ds.sliceTimeResolutionMinutes = "15"


# ----------------------------------------------------------------- runs


@dataclass
class Run:
    flow: pd.DataFrame  # segment x timestep, all loops
    depth: pd.DataFrame
    q0: List[pd.DataFrame]  # every q0 the driver took, in order
    network: object
    runtime: Optional[object]
    loops: int = 0


def run_troute(tmp: Path, mode: str, *, timeslices: Optional[Path], scaling: Optional[dict] = None,
               method: str = "serial", block_off: bool = False) -> Run:
    yaml = pytest.importorskip("yaml")
    main = pytest.importorskip("nwm_routing.__main__")
    from troute.AbstractNetwork import AbstractNetwork
    from troute.scaling_da.runtime import ScalingDA

    sda = {"streamflow_nudging": True, "diffusive_streamflow_nudging": False}
    if scaling is not None:
        sda.update(streamflow_scaling=True, streamflow_scaling_parameters=scaling)
    elif block_off:
        sda.update(streamflow_scaling=False, streamflow_scaling_parameters={"max_reach_km": 50.0})
    no_reservoir_da = {
        "reservoir_persistence_da": {"reservoir_persistence_usgs": False, "reservoir_persistence_usace": False},
        "reservoir_rfc_da": {"reservoir_rfc_forecasts": False},
    }
    da_block = dict(usgs_timeslices_folder=str(timeslices), streamflow_da=sda, reservoir_da=no_reservoir_da)
    topology = dict(
        supernetwork_parameters={"geo_file_path": FIXTURE, "columns": {**COLUMNS[mode], **CHANNEL}},
        waterbody_parameters={"break_network_at_waterbodies": True},
    )
    cfg = {
        "log_parameters": {"showtiming": False, "log_level": "WARNING"},
        "network_topology_parameters": topology,
        "compute_parameters": {
            "parallel_compute_method": method,
            "compute_kernel": "V02-structured",
            "assume_short_ts": True,
            "subnetwork_target_size": 20,
            "cpu_pool": 1,
            "restart_parameters": {"start_datetime": f"{T0:%Y-%m-%d_%H:%M}"},
            "hybrid_parameters": {"run_hybrid_routing": False},
            "forcing_parameters": {
                "qts_subdivisions": 12, "dt": DT, "nts": NTS, "max_loop_size": LOOP_H,
                "qlat_input_folder": str(tmp), "qlat_input_file": "runoff_forcing.nc",
            },
            "data_assimilation_parameters": da_block,
        },
    }
    yfile = tmp / f"run_{len(list(tmp.glob('run_*.yaml')))}.yaml"
    yfile.write_text(yaml.safe_dump(cfg))

    loops: List[pd.DataFrame] = []
    depths: List[pd.DataFrame] = []
    q0s: List[pd.DataFrame] = []
    seen = {}

    def collect(run, run_results, *args, **kwargs):
        ids = np.concatenate([np.asarray(res[0], dtype=np.int64) for res in run_results])
        fvd = np.vstack([res[1] for res in run_results])
        loops.append(pd.DataFrame(fvd[:, 0::3].copy(), index=ids).sort_index())
        depths.append(pd.DataFrame(fvd[:, 2::3].copy(), index=ids).sort_index())

    orig_new_q0 = AbstractNetwork.new_q0
    orig_build = ScalingDA.build.__func__

    def new_q0(self, run_results):
        out = orig_new_q0(self, run_results)
        q0s.append(out.sort_index().copy())
        return out

    def build(cls, network, *a, **k):
        seen["network"] = network
        seen["runtime"] = orig_build(cls, network, *a, **k)
        return seen["runtime"]

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(main, "nwm_output_generator", collect)
        mp.setattr(AbstractNetwork, "new_q0", new_q0)
        mp.setattr(ScalingDA, "build", classmethod(build))
        main.main_v04(["-f", str(yfile)])

    flow = pd.concat(loops, axis=1)
    flow.columns = range(flow.shape[1])
    depth = pd.concat(depths, axis=1)
    depth.columns = range(depth.shape[1])
    return Run(flow=flow, depth=depth, q0=q0s, network=seen["network"], runtime=seen.get("runtime"),
               loops=len(loops))


@dataclass
class Case:
    mode: str
    tmp: Path
    site_seg: Dict[str, int]
    obs: Dict[str, np.ndarray]
    runs: Dict[str, Run] = field(default_factory=dict)

    def get(self, name, **kw):
        if name not in self.runs:
            self.runs[name] = run_troute(self.tmp, self.mode, **kw)
        return self.runs[name]


@pytest.fixture(scope="module", params=["flowpath", "flowline"])
def case(request, tmp_path_factory):
    pytest.importorskip("netCDF4")
    pytest.importorskip("xarray")
    mode = request.param
    tmp = tmp_path_factory.mktemp(f"sda_{mode}")
    write_forcing(FIXTURE, tmp / "runoff_forcing.nc")
    empty = tmp / "ts_empty"
    empty.mkdir()
    # Open loop: nudging on (so the gage crosswalk is built) but no observations.
    ol = run_troute(tmp, mode, timeslices=empty)
    assert ol.loops == NTS // (LOOP_H * 3600 // DT)
    site_seg = {}
    for seg, site in (ol.network.gages or {}).get("gages", {}).items():
        if str(site).isdigit() and int(seg) in ol.flow.index:  # USGS ids, as in real TimeSlices
            site_seg.setdefault(str(site), int(seg))
    obs = {s: FACTOR * ol.flow.loc[g].to_numpy(np.float64) for s, g in site_seg.items()}
    write_timeslices(tmp / "ts_obs", obs)
    c = Case(mode=mode, tmp=tmp, site_seg=site_seg, obs=obs)
    c.runs["open_loop"] = ol
    return c


def nudging(c):
    return c.get("nudging", timeslices=c.tmp / "ts_obs")


def scaling(c, **params):
    return c.get("scaling", timeslices=c.tmp / "ts_obs", scaling=params)


def tree_members(rt):
    return set(rt.trees.seg.tolist())


def tree_gages(rt):
    return set(rt.trees.seg[rt.trees.roots].tolist())


def nested_gage(rt):
    """A source gage that stops another source's tree, and that tree's gage."""
    t = rt.trees
    gage_of_seg = {int(t.seg[r]): t.sites[int(t.tree[r])] for r in t.roots}
    for seg, owner in zip(t.pruned_seg.tolist(), t.pruned_owner.tolist()):
        if seg in gage_of_seg:
            return gage_of_seg[seg], t.sites[int(t.tree[owner])]
    pytest.skip("fixture has no gage directly upstream of another gage's tree")


# ----------------------------------------------------------------- tests


def test_flag_off_is_identical_to_nudging(case):
    off = case.get("block_off", timeslices=case.tmp / "ts_obs", block_off=True)
    nd = nudging(case)
    assert off.runtime is None
    assert np.array_equal(off.flow.to_numpy(), nd.flow.to_numpy())
    assert np.array_equal(off.depth.to_numpy(), nd.depth.to_numpy())


def test_scaling_changes_only_tree_members(case):
    nd, sc = nudging(case), scaling(case)
    rt = sc.runtime
    assert rt is not None and rt.trees.sites, "no trees"
    members = tree_members(rt)
    gages = tree_gages(rt)
    diff = (sc.flow.to_numpy() != nd.flow.to_numpy()).any(axis=1)
    changed = set(sc.flow.index[diff].tolist())
    assert changed, "scaling changed nothing"
    assert changed <= members - gages, f"changed outside trees: {sorted(changed - (members - gages))[:10]}"
    # outside the trees (gages included) every loop matches nudging exactly: the
    # first loop's correction did not enter the second loop's state
    outside = sc.flow.index.difference(sorted(members - gages))
    assert np.array_equal(sc.flow.loc[outside].to_numpy(), nd.flow.loc[outside].to_numpy())
    assert len(rt.loops) == sc.loops and rt.loops[-1]["final"] and not rt.loops[0]["final"]


def test_state_hand_off_only_on_the_final_loop(case):
    nd, sc = nudging(case), scaling(case)
    assert len(nd.q0) == nd.loops and len(sc.q0) == sc.loops + 1
    for k in range(sc.loops):  # q0 taken before the spread equals nudging's, every loop
        pd.testing.assert_frame_equal(sc.q0[k], nd.q0[k])
    final = sc.q0[-1]
    last = sc.flow.columns[-1]
    common = final.index.intersection(sc.flow.index)
    np.testing.assert_array_equal(final.loc[common, "qd0"].to_numpy(), sc.flow.loc[common, last].to_numpy())
    np.testing.assert_array_equal(final.loc[common, "h0"].to_numpy(), sc.depth.loc[common, last].to_numpy())
    assert not final.equals(nd.q0[-1])


def test_holdout_gage_has_no_influence(case):
    sc = scaling(case)
    held, below = nested_gage(sc.runtime)
    hf = case.tmp / "holdout.txt"
    hf.write_text(f"{held}\n")
    perturbed = {s: (v * 2.0 / FACTOR if s == held else v) for s, v in case.obs.items()}
    write_timeslices(case.tmp / "ts_perturbed", perturbed)
    a = case.get("holdout", timeslices=case.tmp / "ts_obs", scaling={"holdout_sites_file": str(hf)})
    b = case.get("holdout_perturbed", timeslices=case.tmp / "ts_perturbed", scaling={"holdout_sites_file": str(hf)})
    assert a.runtime.held_out == {held} and held not in a.runtime.trees.sites
    seg = case.site_seg[held]
    t = a.runtime.trees
    assert seg in t.seg[t.tree == t.sites.index(below)].tolist()  # the gage below now covers it
    assert np.array_equal(a.flow.to_numpy(), b.flow.to_numpy())
    assert np.array_equal(a.depth.to_numpy(), b.depth.to_numpy())
    k = np.arange(900 // DT - 1, NTS, 900 // DT)
    assert not np.allclose(a.flow.loc[seg].to_numpy()[k], case.obs[held][k])  # not nudged
