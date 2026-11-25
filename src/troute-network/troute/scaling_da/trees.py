"""Upstream trees of the source gages, flattened into arrays."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Dict, List, Mapping, Sequence, Set, Tuple

import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class Trees:
    """All gage trees as flat node arrays.

    Trees are stored one after another, each in breadth-first order with its
    gage first, so a node's parent always comes before it.

    Attributes
    ----------
    sites : list of str
        Gage id of each tree.
    seg : ndarray of int64
        Segment id of each node.
    tree : ndarray of int64
        Tree index of each node.
    parent : ndarray of int64
        Node index of the downstream neighbour; -1 at each gage.
    level : ndarray of int64
        Number of steps from the gage (0 at the gage).
    junction : ndarray of bool
        True where the parent receives flow from two or more segments in the
        full network, so the flow split applies instead of the area ratio.
    area : ndarray of float64
        Total drainage area (km2) of each node.
    theta : ndarray of float64
        Area-scaling exponent of each tree.
    pruned_seg, pruned_owner : ndarray of int64
        Upstream neighbours that are not in the tree (stopped or beyond the
        cap), and the node each one flows into. They get no correction, but
        their flow counts in their confluence's denominator.
    """

    sites: List[str]
    seg: np.ndarray
    tree: np.ndarray
    parent: np.ndarray
    level: np.ndarray
    junction: np.ndarray
    area: np.ndarray
    theta: np.ndarray
    pruned_seg: np.ndarray
    pruned_owner: np.ndarray

    @property
    def roots(self) -> np.ndarray:
        """Node index of each tree's gage."""
        return np.flatnonzero(self.parent < 0)


def _trace(gage_seg: int, rconn: Mapping[int, Sequence[int]], stops: Set[int]):
    """Breadth-first walk upstream from gage_seg, ending each branch before a stop."""
    segs, parent, pruned = [gage_seg], [-1], [[]]
    seen = {gage_seg}
    i = 0
    while i < len(segs):
        for up in rconn.get(segs[i], ()):
            up = int(up)
            if up in stops:
                pruned[i].append(up)
            elif up in seen:
                raise ValueError(f"segment {up} is reached twice upstream of segment {gage_seg}")
            else:
                seen.add(up)
                segs.append(up)
                parent.append(i)
                pruned.append([])
        i += 1
    return segs, parent, pruned


def _cap(segs, parent, pruned, length_km: pd.Series, max_reach_km: float):
    """Drop nodes farther than max_reach_km from the gage.

    A node's distance is the summed length of every segment below it, down to
    and including the gage segment. An unusable length puts everything above it
    out of reach. A dropped node whose parent is kept becomes a pruned
    neighbour of that parent.
    """
    lengths = length_km.reindex(segs).to_numpy(dtype=np.float64)
    dist = np.zeros(len(segs), dtype=np.float64)
    for i in range(1, len(segs)):
        p = parent[i]
        step = float(lengths[p])
        dist[i] = dist[p] + step if (math.isfinite(step) and step > 0) else math.inf
    keep = dist <= max_reach_km
    if keep.all():
        return segs, parent, pruned
    for i in np.flatnonzero(~keep):
        if keep[parent[i]]:
            pruned[parent[i]].append(segs[i])
    new_index = np.cumsum(keep) - 1
    kept = np.flatnonzero(keep)
    return (
        [segs[i] for i in kept],
        [-1] + [int(new_index[parent[i]]) for i in kept[1:]],
        [pruned[i] for i in kept],
    )


def build_trees(
    rconn: Mapping[int, Sequence[int]],
    gages: Dict[str, int],
    area_km2: pd.Series,
    length_km: pd.Series,
    stops: Set[int],
    theta: Mapping[str, float],
    max_reach_km: float,
) -> Tuple[Trees, Dict[str, str]]:
    """Build one upstream tree per gage.

    Parameters
    ----------
    rconn : mapping
        Segment id to the ids of the segments flowing into it.
    gages : mapping
        Gage id to segment id, in the order the trees are stored.
    area_km2, length_km : pandas.Series
        Total drainage area and length per segment id.
    stops : set of int
        Segments no tree enters: other gages, lakes, the diffusive domain.
    theta : mapping
        Gage id to exponent.
    max_reach_km : float
        Distance cap.

    Returns
    -------
    trees : Trees
    dropped : dict
        Gage id to the reason it has no tree.
    """
    dropped: Dict[str, str] = {}
    cols: Dict[str, list] = {k: [] for k in ("seg", "tree", "parent", "level", "junction", "area",
                                              "pruned_seg", "pruned_owner")}
    sites: List[str] = []
    offset = 0
    for site in gages:
        gage_seg = gages[site]
        a_gage = area_km2.get(gage_seg)
        if a_gage is None or not (math.isfinite(a_gage) and a_gage > 0):
            dropped[site] = "missing drainage area"
            continue
        segs, parent, pruned = _trace(gage_seg, rconn, stops - {gage_seg})
        segs, parent, pruned = _cap(segs, parent, pruned, length_km, max_reach_km)
        area = area_km2.reindex(segs).to_numpy(dtype=np.float64)
        if np.isnan(area).any() or np.isinf(area).any():
            dropped[site] = "missing drainage area"
            continue
        n = len(segs)
        level = np.zeros(n, dtype=np.int64)
        junction = np.zeros(n, dtype=bool)
        for i in range(1, n):
            level[i] = level[parent[i]] + 1
            junction[i] = len(rconn.get(segs[parent[i]], ())) >= 2
        cols["seg"].append(np.asarray(segs, dtype=np.int64))
        cols["tree"].append(np.full(n, len(sites), dtype=np.int64))
        cols["parent"].append(np.asarray(parent, dtype=np.int64) + np.where(np.arange(n) > 0, offset, 0))
        cols["level"].append(level)
        cols["junction"].append(junction)
        cols["area"].append(area)
        cols["pruned_seg"].append(np.asarray([s for p in pruned for s in p], dtype=np.int64))
        cols["pruned_owner"].append(np.repeat(np.arange(n, dtype=np.int64) + offset, [len(p) for p in pruned]))
        sites.append(site)
        offset += n

    dtypes = {"junction": bool, "area": np.float64}
    flat = {k: (np.concatenate(v) if v else np.empty(0, dtypes.get(k, np.int64))) for k, v in cols.items()}
    if np.unique(flat["seg"]).size != flat["seg"].size:
        raise RuntimeError("scaling DA: gage trees overlap")
    trees = Trees(sites=sites, theta=np.asarray([float(theta[s]) for s in sites], dtype=np.float64), **flat)
    return trees, dropped
