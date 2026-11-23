"""Spread gage corrections over their upstream trees.

Flows are arrays of shape (rows, timesteps). For each node and timestep a
multiplier M is built from the gage up:

    M = 1                                                   at the gage
    M = M_parent * clip((A / A_parent) ** theta, 0, 1)      along a channel
    M = M_parent * clip(Q / max(Q_parent, sum_b Q_b), 0, 1) at a confluence

and the node's correction is M times the gage correction. The branch sum
includes pruned branches. A confluence whose denominator is at or below
MIN_FLOW_CMS, or whose share is not finite, passes on nothing. All arithmetic is
float64.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List, Mapping

import numpy as np

from troute.scaling_da.trees import Trees

LOG = logging.getLogger(__name__)

MIN_FLOW_CMS = 1e-6


@dataclass(frozen=True)
class Plan:
    """Trees mapped to rows of a flow array, with the static factors precomputed."""

    rows: np.ndarray  # flow-array row of each node
    roots: np.ndarray  # node index of each gage
    parent: np.ndarray
    tree: np.ndarray
    channel_factor: np.ndarray  # clipped area factor (1 at gages; unused at confluences)
    junction_nodes: np.ndarray
    junction_parents: np.ndarray
    pruned_rows: np.ndarray
    pruned_owner: np.ndarray
    levels: List[np.ndarray]  # node indices by distance from the gage, from 1 up


def make_plan(trees: Trees, seg_to_row: Mapping[int, int]) -> Plan:
    """Map tree segments to flow-array rows and precompute the area factors."""
    rows = np.asarray([seg_to_row[int(s)] for s in trees.seg], dtype=np.int64)
    pruned_rows = np.asarray([seg_to_row[int(s)] for s in trees.pruned_seg], dtype=np.int64)
    if np.unique(rows).size != rows.size:
        raise ValueError("scaling DA: gage trees share flow rows")

    factor = np.ones(rows.size, dtype=np.float64)
    child = np.flatnonzero(trees.parent >= 0)
    a = trees.area[child]
    a_parent = trees.area[trees.parent[child]]
    has_parent_area = a_parent > 0
    ratio = np.divide(a, a_parent, out=np.zeros_like(a), where=has_parent_area)
    raw = np.zeros(child.size, dtype=np.float64)
    theta = trees.theta[trees.tree[child]]
    for value in np.unique(theta):
        sel = theta == value
        raw[sel] = np.where(has_parent_area[sel], ratio[sel] ** float(value), 0.0)
    factor[child] = np.clip(np.where(np.isfinite(raw), raw, 0.0), 0.0, 1.0)

    junction_nodes = np.flatnonzero(trees.junction)
    order = np.argsort(trees.level, kind="stable")
    depth = trees.level[order]
    levels = [order[depth == d] for d in range(1, int(depth.max(initial=0)) + 1)]
    return Plan(
        rows=rows,
        roots=trees.roots,
        parent=trees.parent,
        tree=trees.tree,
        channel_factor=factor,
        junction_nodes=junction_nodes,
        junction_parents=trees.parent[junction_nodes],
        pruned_rows=pruned_rows,
        pruned_owner=trees.pruned_owner,
        levels=levels,
    )


def _sum_by_owner(values: np.ndarray, owner: np.ndarray, n: int) -> np.ndarray:
    """Sum the rows of values into n rows, by owner index."""
    out = np.zeros((n, values.shape[1]), dtype=np.float64)
    if owner.size:
        order = np.argsort(owner, kind="stable")
        sorted_owner = owner[order]
        starts = np.flatnonzero(np.r_[True, sorted_owner[1:] != sorted_owner[:-1]])
        out[sorted_owner[starts]] = np.add.reduceat(values[order], starts, axis=0)
    return out


def finite_corrections(dq_gage: np.ndarray) -> np.ndarray:
    """dq_gage with non-finite values replaced by 0 (logged)."""
    bad = ~np.isfinite(dq_gage)
    if bad.any():
        LOG.warning("scaling DA: %d non-finite gage correction value(s) treated as 0", int(bad.sum()))
        dq_gage = np.where(bad, 0.0, dq_gage)
    return dq_gage


def corrections(plan: Plan, q: np.ndarray, dq_gage: np.ndarray, gage_flow: np.ndarray) -> np.ndarray:
    """Correction of every tree node, shape (nodes, timesteps), float64.

    Parameters
    ----------
    plan : Plan
    q : ndarray
        Modeled flow, shape (rows, timesteps), read for the confluence split.
    dq_gage : ndarray
        Gage correction per tree, shape (trees, timesteps). Non-finite values count as 0.
    gage_flow : ndarray
        Flow to use at each tree's gage instead of q, shape (trees, timesteps).

    Returns
    -------
    ndarray
        Zero at the gages, which nudging already corrects.
    """
    n, nt = plan.rows.size, q.shape[1]
    dq_gage = finite_corrections(dq_gage)
    if n == 0 or nt == 0:
        return np.zeros((n, nt), dtype=np.float64)

    flow = q[plan.rows].astype(np.float64)
    flow[plan.roots] = gage_flow
    mult = np.repeat(plan.channel_factor[:, None], nt, axis=1)
    jn, jp = plan.junction_nodes, plan.junction_parents
    if jn.size:
        branch_sum = _sum_by_owner(flow[jn], jp, n)
        if plan.pruned_owner.size:
            branch_sum += _sum_by_owner(q[plan.pruned_rows].astype(np.float64), plan.pruned_owner, n)
        denom = np.maximum(flow[jp], branch_sum[jp])
        share = np.divide(flow[jn], denom, out=np.zeros_like(denom), where=denom > MIN_FLOW_CMS)
        mult[jn] = np.clip(np.where(np.isfinite(share), share, 0.0), 0.0, 1.0)
    for nodes in plan.levels:
        mult[nodes] *= mult[plan.parent[nodes]]
    mult *= dq_gage[plan.tree]
    mult[plan.roots] = 0.0
    return mult
