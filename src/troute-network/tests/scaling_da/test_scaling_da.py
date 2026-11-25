"""Unit tests for area-scaled streamflow DA on a small hand-checkable network.

    40 (lake) -> 21 -+
                     +-> 20 -+
                 22 -+       +-> 11 -> 10 (gage A)
             30 (gage B) ----+

theta = 0.5. Areas: 10: 100, 11: 81, 20: 40, 21: 20, 22: 15, 30: 30, 40: 5.
Lengths are 10 km each. Flows: 11: 10, 20: 6, 21: 3, 22: 2, 30: 4.
"""

import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from troute.scaling_da.runtime import ScalingDA, read_theta_csv
from troute.scaling_da.spread import corrections, make_plan
from troute.scaling_da.trees import build_trees

RCONN = {10: [11], 11: [20, 30], 20: [21, 22], 21: [40], 22: [], 30: [], 40: []}
AREA = pd.Series({10: 100.0, 11: 81.0, 20: 40.0, 21: 20.0, 22: 15.0, 30: 30.0, 40: 5.0})
LENGTH = pd.Series(10.0, index=AREA.index)
FLOW = {10: 12.0, 11: 10.0, 20: 6.0, 21: 3.0, 22: 2.0, 30: 4.0, 40: 1.0}
ROW = {s: i for i, s in enumerate(sorted(FLOW))}
STOPS = {10, 30, 40}


def tree_of_a(max_reach_km=200.0, area=AREA):
    trees, dropped = build_trees(RCONN, {"A": 10}, area, LENGTH, STOPS, {"A": 0.5}, max_reach_km)
    return trees, dropped


def spread(trees, flow=FLOW, dq=1.0, nt=1):
    q = np.zeros((len(ROW), nt), dtype=np.float32)
    for s, v in flow.items():
        q[ROW[s]] = v
    plan = make_plan(trees, ROW)
    dq_gage = np.full((len(trees.sites), nt), dq)
    gage_flow = q[plan.rows[plan.roots]].astype(np.float64)
    inc = corrections(plan, q, dq_gage, gage_flow)
    return {int(s): inc[i] for i, s in enumerate(trees.seg)}


def test_tree_stops_at_other_gages_and_lakes():
    trees, dropped = tree_of_a()
    assert not dropped
    assert trees.seg.tolist() == [10, 11, 20, 21, 22]
    assert sorted(trees.pruned_seg.tolist()) == [30, 40]
    assert trees.junction.tolist() == [False, False, True, True, True]


def test_linear_step_scales_by_area():
    inc = spread(tree_of_a()[0])
    assert inc[10][0] == 0.0  # the gage itself is left to nudging
    assert inc[11][0] == pytest.approx(math.sqrt(81 / 100))


def test_confluence_splits_by_flow_and_counts_stopped_branches():
    inc = spread(tree_of_a()[0])
    # 20 gets 6 / max(10, 6 + 4): gage B's flow counts but B gets nothing
    assert inc[20][0] == pytest.approx(0.9 * 0.6)
    assert inc[21][0] == pytest.approx(0.9 * 0.6 * 3 / 6)
    assert inc[22][0] == pytest.approx(0.9 * 0.6 * 2 / 6)


def test_parent_flow_bounds_the_split():
    flow = {**FLOW, 11: 20.0}
    inc = spread(tree_of_a()[0], flow=flow)
    assert inc[20][0] == pytest.approx(0.9 * 6 / 20)


def test_distance_cap_drops_far_nodes_but_keeps_their_flow_in_the_split():
    trees, _ = tree_of_a(max_reach_km=15.0)  # 11 is 10 km above the gage outlet, 20 is 20 km
    assert trees.seg.tolist() == [10, 11]
    assert sorted(trees.pruned_seg.tolist()) == [20, 30]


def test_unusable_length_puts_everything_above_out_of_reach():
    length = LENGTH.copy()
    length[11] = np.nan
    trees, _ = build_trees(RCONN, {"A": 10}, AREA, length, STOPS, {"A": 0.5}, 1e9)
    assert trees.seg.tolist() == [10, 11]


def test_area_increasing_upstream_is_clamped():
    area = AREA.copy()
    area[11] = 150.0
    inc = spread(tree_of_a(area=area)[0])
    assert inc[11][0] == 1.0


def test_tiny_or_missing_flow_passes_nothing_at_a_confluence():
    inc = spread(tree_of_a()[0], flow={**FLOW, 11: 0.0, 20: 0.0, 30: 0.0})
    assert inc[20][0] == 0.0 and inc[21][0] == 0.0
    inc = spread(tree_of_a()[0], flow={**FLOW, 21: np.nan})
    assert inc[21][0] == 0.0 and np.isfinite(inc[20][0])


def test_non_finite_gage_correction_counts_as_zero():
    inc = spread(tree_of_a()[0], dq=np.inf)
    assert all(v[0] == 0.0 for v in inc.values())


def test_gage_without_drainage_area_gets_no_tree():
    area = AREA.copy()
    area[22] = np.nan
    trees, dropped = tree_of_a(area=area)
    assert dropped == {"A": "missing drainage area"} and trees.seg.size == 0


def test_network_that_is_not_a_tree_raises():
    rconn = {**RCONN, 22: [21]}
    with pytest.raises(ValueError, match="reached twice"):
        build_trees(rconn, {"A": 10}, AREA, LENGTH, STOPS, {"A": 0.5}, 200.0)


# ------------------------------------------------------------ ScalingDA

def fake_network():
    attrs = pd.DataFrame({"total_dasqkm": AREA, "vpuid": "01"})
    return SimpleNamespace(
        reverse_network=RCONN, segment_attributes=attrs,
        dataframe=pd.DataFrame({"dx": LENGTH * 1000.0}),
        gages={"gages": {10: "01000010", 30: "01000030"}},
        waterbody_connections={40: 40}, diffusive_network_data={},
    )


def da_config(**params):
    params.setdefault("theta", {"default": 0.5})
    return {"streamflow_da": {"streamflow_nudging": True, "streamflow_scaling": True,
                              "streamflow_scaling_parameters": params}}


def run_results(flow=FLOW, nudge_a=2.0, depth=1.0, nts=2):
    ids = np.asarray(sorted(flow), dtype=np.int64)
    fvd = np.zeros((ids.size, 3 * nts), dtype=np.float32)
    fvd[:, 0::3] = np.asarray([flow[s] for s in ids], dtype=np.float32)[:, None]
    fvd[:, 2::3] = depth
    nudge = np.zeros((2, nts + 1), dtype=np.float32)
    nudge[0, 1:] = nudge_a
    fvd[ids == 10, 0::3] += nudge_a  # published gage flow includes the nudge
    return [[ids, fvd, None, (np.array([10, 30]),), None, None, None, None, nudge]]


def flows(results):
    r = results[0]
    return {int(s): r[1][i, 0::3].copy() for i, s in enumerate(r[0])}


def test_off_by_default():
    cfg = da_config()
    cfg["streamflow_da"]["streamflow_scaling"] = False
    assert ScalingDA.build(fake_network(), cfg, observed={"01000010"}) is None


def test_apply_corrects_only_tree_members_and_depth_follows():
    rt = ScalingDA.build(fake_network(), da_config(), observed={"01000010", "01000030"})
    before = run_results()
    res = run_results()
    assert rt.apply(res, dt=300.0, final=False) is False
    q0, q1 = flows(before), flows(res)
    for s in (10, 30, 40):
        np.testing.assert_array_equal(q1[s], q0[s])
    np.testing.assert_allclose(q1[11], 10.0 + 2.0 * 0.9, rtol=1e-6)
    np.testing.assert_allclose(q1[20], 6.0 + 2.0 * 0.9 * 0.6, rtol=1e-6)
    h = res[0][1][list(res[0][0]).index(11), 2::3]
    np.testing.assert_allclose(h, (11.8 / 10.0) ** 0.6, rtol=1e-6)
    assert rt.loops[-1]["segments"] == 4
    assert rt.loops[-1]["volume_m3"] == pytest.approx(2 * 2.0 * 300.0 * sum(
        [0.9, 0.9 * 0.6, 0.9 * 0.3, 0.9 * 0.2]), rel=1e-5)
    assert rt.apply(run_results(), dt=300.0, final=True) is True


def test_corrected_flow_is_floored_at_zero_and_depth_ratio_capped():
    rt = ScalingDA.build(fake_network(), da_config(), observed={"01000010", "01000030"})
    res = run_results(nudge_a=-50.0)
    rt.apply(res, dt=300.0, final=False)
    q = flows(res)
    assert (q[11] == 0.0).all() and (q[21] == 0.0).all()
    res = run_results(nudge_a=40.0)
    rt.apply(res, dt=300.0, final=False)
    h = res[0][1][list(res[0][0]).index(11), 2::3]
    np.testing.assert_allclose(h, 2.0 ** 0.6, rtol=1e-6)


def test_gage_without_observations_is_not_a_stop():
    rt = ScalingDA.build(fake_network(), da_config(), observed={"01000010"})
    assert rt.no_observations == {"01000030"}
    assert 30 in rt.trees.seg.tolist()


def test_holdout_removes_the_gage_from_nudging_and_trees(tmp_path):
    f = tmp_path / "holdout.txt"
    f.write_text("# withheld\n01000030\n")
    rt = ScalingDA.build(fake_network(), da_config(holdout_sites_file=str(f)), observed={"01000010", "01000030"})
    assert rt.held_out == {"01000030"} and rt.trees.sites == ["01000010"]
    assert 30 in rt.trees.seg.tolist()
    usgs = pd.DataFrame({"t": [1.0, 2.0]}, index=[10, 30])
    kept, _ = rt.nudging_inputs(usgs, usgs)
    assert kept.index.tolist() == [10]


def test_unknown_holdout_id_raises(tmp_path):
    f = tmp_path / "holdout.txt"
    f.write_text("0999\n")
    with pytest.raises(ValueError, match="not gages of this domain"):
        ScalingDA.build(fake_network(), da_config(holdout_sites_file=str(f)), observed={"01000010"})


def test_theta_csv_wins_over_vpu_and_default(tmp_path):
    f = tmp_path / "theta.csv"
    f.write_text("gage_id,theta\n01000010,1.0\n")
    assert read_theta_csv(f) == {"01000010": 1.0}
    obs = {"01000010", "01000030"}
    rt = ScalingDA.build(fake_network(), da_config(theta={"by_vpu": {"01": 0.9}}), observed=obs)
    assert rt.trees.theta.tolist() == [0.9, 0.9]
    rt = ScalingDA.build(fake_network(), da_config(theta={"by_vpu": {"01": 0.9}, "per_tree_csv": str(f)}),
                         observed=obs)
    assert rt.trees.theta.tolist() == [1.0, 0.9]
    with pytest.raises(ValueError, match="none of its 1 key"):
        ScalingDA.build(fake_network(), da_config(theta={"by_vpu": {"1": 0.9}}), observed=obs)


def test_flow_is_rewritten_only_in_blocks_with_a_correction():
    rt = ScalingDA.build(fake_network(), da_config(), observed={"01000010", "01000030"})
    res = run_results(nudge_a=0.0, nts=72)
    ids, fvd, nudge = list(res[0][0]), res[0][1], res[0][8]
    nudge[0, 1:49] = 2.0  # gage correction in the first two 24-step blocks only
    fvd[ids.index(10), 0:144:3] += 2.0
    fvd[ids.index(21), [30 * 3, 60 * 3]] = -0.5  # negative routed flow (no share at those steps)
    rt.apply(res, dt=300.0, final=False)
    q21 = flows(res)[21]
    assert q21[30] == 0.0  # in a corrected block: floored at 0
    assert q21[60] == np.float32(-0.5)  # in a block without correction: left as routed
    assert q21[0] > 3.0 and q21[50] == 3.0
