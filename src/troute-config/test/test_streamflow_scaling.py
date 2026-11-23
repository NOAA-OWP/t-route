"""Configuration of area-scaled streamflow DA (streamflow_da.streamflow_scaling)."""

import pytest
from pydantic import ValidationError

from troute.config.compute_parameters import ComputeParameters, StreamflowDA


def test_off_by_default_with_complete_defaults():
    s = StreamflowDA(streamflow_nudging=True)
    assert s.streamflow_scaling is False
    p = s.streamflow_scaling_parameters
    assert p.theta.default == 0.77 and p.theta.by_vpu == {} and p.theta.per_tree_csv is None
    assert p.max_reach_km == 200.0 and p.holdout_sites_file is None


def test_requires_nudging():
    with pytest.raises(ValidationError, match="requires streamflow_nudging"):
        StreamflowDA(streamflow_scaling=True)
    assert StreamflowDA(streamflow_nudging=True, streamflow_scaling=True).streamflow_scaling


@pytest.mark.parametrize("theta", [0.0, -0.5, float("nan"), float("inf")])
def test_theta_must_be_finite_and_positive(theta):
    with pytest.raises(ValidationError, match="finite and > 0"):
        StreamflowDA(streamflow_nudging=True, streamflow_scaling_parameters={"theta": {"default": theta}})
    with pytest.raises(ValidationError, match="finite and > 0"):
        StreamflowDA(streamflow_nudging=True, streamflow_scaling_parameters={"theta": {"by_vpu": {"01": theta}}})


def test_max_reach_km_must_be_positive():
    with pytest.raises(ValidationError, match="max_reach_km"):
        StreamflowDA(streamflow_nudging=True, streamflow_scaling_parameters={"max_reach_km": 0})


def test_unknown_keys_are_rejected():
    for params in ({"min_flow_cms": 1e-6}, {"innovation_spread_h": 0}, {"theta": {"per_site": {}}}):
        with pytest.raises(ValidationError, match="extra fields not permitted"):
            StreamflowDA(streamflow_nudging=True, streamflow_scaling_parameters=params)


def test_vpu_keys_stay_strings(tmp_path):
    f = tmp_path / "theta.csv"
    f.write_text("gage_id,theta\n01013500,0.8\n")
    p = StreamflowDA(
        streamflow_nudging=True,
        streamflow_scaling_parameters={"theta": {"by_vpu": {"03N": 0.9}, "per_tree_csv": str(f)}},
    ).streamflow_scaling_parameters
    assert p.theta.by_vpu == {"03N": 0.9} and p.theta.per_tree_csv == f


def test_not_with_compute_nodes():
    da = {"streamflow_da": {"streamflow_nudging": True, "streamflow_scaling": True}}
    with pytest.raises(ValidationError, match="compute_nodes"):
        ComputeParameters(compute_nodes=True, data_assimilation_parameters=da)
    assert ComputeParameters(compute_nodes=False, data_assimilation_parameters=da)
    da["streamflow_da"]["streamflow_scaling"] = False
    assert ComputeParameters(compute_nodes=True, data_assimilation_parameters=da)
