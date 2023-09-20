from pydantic import BaseModel, Field, validator
from datetime import datetime

from typing import Optional, List
from typing_extensions import Literal

from .types import FilePath, DirectoryPath
from ._validators import coerce_datetime, coerce_none_to_default


# ---------------------------- Compute Parameters ---------------------------- #

ParallelComputeMethod = Literal[
    "serial",
    "by-network",
    "by-subnetwork-jit",
    "by-subnetwork-jit-clustered",
    "by-subnetwork-diffusive",
    "bmi",
]

ComputeKernel = Literal["V02-structured", "diffusive", "diffusice_cnt"]


class ComputeParameters(BaseModel, extra='forbid'):
    parallel_compute_method: ParallelComputeMethod = "by-network"
    compute_kernel: ComputeKernel = "V02-structured"
    assume_short_ts: bool = False
    subnetwork_target_size: int = 10000
    cpu_pool: Optional[int] = 1
    return_courant: bool = False

    # TODO: default appears to be {}. see nhd_io.read_config_file ~:138
    restart_parameters: "RestartParameters" = Field(default_factory=dict)
    # TODO: default appears to be {}. see nhd_io.read_config_file ~:138
    hybrid_parameters: "HybridParameters" = Field(default_factory=dict)
    # TODO: default appears to be {}. see nhd_io.read_config_file ~:138
    forcing_parameters: "ForcingParameters" = Field(default_factory=dict)
    # TODO: default appears to be {}. see nhd_io.read_config_file ~:138
    data_assimilation_parameters: "DataAssimilationParameters" = Field(
        default_factory=dict
    )


# TODO: determine how to handle context specific required fields
# TODO: consider other ways to handle wrf hydro fields (i.e. subclass)
class RestartParameters(BaseModel, extra='forbid'):
    # NOTE: this is technically optional as it can be derived from the
    # `wrf_hydro_channel_restart_file` if the `start_datetime` is not provided.
    start_datetime: Optional[datetime] = None
    lite_channel_restart_file: Optional[FilePath] = None
    lite_waterbody_restart_file: Optional[FilePath] = None

    wrf_hydro_channel_restart_file: Optional[FilePath] = None
    # NOTE: if `wrf_hydro_channel_restart_file` is given, `wrf_hydro_channel_ID_crosswalk_file` is required
    wrf_hydro_channel_ID_crosswalk_file: Optional[FilePath] = None

    wrf_hydro_channel_ID_crosswalk_file_field_name: Optional[str] = None
    wrf_hydro_channel_restart_upstream_flow_field_name: Optional[str] = None
    wrf_hydro_channel_restart_downstream_flow_field_name: Optional[str] = None
    wrf_hydro_channel_restart_depth_flow_field_name: Optional[str] = None

    wrf_hydro_waterbody_restart_file: Optional[FilePath] = None
    # NOTE: required if `wrf_hydro_waterbody_restart_file`
    wrf_hydro_waterbody_ID_crosswalk_file: Optional[FilePath] = None
    wrf_hydro_waterbody_ID_crosswalk_file_field_name: Optional[str] = None
    wrf_hydro_waterbody_crosswalk_filter_file: Optional[FilePath] = None
    wrf_hydro_waterbody_crosswalk_filter_file_field_name: Optional[str] = None

    # TODO: missing from `v3_doc.yaml`
    # TODO: shorvath: I think we can remove this...
    hyfeature_channel_restart_file: Optional[FilePath] = None

    _coerce_datetime = validator("start_datetime", pre=True, allow_reuse=True)(
        coerce_datetime
    )


# TODO: determine how to handle context specific required fields
class HybridParameters(BaseModel, extra='forbid'):
    # NOTE: required for hybrid simulations
    run_hybrid_routing: bool
    # NOTE: required for hybrid simulations
    diffusive_domain: Optional[FilePath] = None

    use_natl_xsections: bool = False
    # NOTE: required for diffusive routing for natural cross sections
    topobathy_domain: Optional[FilePath] = None

    # TODO: missing from `v3_doc.yaml`
    run_refactored_network: bool = False
    # TODO: missing from `v3_doc.yaml`
    refactored_domain: Optional[FilePath] = None
    # TODO: missing from `v3_doc.yaml`
    refactored_topobathy_domain: Optional[FilePath] = None
    # TODO: missing from `v3_doc.yaml`
    coastal_boundary_domain: Optional[FilePath] = None


class QLateralForcingSet(BaseModel, extra='forbid'):
    nts: "QLateralFiles"


class QLateralFiles(BaseModel, extra='forbid'):
    qlat_files: List[FilePath]


class StreamflowDA(BaseModel, extra='forbid'):
    # NOTE: mandatory for streamflow DA, defaults to False
    streamflow_nudging: bool = False
    # NOTE: mandatory for streamflow DA on NHDNetwork.
    gage_segID_crosswalk_file: Optional[FilePath] = None

    # TODO: not sure if these are dependent on anything
    crosswalk_gage_field: Optional[str] = 'gages'
    crosswalk_segID_field: Optional[str] = 'link'

    # NOTE: required for analysis and
    # TODO: changed the name of this parameter from "wrf_hydro_lastobs_file" to "lastobs_file"
    # Need to update this in t-route as well.
    lastobs_file: Optional[FilePath] = None

    # NOTE: required if lastobs are to be written out during and after simulations
    lastobs_output_folder: Optional[DirectoryPath] = None

    # TODO: missing from `v3_doc.yaml`
    # see troute/DataAssimilation.py :57
    # see troute/nhd_network_utilities_v02.py :765
    diffusive_streamflow_nudging: bool = False


class ReservoirPersistenceDA(BaseModel, extra='forbid'):
    # NOTE: mandatory for USGS reservoir DA, defaults to False
    reservoir_persistence_usgs: bool = False
    # NOTE: mandatory for USACE reservoir DA, defaults to False
    reservoir_persistence_usace: bool = False

    crosswalk_usgs_gage_field: str = "usgs_gage_id"
    crosswalk_usace_gage_field: str = "usace_gage_id"
    crosswalk_usgs_lakeID_field: str = "usgs_lake_id"
    crosswalk_usace_lakeID_field: str = "usace_lake_id"


class ReservoirRfcParameters(BaseModel, extra='forbid'):
    reservoir_rfc_forecasts: bool = False
    reservoir_rfc_forecasts_time_series_path: Optional[FilePath] = None
    reservoir_rfc_forecasts_lookback_hours: int = 28
    reservoir_rfc_forecasts_offset_hours: int = 28
    reservoir_rfc_forecast_persist_days: int = 11


class ReservoirDA(BaseModel, extra='forbid'):
    reservoir_persistence_da: Optional[ReservoirPersistenceDA] = None
    reservoir_rfc_da: Optional[ReservoirRfcParameters] = None
    reservoir_parameter_file: Optional[FilePath] = None


class DataAssimilationParameters(BaseModel, extra='forbid'):
    # NOTE: required for streamflow nudging and/or USGS reservoir DA
    usgs_timeslices_folder: Optional[DirectoryPath] = None
    # NOTE: required for USACE reservoir DA
    usace_timeslices_folder: Optional[DirectoryPath] = None
    # NOTE: required for reservoir DA - suggested value 24 (1 days)
    timeslice_lookback_hours: int = 24

    interpolation_limit_min: int = 59

    wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time: int = 0
    wrf_lastobs_type: str = "obs-based"
    streamflow_da: StreamflowDA
    # NOTE: this appears to be optional. See nwm_routing/input.py ~:439
    reservoir_da: Optional[ReservoirDA] = None

    # NOTE: not present in v3_doc.yaml
    # see troute/nhd_network_utilities_v02.py ~:801
    qc_threshold: float = Field(1, ge=0, le=1)

    _coerce_none_to_default = validator(
        "timeslice_lookback_hours", "qc_threshold", pre=True, allow_reuse=True
    )(coerce_none_to_default)


class ForcingParameters(BaseModel, extra='forbid'):
    qts_subdivisions: int = 12
    dt: int = 300
    # TODO: see note about potentially throwing in v3_doc.yaml
    # aaraney: this is optional if `qlat_forcing_sets` is provided
    qlat_input_folder: Optional[DirectoryPath] = None
    # TODO: mandatory if loop sets will be automatically created
    nts: Optional[int] = 288
    max_loop_size: int = 24
    # NOTE: determine if okay to use this default
    qlat_file_index_col: str = "feature_id"
    qlat_file_value_col: str = "q_lateral"
    qlat_file_gw_bucket_flux_col: str = "qBucket"
    qlat_file_terrain_runoff_col: str = "qSfcLatRunoff"
    qlat_file_pattern_filter: Optional[str] = "*NEXOUT"
    # NOTE:
    # If this variable is not present, make sure nts, qlat_file_pattern_filter, and
    # max_loop_size variables are listed above.
    qlat_forcing_sets: Optional[List[QLateralForcingSet]] = None

    # TODO: shorvath: We might be able to remove binary_nexus_file_folder. 
    # This converts ngen output .csv files into parquet files for t-route.
    binary_nexus_file_folder: Optional[DirectoryPath] = None
    coastal_boundary_input_file: Optional[FilePath] = None
    # NOTE: aaraney: seen as:
    # in code  : "*.NEXOUT", "*NEXOUT*",
    # in config: "*NEXOUT.parquet" "*NEXUS.csv", "nex-*"
    # TODO: shorvath: I belive we no longer use these two arguments...
    # need to double check.
    #nexus_file_pattern_filter: Optional[str] = None
    #nexus_input_folder: Optional[DirectoryPath] = None


ComputeParameters.update_forward_refs()
QLateralForcingSet.update_forward_refs()
ReservoirRfcParameters.update_forward_refs()
