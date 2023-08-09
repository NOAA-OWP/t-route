from pydantic import BaseModel, Field

from typing import Optional, List, Union
from typing_extensions import Literal

from .types import FilePath, DirectoryPath


class NetworkTopologyParameters(BaseModel):
    # TODO: default {}. see nhd_io.read_config_file ~:100
    preprocessing_parameters: "PreprocessingParameters" = Field(default_factory=dict)
    # TODO: not sure if default {}. see nhd_io.read_config_file ~:100
    supernetwork_parameters: "SupernetworkParameters"
    # TODO: default {}. see nhd_io.read_config_file ~:100
    waterbody_parameters: "WaterbodyParameters" = Field(default_factory=dict)

    # TODO: error in v3_doc.yaml; `rfc` is listed as network_topology_parameters parameter.
    # should instead be waterbody_parameters


class PreprocessingParameters(BaseModel):
    preprocess_only: bool = False
    # NOTE: required if preprocess_only = True
    # TODO: determine if str type
    preprocess_output_folder: Optional[DirectoryPath] = None
    preprocess_output_filename: str = "preprocess_output"
    use_preprocessed_data: bool = False
    # NOTE: required if use_preprocessed_data = True
    # TODO: determine if str type
    preprocess_source_file: Optional[FilePath] = None


class SupernetworkParameters(BaseModel):
    title_string: Optional[str] = None
    # TODO: hopefully places in the code can be changed so this is a `Path` instead of a `str`
    geo_file_path: str
    # TODO: not sure if this is correct
    geo_file_type: Literal["HYFeaturesNetwork", "NHDNetwork"] = "HYFeaturesNetwork"
    mask_file_path: Optional[FilePath] = None
    mask_layer_string: str = ""
    # TODO: determine if this is still used
    # TODO: determine what the default for this should be. Not sure if this is right?
    mask_driver_string: Optional[str] = None
    mask_key: int = 0

    # TODO: Not sure if this should default to None
    columns: Optional["Columns"] = None
    waterbody_null_code: int = -9999
    # NOTE: required for CONUS-scale simulations with NWM 2.1 or 3.0 Route_Link.nc data
    synthetic_wb_segments: Optional[List[int]] = Field(
        default_factory=lambda: [
            4800002,
            4800004,
            4800006,
            4800007,
        ]
    )
    synthetic_wb_id_offset: float = 9.99e11
    terminal_code: int = 0
    # TODO: It would be nice if this were a literal / str
    driver_string: Union[str, Literal["NetCDF"]] = "NetCDF"
    layer_string: int = 0
    # TODO: this parameter is in `read_config_file` fn, but ive not seen it in used in a config
    # file
    # flowpath_edge_list

    # TODO: missing from `v3_doc.yaml`
    ngen_nexus_file: Optional[FilePath] = None


# TODO: is it okay to set defaults for these?
class Columns(BaseModel):
    # string, unique segment identifier
    key: str = "id"
    # string, unique identifier of downstream segment
    downstream: str = "toid"
    # string, segment length
    dx: str = "lengthkm"
    # string, manning's roughness of main channel
    n: str = "n"
    # string, mannings roughness of compound channel
    ncc: str = "nCC"
    # string, channel slope
    s0: str = "So"
    # string, channel bottom width
    bw: str = "BtmWdth"
    # string, waterbody identifier
    waterbody: str = "rl_NHDWaterbodyComID"
    # string, channel top width
    tw: str = "TopWdth"
    # string, compound channel top width
    twcc: str = "TopWdthCC"
    # string, channel bottom altitude
    alt: str = "alt"
    # string, muskingum K parameter
    musk: str = "MusK"
    # string, muskingum X parameter
    musx: str = "MusX"
    # string, channel sideslope
    cs: str = "ChSlp"
    # string, gage ID
    gages: str = "rl_gages"


class WaterbodyParameters(BaseModel):
    # NOTE: required, True for simulations with waterbodies.
    break_network_at_waterbodies: bool = False
    level_pool: Optional["LevelPool"] = None
    rfc: Optional["RfcParameters"] = None


# TODO: not sure if this is still used
# TODO: not sure if it is okay to use these defaults
class LevelPool(BaseModel):
    # string, filepath to waterbody parameter file (LAKEPARM.nc)
    level_pool_waterbody_parameter_file_path: Optional[FilePath] = None
    level_pool_waterbody_id: Union[str, Literal["lake_id"]] = "lake_id"
    # NOTE: not sure if the below fields are still used by t-route
    level_pool_waterbody_area: str = "LkArea"
    level_pool_weir_elevation: str = "WeirE"
    level_pool_waterbody_max_elevation: str = "LkMxE"
    level_pool_outfall_weir_coefficient: str = "WeirC"
    level_pool_outfall_weir_length: str = "WeirL"
    level_pool_overall_dam_length: str = "DamL"
    level_pool_orifice_elevation: str = "OrificeE"
    level_pool_orifice_coefficient: str = "OrificeC"
    level_pool_orifice_area: str = "OrificeA"

    # TODO: missing from `v3_doc.yaml`
    # NOTE: not sure if this is required
    reservoir_parameter_file: Optional[FilePath] = None


class RfcParameters(BaseModel):
    # NOTE: required for RFC forecasting
    reservoir_parameter_file: FilePath
    reservoir_rfc_forecasts: bool = False
    # NOTE: required if reservoir_rfc_forecasts = True
    reservoir_rfc_forecasts_time_series_path: Optional[FilePath] = None
    reservoir_rfc_forecasts_lookback_hours: int


NetworkTopologyParameters.update_forward_refs()
PreprocessingParameters.update_forward_refs()
SupernetworkParameters.update_forward_refs()
WaterbodyParameters.update_forward_refs()
LevelPool.update_forward_refs()
RfcParameters.update_forward_refs()
