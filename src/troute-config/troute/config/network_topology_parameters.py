from pydantic import BaseModel, Field

from typing import Optional, List, Union
from typing_extensions import Literal

from .types import FilePath, DirectoryPath


class NetworkTopologyParameters(BaseModel, extra='forbid'):
    # TODO: default {}. see nhd_io.read_config_file ~:100
    preprocessing_parameters: "PreprocessingParameters" = Field(default_factory=dict)
    # TODO: not sure if default {}. see nhd_io.read_config_file ~:100
    supernetwork_parameters: "SupernetworkParameters"
    # TODO: default {}. see nhd_io.read_config_file ~:100
    waterbody_parameters: "WaterbodyParameters" = Field(default_factory=dict)

    # TODO: error in v3_doc.yaml; `rfc` is listed as network_topology_parameters parameter.
    # should instead be waterbody_parameters

# TODO: This is an old parameter but probably worth keeping moving forward. However, it is
# not implemented in V4 at the moment (Aug 11, 2023). Need to add this functionality to t-route.
class PreprocessingParameters(BaseModel, extra='forbid'):
    preprocess_only: bool = False
    # NOTE: required if preprocess_only = True
    # TODO: determine if str type
    preprocess_output_folder: Optional[DirectoryPath] = None
    preprocess_output_filename: str = "preprocess_output"
    use_preprocessed_data: bool = False
    # NOTE: required if use_preprocessed_data = True
    # TODO: determine if str type
    preprocess_source_file: Optional[FilePath] = None


class SupernetworkParameters(BaseModel, extra='forbid'):
    title_string: Optional[str] = None
    # TODO: hopefully places in the code can be changed so this is a `Path` instead of a `str`
    geo_file_path: str
    network_type: Literal["HYFeaturesNetwork", "NHDNetwork"] = "HYFeaturesNetwork"
    mask_file_path: Optional[FilePath] = None
    mask_layer_string: str = ""
    # TODO: determine if this is still used
    # TODO: determine what the default for this should be. Not sure if this is right?
    mask_driver_string: Optional[str] = None
    mask_key: int = 0

    columns: Optional["Columns"] = None
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
    duplicate_wb_segments: Optional[List[int]] = Field(
        default_factory=lambda: [
            717696,
            1311881,
            3133581,
            1010832,
            1023120,
            1813525,
            1531545,
            1304859,
            1320604,
            1233435,
            11816,
            1312051,
            2723765,
            2613174,
            846266,
            1304891,
            1233595,
            1996602,
            2822462,
            2384576,
            1021504,
            2360642,
            1326659,
            1826754,
            572364,
            1336910,
            1332558,
            1023054,
            3133527,
            3053788,
            3101661,
            2043487,
            3056866,
            1296744,
            1233515,
            2045165,
            1230577,
            1010164,
            1031669,
            1291638,
            1637751,
        ]
    )
    terminal_code: int = 0
    # TODO: It would be nice if this were a literal / str
    driver_string: Union[str, Literal["NetCDF"]] = "NetCDF"
    layer_string: int = 0


class Columns(BaseModel, extra='forbid'):
    # string, unique segment identifier
    key: str 
    # string, unique identifier of downstream segment
    downstream: str 
    # string, segment length
    dx: str 
    # string, manning's roughness of main channel
    n: str 
    # string, mannings roughness of compound channel
    ncc: str 
    # string, channel slope
    s0: str 
    # string, channel bottom width
    bw: str 
    # string, waterbody identifier
    waterbody: str 
    # string, channel top width
    tw: str 
    # string, compound channel top width
    twcc: str 
    # string, channel bottom altitude
    alt: str 
    # string, muskingum K parameter
    musk: str 
    # string, muskingum X parameter
    musx: str 
    # string, channel sideslope
    cs: str 
    # string, gage ID
    gages: str 


class WaterbodyParameters(BaseModel, extra='forbid'):
    # NOTE: required, True for simulations with waterbodies.
    break_network_at_waterbodies: bool = False
    level_pool: Optional["LevelPool"] = None
    waterbody_null_code: int = -9999


class LevelPool(BaseModel, extra='forbid'):
    # string, filepath to waterbody parameter file (LAKEPARM.nc)
    level_pool_waterbody_parameter_file_path: Optional[FilePath] = None
    level_pool_waterbody_id: Union[str, Literal["lake_id"]] = "lake_id"


NetworkTopologyParameters.update_forward_refs()
PreprocessingParameters.update_forward_refs()
SupernetworkParameters.update_forward_refs()
WaterbodyParameters.update_forward_refs()
LevelPool.update_forward_refs()

