from pydantic import BaseModel, Field, validator

from typing import Optional, List, Union, Dict, Any
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

# TODO: This is an old parameter but probably worth keeping moving forward. However, it is
# not implemented in V4 at the moment (Aug 11, 2023). Need to add this functionality to t-route.
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
    geo_file_path: FilePath
    network_type: Literal["HYFeaturesNetwork", "NHDNetwork"] = "HYFeaturesNetwork"
    flowpath_edge_list: Optional[str] = None
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

    terminal_code: int = 0
    # TODO: It would be nice if this were a literal / str
    driver_string: Union[str, Literal["NetCDF"]] = "NetCDF"
    layer_string: int = 0
    
    @validator("columns", always=True)
    def get_columns(cls, columns: dict, values: Dict[str, Any]) -> dict:
        if columns is None:
            if values['network_type']=="HYFeaturesNetwork":
                default_columns = {
                    'key'       : 'id',
                    'downstream': 'toid',
                    'dx'        : 'length_m',
                    'n'         : 'n',
                    'ncc'       : 'nCC',
                    's0'        : 'So',
                    'bw'        : 'BtmWdth',
                    'waterbody' : 'rl_NHDWaterbodyComID',
                    'gages'     : 'rl_gages',
                    'tw'        : 'TopWdth',
                    'twcc'      : 'TopWdthCC',
                    'musk'      : 'MusK',
                    'musx'      : 'MusX',
                    'cs'        : 'ChSlp',
                    'alt'       : 'alt',
                    'mainstem'  : 'mainstem',
                    }
            else:
                default_columns = {
                    'key'       : 'link',
                    'downstream': 'to',
                    'dx'        : 'Length',
                    'n'         : 'n',
                    'ncc'       : 'nCC',
                    's0'        : 'So',
                    'bw'        : 'BtmWdth',
                    'waterbody' : 'NHDWaterbodyComID',
                    'gages'     : 'gages',
                    'tw'        : 'TopWdth',
                    'twcc'      : 'TopWdthCC',
                    'alt'       : 'alt',
                    'musk'      : 'MusK',
                    'musx'      : 'MusX',
                    'cs'        : 'ChSlp',
                    }
        else:
            default_columns = columns
        return default_columns


class Columns(BaseModel):
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
    waterbody: Optional[str] 
    # string, channel top width
    tw: str 
    # string, compound channel top width
    twcc: str 
    # string, channel bottom altitude
    alt: Optional[str] 
    # string, muskingum K parameter
    musk: str 
    # string, muskingum X parameter
    musx: str 
    # string, channel sideslope
    cs: str 
    # string, gage ID
    gages: Optional[str]
    # string, mainstem ID
    mainstem: Optional[str]


class WaterbodyParameters(BaseModel):
    # NOTE: required, True for simulations with waterbodies.
    break_network_at_waterbodies: bool = False
    level_pool: Optional["LevelPool"] = None
    waterbody_null_code: int = -9999


class LevelPool(BaseModel):
    # string, filepath to waterbody parameter file (LAKEPARM.nc)
    level_pool_waterbody_parameter_file_path: Optional[FilePath] = None
    level_pool_waterbody_id: Union[str, Literal["lake_id"]] = "lake_id"


NetworkTopologyParameters.update_forward_refs()
PreprocessingParameters.update_forward_refs()
SupernetworkParameters.update_forward_refs()
WaterbodyParameters.update_forward_refs()
LevelPool.update_forward_refs()

