from pathlib import Path
from pydantic import BaseModel, Field, validator

from typing import Optional, List
from typing_extensions import Annotated, Literal
from .types import FilePath, DirectoryPath

streamOutput_allowedTypes = Literal['.csv', '.nc', '.pkl']

class OutputParameters(BaseModel):
    chanobs_output: Optional["ChanobsOutput"] = None
    # NOTE: this appears to be optional. See nwm_routing/input.py ~:477
    csv_output: Optional["CsvOutput"] = None
    # NOTE: this appears to be optional. See nwm_routing/input.py ~:496
    chrtout_output: Optional["ChrtoutOutput"] = None
    lite_restart: Optional["LiteRestart"] = None
    # NOTE: this appears to be optional. See nwm_routing/input.py ~:520
    hydro_rst_output: Optional["HydroRstOutput"] = None
    # TODO: default appears to be {}. see nhd_io.read_config_file ~:141
    # shorvath: parity_parameters defaults to {}, but omitting 'wrf_hydro_parity_check'
    # from output_parameters will successfully skip lines~112-115 in __main__.py if this
    # parameter is left blank.
    wrf_hydro_parity_check: Optional["WrfHydroParityCheck"] = None
    # NOTE: mandatory if writing results to lakeout.
    lakeout_output: Optional[DirectoryPath] = None

    # NOTE: assuming this should be removed
    # TODO: missing from `v3_doc.yaml`
    # see nwm_routing/output.py :114
    test_output: Optional[FilePath] = None
    stream_output: Optional["StreamOutput"] = None
    # NOTE: mandatory if writing results to lastobs
    lastobs_output: Optional[DirectoryPath] = None

class ChanobsOutput(BaseModel):
    # NOTE: required if writing chanobs files
    chanobs_output_directory: Optional[DirectoryPath] = None
    # NOTE: required if writing chanobs files
    # NOTE: is `Path` b.c. is output file
    chanobs_filepath: Optional[Path] = None


class CsvOutput(BaseModel):
    # NOTE: required if writing results to csv
    csv_output_folder: Optional[DirectoryPath] = None
    csv_output_segments: Optional[List[str]] = None


class ChrtoutOutput(BaseModel):
    # NOTE: mandatory if writing results to CHRTOUT.
    wrf_hydro_channel_output_source_folder: Optional[DirectoryPath] = None


class LiteRestart(BaseModel):
    # NOTE: required if writing restart data lite files.
    lite_restart_output_directory: Optional[DirectoryPath] = None


class HydroRstOutput(BaseModel):
    # NOTE: required if writing restart data to HYDRO_RST
    wrf_hydro_restart_dir: Optional[DirectoryPath] = None
    wrf_hydro_channel_restart_pattern_filter: str = "HYDRO_RST.*"

    wrf_hydro_channel_restart_source_directory: Optional[DirectoryPath] = None
    wrf_hydro_channel_output_source_folder: Optional[DirectoryPath] = None


class WrfHydroParityCheck(BaseModel):
    # NOTE: required for parity check to occur
    # TODO: not sure if this should be optional?
    # shorvath: I'm ok with removing parity_checks for t-routeV4...
    parity_check_input_folder: Optional[DirectoryPath] = None
    parity_check_file_index_col: str
    parity_check_file_value_col: str
    parity_check_compare_node: str
    parity_check_compare_file_sets: Optional[List["ParityCheckCompareFileSet"]] = None


class ParityCheckCompareFileSet(BaseModel):
    validation_files: List[FilePath]

class StreamOutput(BaseModel):
    # NOTE: required if writing StreamOutput files
    stream_output_directory: Optional[DirectoryPath] = None
    stream_output_time: int = 1
    stream_output_type:streamOutput_allowedTypes = ".nc"
    stream_output_internal_frequency: Annotated[int, Field(strict=True, ge=5)] = 5
    @validator('stream_output_internal_frequency')
    def validate_stream_output_internal_frequency(cls, value, values):
        if value is not None:
            if value % 5 != 0:
                raise ValueError("stream_output_internal_frequency must be a multiple of 5.")
            if values.get('stream_output_time') != -1 and value / 60 > values['stream_output_time']:
                raise ValueError("stream_output_internal_frequency should be less than or equal to stream_output_time in minutes.")
        return value
 

OutputParameters.update_forward_refs()
WrfHydroParityCheck.update_forward_refs()
