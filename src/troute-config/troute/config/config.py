import os

from pydantic import Field, BaseModel, root_validator, validator
from pathlib import Path
from .types import FilePath, DirectoryPath

from typing import Any, Dict, Optional
from typing_extensions import Self

from .logging_parameters import LoggingParameters
from .network_topology_parameters import NetworkTopologyParameters
from .compute_parameters import ComputeParameters
from .output_parameters import OutputParameters
from .bmi_parameters import BMIParameters
from ._utils import use_strict

def classCheck(var):
    return str(type(var)).startswith("<class 'troute")

def varIsPath(var):
    return ('/' in str(var) or 'Path' in str(type(var)))

class Config(BaseModel, extra='forbid'):
    log_parameters: LoggingParameters = Field(default_factory=LoggingParameters)
    # TODO: not sure if default is None or {}. see nhd_io.read_config_file ~:100
    network_topology_parameters: Optional[NetworkTopologyParameters] = None
    # TODO: default appears to be {}. see nhd_io.read_config_file ~:138
    compute_parameters: ComputeParameters = Field(default_factory=dict)
    # TODO: default appears to be {}. see nhd_io.read_config_file ~:141
    output_parameters: OutputParameters = Field(default_factory=dict)
    bmi_parameters: Optional[BMIParameters] = None

    # TODO: `reservoir_data_assimilation_parameters` missing from `v3_doc.yaml` and only included
    # test/BMI/bmi_reservoir_example.yaml

    @classmethod
    def with_strict_mode(cls, **data: Dict[str, Any]) -> Self:
        """
        Create a `Config` instance that verifies existence of file and directory field types. If any
        do not exist, a `pydantic.ValidationError` is raised.

        Returns
        -------
        Self
            `Config` instance
        """
        with use_strict():
            return cls(**data)
    

    # verifies all fields that are either "Path" types, or those that contain a slash ('/') in the name as valid paths
    @validator("*") # validates all fields
    def validatePathGlobal(cls, inputConfig):
        
        vars1 = vars(inputConfig)
        vars1Values = [i for i in vars1.values()]
        
        for var1 in vars1Values:

            if (classCheck(var1)):
                vars2 = vars(var1)
                vars2Values = [i for i in vars2.values()]

                for var2 in vars2Values:

                    if (classCheck(var2)):
                        vars3 = vars(var2)
                        vars3Values = [i for i in vars3.values()]

                        for var3 in vars3Values:

                            if (classCheck(var3)):
                                vars4 = vars(var3)
                                vars4Values = [i for i in vars4.values()]

                                for var4 in vars4Values:

                                    if (varIsPath(var4)):
                                        assert os.path.exists(var4), 'Path does not exist: '+str(var4)  

                            elif (varIsPath(var3)):
                                assert os.path.exists(var3), 'Path does not exist: '+str(var3)  

                    elif (varIsPath(var2)):
                        assert os.path.exists(var2), 'Path does not exist: '+str(var2)  

            elif (varIsPath(var1)):
                assert os.path.exists(var1), 'Path does not exist: '+str(var1)  

        return inputConfig


    ##########################################################################################
    # Validate that, given certain configuration inputs, other inputs are provided.
    ##########################################################################################

    @root_validator(skip_on_failure=True)
    def check_levelpool_filepath(cls, values):
        network_type = values['network_topology_parameters'].supernetwork_parameters.network_type
        waterbody_parameters = values['network_topology_parameters'].waterbody_parameters
        if waterbody_parameters:
            simulate_waterbodies = waterbody_parameters.break_network_at_waterbodies
            levelpool = waterbody_parameters.level_pool

            if simulate_waterbodies and network_type=='NHDNetwork':
                assert levelpool, 'Waterbody simulation is enabled for NHDNetwork, but levelpool parameters are missing.'
                levelpool_file = levelpool.level_pool_waterbody_parameter_file_path
                assert levelpool_file, 'Waterbody simulation is enabled for NHDNetwork, but no levelpool parameter file is provided.'
                assert levelpool_file.exists(), 'Levelpool file specified, but does not exist in path'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_diffusive_domain(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        if hybrid_parameters:
            run_hybrid = hybrid_parameters.run_hybrid_routing
            if run_hybrid:
                path=hybrid_parameters.diffusive_domain
                assert path, 'Diffusive routing is enabled, but diffusive domain file is missing.'
                assert path.exists(), 'Diffusive domain file specified, but does not exist in path'

        return values

    @root_validator(skip_on_failure=True)
    def check_topobathy_domain(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        if hybrid_parameters:
            use_natl_xsections = hybrid_parameters.use_natl_xsections
            if use_natl_xsections:
                path = hybrid_parameters.topobathy_domain
                assert path, 'Use natural cross-sections is enabled, but topobathy domain file is missing.'
                assert path.exists(), 'Topobathy domain file specified, but does not exist in path'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_refactored(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        if hybrid_parameters:
            run_refactored_network = hybrid_parameters.run_refactored_network
            if run_refactored_network:
                assert hybrid_parameters.refactored_domain, 'Run refactored network is enabled, but refactored domain file is missing.'
                path = hybrid_parameters.refactored_topobathy_domain
                assert path, 'Run refactored network is enabled, but refactored topobathy domain file is missing.'
                assert path.exists(), 'Refactored topobathy domain file specified, but does not exist in path'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_coastal_domain(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        forcing_parameters = values['compute_parameters'].forcing_parameters
        if hybrid_parameters:
            coastal_boundary_input_file = forcing_parameters.coastal_boundary_input_file
            if coastal_boundary_input_file:
                assert coastal_boundary_input_file.exists(), 'Coastal boundary input file specified, but does not exist in path'
                path = hybrid_parameters.coastal_boundary_domain
                assert path, 'Coastal boundary forcing files specified, but coastal boundary domain file is missing.'
                assert path.exists(), 'Coastal boundary domain file specified, but does not exist in path'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_gage_segID_crosswalk_file(cls, values):
        streamflow_DA = values['compute_parameters'].data_assimilation_parameters.streamflow_da
        streamflow_nudging = streamflow_DA.streamflow_nudging
        network_type = values['network_topology_parameters'].supernetwork_parameters.network_type
        if streamflow_nudging and network_type=='NHDNetwork':
            path = streamflow_DA.gage_segID_crosswalk_file
            assert path, 'Streamflow nuding is enabled on NHDNetwork, but gage_segID_crosswalk_file is missing.'
            assert path.exists(), 'gage_segID_crosswalk_file is specified, but does not exist in path'

        return values

    @root_validator(skip_on_failure=True)
    def check_rfc_parameters(cls, values):
        reservoir_da = values['compute_parameters'].data_assimilation_parameters.reservoir_da
        if reservoir_da:
            reservoir_rfc_da = reservoir_da.reservoir_rfc_da
            reservoir_rfc_forecasts = False
            if reservoir_rfc_da:
                reservoir_rfc_forecasts = reservoir_rfc_da.reservoir_rfc_forecasts
                reservoir_rfc_forecasts_time_series_path = reservoir_rfc_da.reservoir_rfc_forecasts_time_series_path
            if reservoir_rfc_forecasts:
                error_message = ''
                network_type = values['network_topology_parameters'].supernetwork_parameters.network_type
                reservoir_parameter_file = reservoir_da.reservoir_parameter_file
                if not reservoir_parameter_file and network_type=='NHDNetwork':
                    error_message += ' Reservoir_parameter_file is missing (and network type is NHDNetwork).'
                if not reservoir_rfc_forecasts_time_series_path:
                    error_message  += ' RFC timeseries path is missing.'
                else:
                    if not os.path.exists(reservoir_rfc_forecasts_time_series_path):
                        error_message += ' reservoir_rfc_forecasts_time_series_path provided does not exist. '                    
                assert not error_message, 'RFC forecast is enabled, but:' + error_message

        return values
    
    @root_validator(skip_on_failure=True)
    def check_usgs_reservoir_da_parameters(cls, values):
        reservoir_da = values['compute_parameters'].data_assimilation_parameters.reservoir_da
        if reservoir_da:
            reservoir_persistence_da = reservoir_da.reservoir_persistence_da
            reservoir_persistence_usgs = False
            if reservoir_persistence_da:
                reservoir_persistence_usgs = reservoir_persistence_da.reservoir_persistence_usgs
                usgs_timeslices_folder = values['compute_parameters'].data_assimilation_parameters.usgs_timeslices_folder
            if reservoir_persistence_usgs:
                error_message = ''
                network_type = values['network_topology_parameters'].supernetwork_parameters.network_type
                reservoir_parameter_file = reservoir_da.reservoir_parameter_file
                if not reservoir_parameter_file and network_type=='NHDNetwork':
                    error_message += ' Reservoir_parameter_file is missing (and network type is NHDNetwork).'
                if not usgs_timeslices_folder:
                    error_message  += ' USGS_timeslices_folder is missing.'
                else:
                    if not os.path.exists(usgs_timeslices_folder):
                        error_message += ' USGS_timeslices_folder path provided does not exist. '
                assert not error_message, 'USGS reservoir DA is enabled, but:' + error_message

        return values
    
    @root_validator(skip_on_failure=True)
    def check_usace_reservoir_da_parameters(cls, values):
        reservoir_da = values['compute_parameters'].data_assimilation_parameters.reservoir_da
        if reservoir_da:
            reservoir_persistence_da = reservoir_da.reservoir_persistence_da
            reservoir_persistence_usace = False
            if reservoir_persistence_da:
                reservoir_persistence_usace = reservoir_persistence_da.reservoir_persistence_usace
                usace_timeslices_folder = values['compute_parameters'].data_assimilation_parameters.usace_timeslices_folder
            if reservoir_persistence_usace:
                error_message = ''
                network_type = values['network_topology_parameters'].supernetwork_parameters.network_type
                reservoir_parameter_file = reservoir_da.reservoir_parameter_file
                if not reservoir_parameter_file and network_type=='NHDNetwork':
                    error_message += ' Reservoir_parameter_file is missing (and network type is NHDNetwork).'
                if not usace_timeslices_folder:
                    error_message  += ' USACE_timeslices_folder is missing.'
                else:
                    if not os.path.exists(usace_timeslices_folder):
                        error_message += ' USACE_timeslices_folder path provided does not exist. '                    
                assert not error_message, 'USACE reservoir DA is enabled, but:' + error_message

        return values
    
    @root_validator(skip_on_failure=True)
    def check_qlat_inputs(cls, values):
        forcing_parameters = values['compute_parameters'].forcing_parameters
        if forcing_parameters:
            qlat_forcing_sets = forcing_parameters.qlat_forcing_sets
            qlat_input_folder = forcing_parameters.qlat_input_folder
            if not qlat_forcing_sets:
                assert qlat_input_folder, 'No qlat_input_folder is specified in the forcing_parameters'
                assert qlat_input_folder.exists(), 'qlat_input_folder is specified, but does not exist in path'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_wrf_hydro_restart_files(cls, values):
        restart_parameters = values['compute_parameters'].restart_parameters
        if restart_parameters:
            wrf_hydro_channel_restart_file = restart_parameters.wrf_hydro_channel_restart_file
            wrf_hydro_waterbody_restart_file = restart_parameters.wrf_hydro_waterbody_restart_file
            if wrf_hydro_channel_restart_file:
                assert wrf_hydro_channel_restart_file.exists(), 'WRF-Hydro channel restart file specified, but does not exist in path'
                path = restart_parameters.wrf_hydro_channel_ID_crosswalk_file
                assert path, 'WRF-Hydro channel restart file provided, but wrf_hydro_channel_ID_crosswalk_file file is missing.'
                assert path.exists(), 'wrf_hydro_channel_ID_crosswalk_file is specified, but does not exist in path'

            if wrf_hydro_waterbody_restart_file:
                error_message = ''
                if not restart_parameters.wrf_hydro_waterbody_ID_crosswalk_file:
                    error_message += ' wrf_hydro_channel_ID_crosswalk_file is missing.'
                if not restart_parameters.wrf_hydro_waterbody_crosswalk_filter_file:
                    error_message += ' wrf_hydro_waterbody_crosswalk_filter_file is missing.'
                assert not error_message, 'WRF-Hydro waterbody_restart file is provided, but:' + error_message

        return values
    
    @root_validator(skip_on_failure=True)
    def check_start_datetime(cls, values):
        restart_parameters = values['compute_parameters'].restart_parameters
        if restart_parameters:
            wrf_hydro_channel_restart_file = restart_parameters.wrf_hydro_channel_restart_file
            lite_channel_restart_file = restart_parameters.lite_channel_restart_file
            if wrf_hydro_channel_restart_file:
                assert wrf_hydro_channel_restart_file.exists(), 'wrf_hydro_channel_restart_file specified, but does not exist in path'
            if lite_channel_restart_file:
                assert lite_channel_restart_file.exists(), 'lite_channel_restart_file specified, but does not exist in path'
            if not (wrf_hydro_channel_restart_file or lite_channel_restart_file):
                assert restart_parameters.start_datetime, 'No start_datetime provided in config file for cold start (no restart files).'

        return values

    @root_validator(skip_on_failure=True)
    def check_flowpath_edge_list(cls, values):
        geo_file_path = values['network_topology_parameters'].supernetwork_parameters.geo_file_path
        flowpath_edge_list = values['network_topology_parameters'].supernetwork_parameters.flowpath_edge_list
        if Path(geo_file_path).suffix=='.json':
            assert flowpath_edge_list, "geo_file_path is json, but no flowpath_edge_list is provided."
            assert Path(flowpath_edge_list).suffix=='.json', "geo_file_path is json, but flowpath_edge_list is a different file type."

        return values
    
    @root_validator(skip_on_failure=True)
    def check_lite_restart_directory(cls, values):
        if values['output_parameters']:
            lite_restart = values['output_parameters'].lite_restart
            if lite_restart is not None:
                lite_restart_directory = lite_restart.lite_restart_output_directory
                assert lite_restart_directory, "lite_restart is present in output parameters, but no lite_restart_output_directory is provided."
                assert lite_restart_directory.exists(), "lite_restart_output_directory is specified, but does not exist in path"

        return values

    @root_validator(skip_on_failure=True)
    def check_nts_dt_stream_output_internal_frequency(cls, values):
        compute_params = values.get('compute_parameters')
        output_params = values.get('output_parameters')
        
        if compute_params and output_params:
            # Directly access ForcingParameters
            forcing_params = compute_params.forcing_parameters
            stream_output = output_params.stream_output

            if forcing_params and stream_output:
                nts = forcing_params.nts
                dt = forcing_params.dt
                internal_freq = stream_output.stream_output_internal_frequency

                # Perform the check
                if nts and dt and internal_freq:
                    result = nts * dt / (internal_freq * 60)
                    if not result.is_integer():
                        raise ValueError("UPDATE nts. Make sure 'nts' times 'dt' divided by ('stream_output_internal_frequency' times 60) is a whole number in your configuration.")

        return values    