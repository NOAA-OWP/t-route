from pydantic import Field, BaseModel, root_validator
from pathlib import Path

from typing import Any, Dict, Optional
from typing_extensions import Self

from .logging_parameters import LoggingParameters
from .network_topology_parameters import NetworkTopologyParameters
from .compute_parameters import ComputeParameters
from .output_parameters import OutputParameters
from .bmi_parameters import BMIParameters
from ._utils import use_strict


class Config(BaseModel):
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
    

    ##########################################################################################
    # Validate that, given certain configuration inputs, other inputs are provided.
    ##########################################################################################

    @root_validator(skip_on_failure=True)
    def check_levelpool_filepath(cls, values):
        network_type = values['network_topology_parameters'].supernetwork_parameters.geo_file_type
        waterbody_parameters = values['network_topology_parameters'].waterbody_parameters
        if waterbody_parameters:
            simulate_waterbodies = waterbody_parameters.break_network_at_waterbodies
            levelpool = waterbody_parameters.level_pool

            if simulate_waterbodies and network_type=='NHDNetwork':
                assert levelpool, 'Waterbody simulation is enabled for NHDNetwork, but levelpool parameters are missing.'
                levelpool_file = levelpool.level_pool_waterbody_parameter_file_path
                assert levelpool_file, 'Waterbody simulation is enabled for NHDNetwork, but no levelpool parameter file is provided.'
            
        return values
    
    @root_validator(skip_on_failure=True)
    def check_rfc_filepath(cls, values):
        network_type = values['network_topology_parameters'].supernetwork_parameters.geo_file_type
        waterbody_parameters = values['network_topology_parameters'].waterbody_parameters
        if waterbody_parameters:
            rfc_parameters = waterbody_parameters.rfc
            if rfc_parameters:
                rfc_forecasts = rfc_parameters.reservoir_rfc_forecasts
                rfc_parameter_file = rfc_parameters.reservoir_parameter_file
                rfc_timeseries_path = rfc_parameters.reservoir_rfc_forecasts_time_series_path

                if rfc_forecasts:
                    assert rfc_timeseries_path, 'RFC forecasts are enabled, but RFC timeseries path is missing.'
                    if network_type=='NHDNetwork':
                        assert rfc_parameter_file, 'RFC forecasts are enabled for NHDNetwork, but no RFC parameter file is provided.'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_diffusive_domain(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        if hybrid_parameters:
            run_hybrid = hybrid_parameters.run_hybrid_routing
            if run_hybrid:
                assert hybrid_parameters.diffusive_domain, 'Diffusive routing is enabled, but diffusive domain file is missing.'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_topobathy_domain(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        if hybrid_parameters:
            use_natl_xsections = hybrid_parameters.use_natl_xsections
            if use_natl_xsections:
                assert hybrid_parameters.topobathy_domain, 'Use natural cross-sections is enabled, but topobathy domain file is missing.'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_refactored(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        if hybrid_parameters:
            run_refactored_network = hybrid_parameters.run_refactored_network
            if run_refactored_network:
                assert hybrid_parameters.refactored_domain, 'Run refactored network is enabled, but refactored domain file is missing.'
                assert hybrid_parameters.refactored_topobathy_domain, 'Run refactored network is enabled, but refactored topobathy domain file is missing.'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_coastal_domain(cls, values):
        hybrid_parameters = values['compute_parameters'].hybrid_parameters
        forcing_parameters = values['compute_parameters'].forcing_parameters
        if hybrid_parameters:
            coastal_boundary_input_file = forcing_parameters.coastal_boundary_input_file
            if coastal_boundary_input_file:
                assert hybrid_parameters.coastal_boundary_domain, 'Coastal boundary forcing files provided, but coastal boundary domain file is missing.'

        return values
    
    @root_validator(skip_on_failure=True)
    def check_gage_segID_crosswalk_file(cls, values):
        streamflow_DA = values['compute_parameters'].data_assimilation_parameters.streamflow_da
        streamflow_nudging = streamflow_DA.streamflow_nudging
        network_type = values['network_topology_parameters'].supernetwork_parameters.geo_file_type
        if streamflow_nudging and network_type=='NHDNetwork':
            assert streamflow_DA.gage_segID_crosswalk_file, 'Streamflow nuding is enabled on NHDNetwork, but gage_segID_crosswalk_file is missing.'

        return values

