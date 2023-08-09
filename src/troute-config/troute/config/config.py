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

    @root_validator
    def check_levelpool_filepath(cls, values):
        network_type = values['network_topology_parameters'].supernetwork_parameters.geo_file_type
        waterbody_parameters = values['network_topology_parameters'].waterbody_parameters
        simulate_waterbodies = waterbody_parameters.break_network_at_waterbodies
        levelpool = waterbody_parameters.level_pool

        if simulate_waterbodies and network_type=='NHDNetwork':
            assert levelpool, 'Waterbody simulation is enabled for NHDNetwork, but levelpool parameters are missing'
            levelpool_file = levelpool.level_pool_waterbody_parameter_file_path
            assert levelpool_file, 'Waterbody simulation is enabled for NHDNetwork, but no levelpool file is provided'
            
        return values

