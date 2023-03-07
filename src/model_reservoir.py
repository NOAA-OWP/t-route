import numpy as np
import pandas as pd
from pathlib import Path
import yaml

from array import array
#from troute.network.reservoirs.levelpool.levelpool import MC_Levelpool
from troute.network.reservoirs.levelpool.levelpool import MC_Levelpool #, run_lp_c #, update_lp_c
#from troute.network.reservoirs.hybrid.hybrid import MC_Hybrid
#from troute.network.reservoirs.rfc.rfc import MC_RFC


class reservoir_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_levelpool','_inflow','_outflow','_water_elevation','_time', '_time_step',]
        
        if bmi_cfg_file: #TODO: Do we need a config file for this??
            pass
        else:
            pass

        self._time = 0.0
        self._time_step = 300.0
        self._nts = 1
    
    def preprocess_static_vars(self, values: dict):

        lake_number = values['lake_number']
        lake_area = values['lake_area']
        max_depth = values['max_depth']
        orifice_area = values['orifice_area']
        orifice_coefficient = values['orifice_coefficient']
        orifice_elevation = values['orifice_elevation']
        weir_coefficient = values['weir_coefficient']
        weir_elevation = values['weir_elevation']
        weir_length = values['weir_length']
        initial_fractional_depth = values['initial_fractional_depth']
        water_elevation = values['water_elevation']
        
        args = [lake_area, max_depth, orifice_area,
                orifice_coefficient, orifice_elevation,
                weir_coefficient, weir_elevation, weir_length,
                initial_fractional_depth, 0.0, water_elevation]
        
        upstream_ids = array('l', values['upstream_ids'])
        res_type = values['res_type']

        self._levelpool = MC_Levelpool(0, lake_number, upstream_ids, args, res_type)

        # Set inflow
        self._inflow = values['lake_water~incoming__volume_flow_rate']


    def run(self, values: dict,):
        """
        Run this model into the future.
        Run this model into the future, updating the state stored in the provided model dict appropriately.
        Note that the model assumes the current values set for input variables are appropriately for the time
        duration of this update (i.e., ``dt``) and do not need to be interpolated any here.
        Parameters
        ----------
        model: dict
            The model state data structure.
        dt: int
            The number of seconds into the future to advance the model.
        Returns
        -------
        """
        # Run routing
        self._outflow, self._water_elevation = self._levelpool.run(self._inflow, 0.0, self._time_step)
        
        values['lake_water~outgoing__volume_flow_rate'] = self._outflow
        values['lake_surface__elevation'] = self._water_elevation
        
        # update model time
        self._time += self._time_step


# Utility functions -------
"""
def _read_config_file(custom_input_file):
    '''
    Read-in data from user-created configuration file.
    
    Arguments
    ---------
    custom_input_file (str): configuration filepath, .yaml
    
    Returns
    -------
    preprocessing_parameters     (dict): Input parameters re preprocessing
    supernetwork_parameters      (dict): Input parameters re network extent
    waterbody_parameters         (dict): Input parameters re waterbodies
    compute_parameters           (dict): Input parameters re computation settings
    forcing_parameters           (dict): Input parameters re model forcings
    restart_parameters           (dict): Input parameters re model restart
    hybrid_parameters            (dict): Input parameters re diffusive wave model
    output_parameters            (dict): Input parameters re output writing
    parity_parameters            (dict): Input parameters re parity assessment
    data_assimilation_parameters (dict): Input parameters re data assimilation
    '''
    with open(custom_input_file) as custom_file:
        data = yaml.load(custom_file, Loader=yaml.SafeLoader)

    network_topology_parameters = data.get("network_topology_parameters", None)
    supernetwork_parameters = network_topology_parameters.get(
        "supernetwork_parameters", None
    )
    # add attributes when HYfeature network is selected
    if supernetwork_parameters['geo_file_path'][-4:] == "gpkg":
        supernetwork_parameters["title_string"]       = "HY_Features Test"
        supernetwork_parameters["geo_file_path"]      = supernetwork_parameters['geo_file_path']
        supernetwork_parameters["flowpath_edge_list"] = None    
        routelink_attr = {
                        #link????
                        "key": "id",
                        "downstream": "toid",
                        "dx": "length_m",
                        "n": "n",  # TODO: rename to `manningn`
                        "ncc": "nCC",  # TODO: rename to `mannningncc`
                        "s0": "So",
                        "bw": "BtmWdth",  # TODO: rename to `bottomwidth`
                        #waterbody: "NHDWaterbodyComID",
                        "tw": "TopWdth",  # TODO: rename to `topwidth`
                        "twcc": "TopWdthCC",  # TODO: rename to `topwidthcc`
                        "alt": "alt",
                        "musk": "MusK",
                        "musx": "MusX",
                        "cs": "ChSlp"  # TODO: rename to `sideslope`
                        }
        supernetwork_parameters["columns"]             = routelink_attr 
        supernetwork_parameters["waterbody_null_code"] = -9999
        supernetwork_parameters["terminal_code"]       =  0
        supernetwork_parameters["driver_string"]       = "NetCDF"
        supernetwork_parameters["layer_string"]        = 0
        
    preprocessing_parameters = network_topology_parameters.get(
        "preprocessing_parameters", {}
    )        
    #waterbody_parameters = network_topology_parameters.get(
    #    "waterbody_parameters", None
    #)
    waterbody_parameters = network_topology_parameters.get(
        "waterbody_parameters", {}
    )
    compute_parameters = data.get("compute_parameters", {})
    forcing_parameters = compute_parameters.get("forcing_parameters", {})
    restart_parameters = compute_parameters.get("restart_parameters", {})
    hybrid_parameters = compute_parameters.get("hybrid_parameters", {})
    data_assimilation_parameters = compute_parameters.get(
        "data_assimilation_parameters", {}
    )
    output_parameters = data.get("output_parameters", {})
    parity_parameters = output_parameters.get("wrf_hydro_parity_check", {})

    return (
        preprocessing_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        hybrid_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    )
"""
