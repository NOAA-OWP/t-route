
import pandas as pd
import numpy as np
from datetime import datetime

from array import array
from troute.network.reservoirs.levelpool.levelpool import MC_Levelpool
from troute.routing.fast_reach.reservoir_hybrid_da import reservoir_hybrid_da

#TODO: Once PR#605 is merged, uncomment and use the updated DA objects
#from troute.DataAssimilation import PersistenceDA, RFCDA

class reservoir_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_levelpool','_inflow','_outflow','_water_elevation', '_res_type',
                     '_update_time', '_prev_persisted_flow', '_persistence_update_time',
                     '_persistence_index', '_time', '_time_step',]
        
        if bmi_cfg_file: #TODO: Do we need a config file for this??
            pass
        else:
            pass

        self._time = 0.0
        self._time_step = 300.0
        self._t0 = datetime(2018, 6, 1, 0, 0, 0) #placeholder...read from config file?
    
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
        self._res_type = values['res_type']

        self._levelpool = MC_Levelpool(0, lake_number, upstream_ids, args, self._res_type)

        # Set inflow
        self._inflow = values['lake_water~incoming__volume_flow_rate']

        # Set data assimilation parameters
        if self._res_type==2 or self._res_type==3:
            self._update_time = 0
            self._prev_persisted_outflow = np.nan
            self._persistence_index = 0
            self._persistence_update_time = 0


    def run(self, values: dict,):
        """
        Run this model into the future, updating the state stored in the provided model dict appropriately.
    
        Parameters
        ----------
        values: dict
            The model state data structure.
            
        Returns
        -------
        """
        # Get water elevation before levelpool calculation
        initial_water_elevation = self._levelpool.water_elevation

        # Run routing
        self._outflow, water_elevation = self._levelpool.run(self._inflow, 0.0, self._time_step)
        
        # Data Assimilation
        if self._res_type==2 or self._res_type==3:
            #TODO: figure out what format 'gage_time' will be when passed from model engine...
            #gage_time = np.array((values['gage_time']-self._t0).total_seconds())

            (
                new_outflow,
                new_persisted_outflow,
                new_water_elevation, 
                new_update_time, 
                new_persistence_index, 
                new_persistence_update_time
            ) = reservoir_hybrid_da(
                values['lake_number'],        # lake identification number
                values['gage_observations'],  # gage observation values (cms)
                values['gage_time'], #gage_time,                    # gage observation times (sec)
                self._time,                   # model time (sec)
                self._prev_persisted_outflow, # previously persisted outflow (cms)
                self._persistence_update_time,
                self._persistence_index,      # number of sequentially persisted update cycles
                self._outflow,                # levelpool simulated outflow (cms)
                self._inflow,                 # waterbody inflow (cms)
                self._time_step,              # model timestep (sec)
                values['lake_area'],          # waterbody surface area (km2)
                values['max_depth'],          # max waterbody depth (m)
                values['orifice_elevation'],  # orifice elevation (m)
                initial_water_elevation,      # water surface el., previous timestep (m)
                48.0,                         # gage lookback hours (hrs)
                self._update_time,            # waterbody update time (sec)
            )
            
            # update levelpool water elevation state
            water_elevation = self._levelpool.assimilate_elevation(new_water_elevation)
            
            # change reservoir_outflow
            self._outflow = new_outflow
            
            # update USGS DA reservoir state arrays
            self._update_time = new_update_time
            self._prev_persisted_outflow = new_persisted_outflow
            self._persistence_index = new_persistence_index
            self._persistence_update_time = new_persistence_update_time

        # Set output variables
        values['lake_water~outgoing__volume_flow_rate'] = self._outflow
        values['lake_surface__elevation'] = water_elevation

        # update model time
        self._time += self._time_step

