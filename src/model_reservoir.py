
import pandas as pd
import numpy as np
from datetime import datetime

from array import array
from troute.network.reservoirs.levelpool.levelpool import MC_Levelpool
from troute.routing.fast_reach.reservoir_hybrid_da import reservoir_hybrid_da
from troute.routing.fast_reach.reservoir_RFC_da import reservoir_RFC_da

class reservoir_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_levelpool','_inflow','_outflow','_water_elevation', '_res_type',
                     '_update_time', '_prev_persisted_flow', '_persistence_update_time',
                     '_persistence_index', '_time', '_time_step','_dynamic_res_type',
                     '_assimilated_value', '_assimilated_source_file']
        
        if bmi_cfg_file: #TODO: Do we need a config file for this??
            pass
        else:
            pass

        self._time = 0.0
        self._time_step = 300.0
        self._t0 = datetime(2018, 6, 1, 0, 0, 0) #placeholder...read from config file?
    
    def preprocess_static_vars(self, values: dict):

        lake_number = values['waterbody_id']
        lake_area = values['LkArea']# * 1e6
        max_depth = values['LkMxE']
        orifice_area = values['OrificeA']
        orifice_coefficient = values['OrificeC']
        orifice_elevation = values['OrificeE']
        weir_coefficient = values['WeirC']
        weir_elevation = values['WeirE']
        weir_length = values['WeirL']
        initial_fractional_depth = values['ifd']
        water_elevation = values['lake_surface__elevation']
        
        args = [lake_area, max_depth, orifice_area,
                orifice_coefficient, orifice_elevation,
                weir_coefficient, weir_elevation, weir_length,
                initial_fractional_depth, 0.0, water_elevation]
        
        upstream_ids = array('l', values['upstream_ids'])
        self._res_type = values['reservoir_type']

        self._levelpool = MC_Levelpool(0, lake_number, upstream_ids, args, 1)
        
        # Set inflow
        self._inflow = values['lake_water~incoming__volume_flow_rate']
        self._inflow_list = self._inflow.tolist()
        self._outflow_list = []
        self._water_elevation_list = []

        # Set data assimilation parameters
        if self._res_type==2 or self._res_type==3:
            self._update_time = 0
            self._prev_persisted_outflow = np.nan
            self._persistence_index = 0
            self._persistence_update_time = 0
        
        if self._res_type==4 or self._res_type==5:
            self._update_time = 0
            self._timeseries_idx = values['da_idx'] - 1
            self._da_time_step = values['time_step']
            self._persist_seconds = values['rfc_forecast_persist_seconds']


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
        # Check if new inflow values have been provided
        #TODO Come up with more clever way to do this, this seems breakable...
        if len(self._inflow_list)==0:
            self._inflow_list = self._inflow.tolist()

        # Get water elevation before levelpool calculation
        initial_water_elevation = self._levelpool.water_elevation

        # Run routing
        inflow = self._inflow_list.pop(0)
        self._outflow, water_elevation = self._levelpool.run(inflow, 0.0, self._time_step)

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
                self._levelpool.lake_number,        # lake identification number
                values['gage_observations'],        # gage observation values (cms)
                values['gage_time'], #gage_time,    # gage observation times (sec)
                self._time,                         # model time (sec)
                self._prev_persisted_outflow,       # previously persisted outflow (cms)
                self._persistence_update_time,
                self._persistence_index,            # number of sequentially persisted update cycles
                self._outflow,                      # levelpool simulated outflow (cms)
                inflow,                             # waterbody inflow (cms)
                self._time_step,                    # model timestep (sec)
                self._levelpool.lake_area,          # waterbody surface area (km2)
                self._levelpool.max_depth,          # max waterbody depth (m)
                self._levelpool.orifice_elevation,  # orifice elevation (m)
                initial_water_elevation,            # water surface el., previous timestep (m)
                24.0,                               # gage lookback hours (hrs)
                self._update_time,                  # waterbody update time (sec)
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
        
        elif self._res_type==4 or self._res_type==5:
            (
                new_outflow, 
                new_water_elevation, 
                new_update_time,
                new_timeseries_idx,
                dynamic_reservoir_type, 
                assimilated_value, 
                #assimilated_source_file,
            ) = reservoir_RFC_da(
                self._levelpool.lake_number,        # lake identification number
                values['gage_observations'],  # gage observation values (cms)
                self._timeseries_idx,         # index of for current time series observation
                values['synthetic_flag'],     # boolean flags indicating synthetic values
                self._time_step,              # routing period (sec)
                self._time,                   # model time (sec)
                self._update_time,                  # time to advance to next time series index
                self._da_time_step,           # frequency of DA observations (sec)
                self._persist_seconds,        # max seconds RFC forecasts will be used/persisted
                self._res_type,               # reservoir type
                self._inflow,                 # waterbody inflow (cms)
                initial_water_elevation,      # water surface el., previous timestep (m)
                self._outflow,                # levelpool simulated outflow (cms)
                water_elevation,              # levelpool simulated water elevation (m)
                self._levelpool.lake_area,          # waterbody surface area (km2)
                self._levelpool.max_depth,          # max waterbody depth (m)
            )
            
            # update levelpool water elevation state
            water_elevation = self._levelpool.assimilate_elevation(new_water_elevation)
            
            # change reservoir_outflow
            self._outflow = new_outflow

            # update DA reservoir state parameters
            self._update_time = new_update_time
            self._timeseries_idx = new_timeseries_idx
            self._dynamic_res_type = dynamic_reservoir_type
            self._assimilated_value = assimilated_value
            #self._assimilated_source_file = assimilated_source_file

        # Set output variables
        values['lake_water~outgoing__volume_flow_rate'] = self._outflow
        values['lake_surface__elevation'] = water_elevation

        # Append outflow list
        self._outflow_list.append(self._outflow)
        self._water_elevation_list.append(water_elevation)

        # update model time
        self._time += self._time_step

