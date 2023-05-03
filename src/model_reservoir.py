
import pandas as pd
import numpy as np
from datetime import datetime

from array import array
from troute.network.reservoirs.levelpool.levelpool import MC_Levelpool
from troute.routing.fast_reach.reservoir_hybrid_da import reservoir_hybrid_da
from troute.routing.fast_reach.reservoir_RFC_da import reservoir_RFC_da, reservoir_RFC_da_v2

class reservoir_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_levelpool','_inflow','_outflow','_water_elevation', '_res_type',
                     '_update_time', '_prev_persisted_flow', '_persistence_update_time',
                     '_persistence_index', '_time', '_time_step','_dynamic_res_type',
                     '_assimilated_value', '_assimilated_source_file',
                     '_levelpool_outflow,', 
                     '_current_time', '_timeseries_update_time',                     # RFC variables
                     '_timeseries_idx', '_rfc_gage_id', '_rfc_timeseries_folder',    # RFC variables
                     '_time_step_seconds', '_total_counts','_timeseries_discharges', '_use_RFC'] # RFC variables
        
        if bmi_cfg_file: #TODO: Do we need a config file for this??
            pass
        else:
            pass

        self._time = 0.0
        self._time_step = 300.0
        self._t0 = datetime(2021, 10, 20, 13, 0, 0) #placeholder...read from config file?
    
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
            #self._update_time = 0
            #self._timeseries_idx = values['da_idx'] - 1
            #self._da_time_step = values['time_step']
            #self._persist_seconds = values['rfc_forecast_persist_seconds']
            self._current_time=0
            self._timeseries_update_time=0 
            self._timeseries_idx = 0
            self._rfc_gage_id = "KNFC1"
            self._rfc_timeseries_folder = "./../test/BMI/rfc_timeseries/"
            self._time_step_seconds = 0
            self._total_counts = 0
            self._timeseries_discharges =np.zeros(432) # maximum numbers of hourly data (3d obs + 15d forecasts)
            self._use_RFC = True
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
        # ISSUE: isn't it used only for water elevation at model start time at t0?
        initial_water_elevation = self._levelpool.water_elevation

        # TODO: For RFC DA, shouldn't water elevation at t0 be computed value from the previous time step? 
        if self._time == 0:
            self._water_elevation = initial_water_elevation

        # Run routing
        inflow = self._inflow_list.pop(0)
        #TOCONSIDER: not to cofused between levelpool and res DA outputs
        self._levelpool_outflow, levelpool_water_elevation = self._levelpool.run(inflow, 0.0, self._time_step)

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
                self._levelpool_outflow, #self._outflow   # levelpool simulated outflow (cms)
                inflow,                             # waterbody inflow (cms)
                self._time_step,                    # model timestep = time step of inflow to reservoir (sec)
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
            '''
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
            '''
            model_start_date = self._t0.strftime("%Y-%m-%d_%H:%M:%S")
            
            (
                new_outflow, 
                new_water_elevation,
                new_current_time, 
                new_timeseries_update_time,
                new_timeseries_idx,
                dynamic_reservoir_type, 
                assimilated_value, 
                time_step_seconds,
                total_counts,
                timeseries_discharges,
                use_RFC,
            ) = reservoir_RFC_da_v2(
                self._levelpool.lake_number,                     # lake identification number
                self._rfc_gage_id,                               # RFC gage ID
                model_start_date,                                # simulation start date
                self._time_step,                                 # model time step = time step of inflow to reservoir = routing period (sec)
                self._current_time,                              # rfc DA simulation time starting from 0 with increments by ch.routing time step
                self._timeseries_update_time,                    # defines location of time interval of timeseries corresponding to _current_time
                self._timeseries_idx,                            # locates discharge in timeseries corresponding to _current_time
                values['rfc_timeseries_offset_hours'],           # offset hours ahead of model start date for searching RFC files backward in time
                values['rfc_forecast_persist_days'],             # max number of days that RFC-supplied forecast will be used/persisted in simulation
                self._rfc_timeseries_folder,                     # path to folder containing RFC timeseires files
                self._res_type,                                  # reservoir type
                inflow,                                          # waterbody inflow (a single value) (cms)
                self._water_elevation, # initial_water_elevation # water surface el., previous timestep (m)
                self._levelpool_outflow,                         # levelpool simulated outflow (cms)
                levelpool_water_elevation,                       # levelpool simulated water elevation (m)
                self._levelpool.lake_area,                       # waterbody surface area (km2)
                self._levelpool.max_depth,                       # ** actually max elevation, not depth, of waterbody (m)
                self._time_step_seconds,                          # seconds of 'timeSteps' of selected RFCTimeSeries file
                self._total_counts,                               # total number of discharge data of selected RFCTimeSeries file
                self._timeseries_discharges,                      # discharge timeseries from a select RFCTimeSeries file
                self._use_RFC,                                    # True: use RFC DA
            )    

            # update levelpool water elevation state
            #TODO: Isn't using the second line below more explicit? 
            #water_elevation = self._levelpool.assimilate_elevation(new_water_elevation)
            water_elevation = new_water_elevation
            
            # change reservoir_outflow
            self._outflow = new_outflow

            # update DA reservoir state parameters
            #self._update_time = new_update_time
            self._current_time           = new_current_time
            self._timeseries_update_time = new_timeseries_update_time
            self._timeseries_idx         = new_timeseries_idx
            self._water_elevation        = new_water_elevation 
            self._dynamic_res_type       = dynamic_reservoir_type
            self._assimilated_value      = assimilated_value
            #self._assimilated_source_file = assimilated_source_file
            self._time_step_seconds        = time_step_seconds
            self._total_counts             = total_counts
            self._timeseries_discharges    = timeseries_discharges
            self._use_RFC                  = use_RFC

        # Set output variables
        values['lake_water~outgoing__volume_flow_rate'] = self._outflow
        values['lake_surface__elevation'] = water_elevation

        # Append outflow list
        self._outflow_list.append(self._outflow)
        self._water_elevation_list.append(water_elevation)

        # update model time
        self._time += self._time_step

