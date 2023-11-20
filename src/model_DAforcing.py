import numpy as np
import pandas as pd
import yaml
from datetime import datetime, timedelta
from joblib import delayed, Parallel
import xarray as xr
import glob
import pathlib

from bmi_array2df import *
import bmi_array2df as a2df

from bmi_df2array import *
import bmi_df2array as df2a

from troute.routing.fast_reach.reservoir_RFC_da import _validate_RFC_data

from nwm_routing.log_level_set import log_level_set
from troute.config import Config

class DAforcing_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_data_assimilation_parameters', '_forcing_parameters', '_compute_parameters',
                     '_output_parameters', '_usgs_df', 'reservoir_usgs_df', 'reservoir_usace_df', 
                     '_rfc_timeseries_df', '_lastobs_df', '_t0', '_q0', '_waterbody_df',
                     '_dateNull', 
                     '_datesSecondsArray_usgs', '_nDates_usgs', '_stationArray_usgs', 
                     '_stationStringLengthArray_usgs', '_nStations_usgs',
                     '_usgs_Array',
                     '_datesSecondsArray_reservoir_usgs', '_nDates_reservoir_usgs',
                     '_stationArray_reservoir_usgs', '_stationStringLengthArray_reservoir_usgs',
                     '_nStations_reservoir_usgs',
                     '_usgs_reservoir_Array',
                     '_datesSecondsArray_reservoir_usace', '_nDates_reservoir_usace',
                     '_stationArray_reservoir_usace', '_stationStringLengthArray_reservoir_usace',
                     '_nStations_reservoir_usace', 
                     '_usace_reservoir_Array',
                     '_rfc_da_timestep', '_rfc_totalCounts', '_rfc_synthetic_values',
                     '_rfc_discharges', '_rfc_timeseries_idx', '_rfc_use_rfc',
                     '_rfc_Datetime', '_rfc_timeSteps', '_rfc_StationId_array',
                     '_rfc_StationId_stringLengths', '_rfc_List_array',
                     '_rfc_List_stringLengths',
                     '_lastObs_gageArray', '_lastObs_gageStringLengths', '_lastObs_timeSince',
                     '_lastObs_discharge',
                     '_q0_columnArray', '_q0_columnLengthArray', '_q0_nCol', '_q0_indexArray',
                     '_q0_nIndex', '_q0_Array',
                     '_waterbodyLR_columnArray', '_waterbodyLR_columnLengthArray', 
                     '_waterbodyLR_nCol', '_waterbodyLR_indexArray', '_waterbodyLR_nIndex',
                     '_waterbodyLR_Array'

                     
                     ]

        if bmi_cfg_file:
            (compute_parameters,
             forcing_parameters,
             data_assimilation_parameters,
             output_parameters,
            ) = _read_config_file(bmi_cfg_file)
            
            self._compute_parameters = compute_parameters
            self._forcing_parameters = forcing_parameters
            self._data_assimilation_parameters = data_assimilation_parameters
            self._output_parameters = output_parameters

            self._t0 = self._compute_parameters['restart_parameters']['start_datetime']

            #############################
            # Read DA files:
            #############################
            nudging = data_assimilation_parameters.get('streamflow_da', {}).get('streamflow_nudging', False)
            
            usgs_persistence = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_persistence_da', {}).get('reservoir_persistence_usgs', False)
            usace_persistence = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_persistence_da', {}).get('reservoir_persistence_usace', False)
            rfc = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_rfc_da', {}).get('reservoir_rfc_forecasts', False)

            qc_threshold = data_assimilation_parameters.get('qc_threshold')
            cpu_pool = compute_parameters.get('cpu_pool')

            # Produce list of datetimes to search for timeslice files
            lookback_hrs = data_assimilation_parameters.get('timeslice_lookback_hours')
            start_datetime = compute_parameters.get('restart_parameters').get('start_datetime')
            dt = compute_parameters.get('forcing_parameters').get('dt')
            nts = compute_parameters.get('forcing_parameters').get('nts')
            timeslice_start = start_datetime - timedelta(hours=lookback_hrs)
            timeslice_end = start_datetime + timedelta(seconds=dt*nts)
            delta = timedelta(minutes=15)
            timeslice_dates = []
            while timeslice_start <= timeslice_end:
                timeslice_dates.append(timeslice_start.strftime('%Y-%m-%d_%H:%M:%S'))
                timeslice_start += delta
            
            # Create empty default dataframes:
            self._usgs_df = pd.DataFrame()
            self._reservoir_usgs_df = pd.DataFrame()
            self._reservoir_usace_df = pd.DataFrame()
            self._rfc_timeseries_df = pd.DataFrame()
            self._lastobs_df = pd.DataFrame()

            # USGS Observations
            if nudging or usgs_persistence:
                usgs_timeslice_path = str(data_assimilation_parameters.get('usgs_timeslices_folder'))
                if nudging:
                    self._usgs_df = _read_timeslice_files(usgs_timeslice_path,
                                                          timeslice_dates,
                                                          qc_threshold,
                                                          dt,
                                                          cpu_pool,)
                    self._reservoir_usgs_df = (
                        self._usgs_df.
                        transpose().
                        resample('15min').asfreq().
                        transpose()
                        )
                else:
                    self._usgs_df = pd.DataFrame()
                    self._reservoir_usgs_df = _read_timeslice_files(usgs_timeslice_path,
                                                                    timeslice_dates,
                                                                    qc_threshold,
                                                                    900, #15 minutes
                                                                    cpu_pool,)

            # USACE Observations        
            if usace_persistence:
                usace_timeslice_path = str(data_assimilation_parameters.get('usace_timeslices_folder'))
                self._reservoir_usace_df = _read_timeslice_files(usace_timeslice_path, 
                                                                 timeslice_dates,
                                                                 qc_threshold,
                                                                 900, #15 minutes
                                                                 cpu_pool,)

            # Produce list of datetimes to search for timeseries files
            rfc_parameters = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_rfc_da', {})
            lookback_hrs = rfc_parameters.get('reservoir_rfc_forecasts_lookback_hours')
            offset_hrs = rfc_parameters.get('reservoir_rfc_forecasts_offset_hours')
            timeseries_end = start_datetime + timedelta(hours=offset_hrs)
            timeseries_start = timeseries_end - timedelta(hours=lookback_hrs)
            delta = timedelta(hours=1)
            timeseries_dates = []
            while timeseries_start <= timeseries_end:
                timeseries_dates.append(timeseries_start.strftime('%Y-%m-%d_%H'))
                timeseries_start += delta
            rfc_forecast_persist_days = rfc_parameters.get('reservoir_rfc_forecast_persist_days')
            final_persist_datetime = start_datetime + timedelta(days=rfc_forecast_persist_days)

            # RFC Observations
            if rfc:
                rfc_timeseries_path = str(rfc_parameters.get('reservoir_rfc_forecasts_time_series_path'))
                self._rfc_timeseries_df = _read_timeseries_files(rfc_timeseries_path, timeseries_dates, start_datetime, final_persist_datetime)

            # Lastobs
            lastobs_file = data_assimilation_parameters.get('streamflow_da', {}).get('lastobs_file', False)

            if lastobs_file:
                self._lastobs_df = _read_lastobs_file(lastobs_file)

            # save time reference: start_datetime
            self._dateNull = start_datetime

            # read in metadata for BMI compliant arrays:

            self._values = {} 

            self._lastObs_gageArray = np.zeros(0)
            self._lastObs_gageStringLengths = np.zeros(0)
            self._lastObs_timeSince = np.zeros(0)
            self._lastObs_discharge = np.zeros(0)

            if not self._lastobs_df.empty:

                if lastobs_file:

                    (_lastObs_gageArray, _lastObs_gageStringLengths, \
                    _lastObs_timeSince, _lastObs_discharge) = \
                    df2a._bmi_disassemble_lastObs (self._lastobs_df) 

                    self._lastObs_gageArray = _lastObs_gageArray
                    self._lastObs_gageStringLengths = _lastObs_gageStringLengths
                    self._lastObs_timeSince = _lastObs_timeSince
                    self._lastObs_discharge = _lastObs_discharge

            # USGS Observations
            self._datesSecondsArray_usgs = np.zeros(0)
            self._nDates_usgs = np.zeros(0)
            self._stationArray_usgs = np.zeros(0)
            self._stationStringLengthArray_usgs = np.zeros(0)
            self._nStations_usgs = np.zeros(0)
            self._usgsArray = np.zeros(0)

            if not self._usgs_df.empty:
                # 
                # USGS dataframe: 
                # 
                # Dates related:
                # - datesSecondsArray_usgs: "column names" of usgs dataframe: dates
                #        in seconds relative to start_datetime
                # - nDates_usgs: number of "dates" columns
                # 
                # Station ID related:
                # - stationArray_usgs: "row names" of usgs dataframe: all station
                #        names decoded into arrays of one ASCII code per string,
                #        attached together
                # - stationStringLengthArray_usgs: array of lengths of station IDs,
                #        in case of mixed arrays (not all station IDs have equal lengths)
                # - nStations_usgs: number of stations
                #
                ( _datesSecondsArray_usgs, _nDates_usgs, _stationArray_usgs, \
                    _stationStringLengthArray_usgs, _nStations_usgs) \
                    = df2a._time_stations_from_df(self._usgs_df,start_datetime)
                # save metadata in class instances
                self._datesSecondsArray_usgs = _datesSecondsArray_usgs
                self._nDates_usgs = _nDates_usgs
                self._stationArray_usgs = _stationArray_usgs
                self._stationStringLengthArray_usgs = _stationStringLengthArray_usgs
                self._nStations_usgs = _nStations_usgs
                # flatten the actual USGS dataframe into a numpy ndarray
                _usgsArray = df2a._flatten_array(self._usgs_df, np.float32)
                # ... and save it with the class instance
                self._usgsArray = _usgsArray

            # USGS Reservoir Observations
            self._datesSecondsArray_reservoir_usgs = np.zeros(0)
            self._nDates_reservoir_usgs = np.zeros(0)
            self._stationArray_reservoir_usgs = np.zeros(0)
            self._stationStringLengthArray_reservoir_usgs = np.zeros(0)
            self._nStations_reservoir_usgs = np.zeros(0)
            self._reservoirUsgsArray = np.zeros(0)

            if not self._reservoir_usgs_df.empty:
                # 
                # Reservoir USGS and USACE dataframes: Array structures of 
                # _reservoir_usgs_df and _reservoir_usace_df are
                # fundamentally the same as usgs_df: input/output analogous.
                # 
                # In the following comments, usgs shall also represent usace
                #
                # Dates related:
                # - datesSecondsArray_reservoir_usgs
                # - nDates_usgs
                # 
                # Station ID related:
                # - stationArray_usgs 
                # - stationStringLengthArray_usgs 
                # - nStations_usgs 
                #
                ( _datesSecondsArray_reservoir_usgs, _nDates_reservoir_usgs, \
                     _stationArray_reservoir_usgs, _stationStringLengthArray_reservoir_usgs, \
                    _nStations_reservoir_usgs) \
                    = df2a._time_stations_from_df(self._reservoir_usgs_df,start_datetime)
                # save metadata in class instances
                self._datesSecondsArray_reservoir_usgs = _datesSecondsArray_reservoir_usgs
                self._nDates_reservoir_usgs = _nDates_reservoir_usgs
                self._stationArray_reservoir_usgs = _stationArray_reservoir_usgs
                self._stationStringLengthArray_reservoir_usgs = _stationStringLengthArray_reservoir_usgs
                self._nStations_reservoir_usgs = _nStations_reservoir_usgs
                # flatten the actual USGS datafrane into a numpy ndarray
                _reservoirUsgsArray = df2a._flatten_array(self._reservoir_usgs_df, np.float32)
                # ... and save it with the class instance
                self._reservoirUsgsArray = _reservoirUsgsArray            

            # USACE Reservoir Observations    
            self._datesSecondsArray_reservoir_usace = np.zeros(0)
            self._nDates_reservoir_usace = np.zeros(0)
            self._stationArray_reservoir_usace = np.zeros(0)
            self._stationStringLengthArray_reservoir_usace = np.zeros(0)
            self._nStations_reservoir_usace = np.zeros(0)
            self._reservoirUsaceArray = np.zeros(0)
    
            if not self._reservoir_usace_df.empty:

                # see detailed comments in USGS branch
                ( _datesSecondsArray_reservoir_usace, _nDates_reservoir_usace, \
                    _stationArray_reservoir_usace, _stationStringLengthArray_reservoir_usace, \
                    _nStations_reservoir_usace) \
                    = df2a._time_stations_from_df(self._reservoir_usace_df,start_datetime)
                # save metadata in class instances
                self._datesSecondsArray_reservoir_usace = _datesSecondsArray_reservoir_usace
                self._nDates_reservoir_usace = _nDates_reservoir_usace
                self._stationArray_reservoir_usace = _stationArray_reservoir_usace
                self._stationStringLengthArray_reservoir_usace = _stationStringLengthArray_reservoir_usace
                self._nStations_reservoir_usace = _nStations_reservoir_usace
                # flatten the actual USACE datafrane into a numpy ndarray
                _reservoirUsaceArray = df2a._flatten_array(self._reservoir_usace_df, np.float32)
                # ... and save it with the class instance
                self._reservoirUsaceArray = _reservoirUsaceArray  


            # RFC Timeseries    
            self._rfc_da_timestep = np.zeros(0)    
            self._rfc_totalCounts = np.zeros(0)   
            self._rfc_synthetic_values = np.zeros(0)   
            self._rfc_discharges = np.zeros(0)   
            self._rfc_timeseries_idx = np.zeros(0)   
            self._rfc_use_rfc = np.zeros(0)   
            self._rfc_Datetime = np.zeros(0)   
            self._rfc_timeSteps = np.zeros(0)   
            self._rfc_StationId_array = np.zeros(0) 
            self._rfc_StationId_stringLengths = np.zeros(0) 
            self._rfc_List_array = np.zeros(0) 
            self._rfc_List_stringLengths = np.zeros(0) 

            if not self._rfc_timeseries_df.empty:

                (_rfc_da_timestep, _rfc_totalCounts, _rfc_synthetic_values, _rfc_discharges, \
                    _rfc_timeseries_idx, _rfc_use_rfc, _rfc_Datetime, _rfc_timeSteps, \
                    _rfc_StationId_array, _rfc_StationId_stringLengths, _rfc_List_array, \
                    _rfc_List_stringLengths) = \
                    df2a._bmi_disassemble_rfc_timeseries (self._rfc_timeseries_df, start_datetime)
                # save all data in class instance
                self._rfc_da_timestep = _rfc_da_timestep
                self._rfc_totalCounts = _rfc_totalCounts
                self._rfc_synthetic_values= _rfc_synthetic_values
                self._rfc_discharges = _rfc_discharges
                self._rfc_timeseries_idx = _rfc_timeseries_idx
                self._rfc_use_rfc = _rfc_use_rfc
                self._rfc_Datetime = _rfc_Datetime
                self._rfc_timeSteps = _rfc_timeSteps
                self._rfc_StationId_array = _rfc_StationId_array
                self._rfc_StationId_stringLengths = _rfc_StationId_stringLengths
                self._rfc_List_array = _rfc_List_array
                self._rfc_List_stringLengths = _rfc_List_stringLengths

            #############################
            # Read Restart files:
            #############################

            # Create empty default dataframes:
            self._q0 = pd.DataFrame()
            self._waterbody_df = pd.DataFrame()

            # q0    
            self._q0_columnArray = np.zeros(0)    
            self._q0_columnLengthArray = np.zeros(0)   
            self._q0_nCol = np.zeros(0)   
            self._q0_indexArray = np.zeros(0)   
            self._q0_nIndex = np.zeros(0)   
            self._q0_Array = np.zeros(0) 

            # waterbody_df
            self._waterbodyLR_columnArray = np.zeros(0)    
            self._waterbodyLR_columnLengthArray = np.zeros(0)    
            self._waterbodyLR_nCol = np.zeros(0)    
            self._waterbodyLR_indexArray = np.zeros(0)    
            self._waterbodyLR_nIndex = np.zeros(0)    
            self._waterbodyLR_Array = np.zeros(0)    

            lite_restart_file = self._compute_parameters['restart_parameters']['lite_channel_restart_file']
            if lite_restart_file:

                self._q0, self._t0 = _read_lite_restart(lite_restart_file)

                if not self._q0.empty:

                    (_q0_columnArray, _q0_columnLengthArray, _q0_nCol, \
                        _q0_indexArray, _q0_nIndex, _q0_Array) = \
                        df2a._bmi_disassemble_lite_restart (self._q0,np.float32)

                    self._q0_columnArray = _q0_columnArray
                    self._q0_columnLengthArray = _q0_columnLengthArray
                    self._q0_nCol = _q0_nCol
                    self._q0_indexArray = _q0_indexArray
                    self._q0_nIndex = _q0_nIndex
                    self._q0_Array = _q0_Array

            lite_restart_file = self._compute_parameters['restart_parameters']['lite_waterbody_restart_file']
            
            if lite_restart_file:

                self._waterbody_df, _ = _read_lite_restart(lite_restart_file)

                if not self._waterbody_df.empty:

                    (_waterbodyLR_columnArray, _waterbodyLR_columnLengthArray, _waterbodyLR_nCol, \
                        _waterbodyLR_indexArray, _waterbodyLR_nIndex, _waterbodyLR_Array) = \
                        df2a._bmi_disassemble_lite_restart (self._waterbody_df,np.float64)

                    self._waterbodyLR_columnArray = _waterbodyLR_columnArray
                    self._waterbodyLR_columnLengthArray = _waterbodyLR_columnLengthArray
                    self._waterbodyLR_nCol = _waterbodyLR_nCol
                    self._waterbodyLR_indexArray = _waterbodyLR_indexArray
                    self._waterbodyLR_nIndex = _waterbodyLR_nIndex
                    self._waterbodyLR_Array = _waterbodyLR_Array

        else:

            raise(RuntimeError("No config file provided."))


    def run(self, values: dict):
        """
        Write t-route output values to files. Namely, create restart files (channel and waterbody), 
        lastobs files, and flowveldepth files to pass to coastal hydraulics models.
        Parameters
        ----------
        values: dict
            The static and dynamic values for the model.
        Returns
        -------
        """
        output_parameters = self._output_parameters
        da_parameters = self._data_assimilation_parameters

        if output_parameters['lite_restart'] is not None:
            _write_lite_restart(
                values,
                self._t0,
                output_parameters['lite_restart']
            )   

        lastobs_output_folder = da_parameters.get('streamflow_da',{}).get('lastobs_output_folder', None)
        if lastobs_output_folder:
            _write_lastobs(
                values,
                self._t0,
                lastobs_output_folder
            )
            
        if output_parameters.get('stream_output',{}).get('qvd_nudge'):
            #write flowveldepth file
            pass


# Utility functions -------
def _read_config_file(custom_input_file):
    '''
    Read-in data from user-created configuration file.
    
    Arguments
    ---------
    custom_input_file (str): configuration filepath, .yaml
    
    Returns
    -------
    bmi_parameters               (dict): Input parameters re bmi configuration
    log_parameters               (dict): Input parameters re logging
    compute_parameters           (dict): Input parameters re computation settings
    data_assimilation_parameters (dict): Input parameters re data assimilation

    '''
    with open(custom_input_file) as custom_file:
        data = yaml.load(custom_file, Loader=yaml.SafeLoader)

    troute_configuration = Config(**data)
    config_dict = troute_configuration.dict()

    compute_parameters = config_dict.get("compute_parameters")
    forcing_parameters = compute_parameters.get("forcing_parameters")
    data_assimilation_parameters = compute_parameters.get("data_assimilation_parameters")
    output_parameters = config_dict.get('output_parameters')

    return (
        compute_parameters,
        forcing_parameters,
        data_assimilation_parameters,
        output_parameters,
        )

def _read_timeslice_files(filepath, 
                          dates, 
                          qc_threshold, 
                          frequency_secs, 
                          cpu_pool=1,
                          interpolation_limit=59,
                          ):
    #Read files
    observation_df = pd.DataFrame()
    for d in dates:
        f = glob.glob(filepath + '/' + d + '*')

        if f:
            temp_df = xr.open_dataset(f[0])[['stationId','time','discharge','discharge_quality']].to_dataframe()
            observation_df = pd.concat([observation_df, temp_df])

    if not observation_df.empty:
        observation_df['stationId'] = observation_df['stationId'].str.decode('utf-8').str.strip()
        observation_df['time'] = observation_df['time'].str.decode('utf-8')
        observation_df['discharge_quality'] = observation_df['discharge_quality']/100

        #QC/QA and interpolation
        observation_df.loc[observation_df['discharge_quality']<0, 'discharge'] = np.nan
        observation_df.loc[observation_df['discharge_quality']>1, 'discharge'] = np.nan
        observation_df.loc[observation_df['discharge_quality']<qc_threshold, 'discharge'] = np.nan
        observation_df.loc[observation_df['discharge']<=0, 'discharge'] = np.nan

        observation_df = observation_df[['stationId','time','discharge']].set_index(['stationId', 'time']).unstack(1, fill_value = np.nan)['discharge']

        # ---- Interpolate USGS observations to the input frequency (frequency_secs)
        observation_df_T = observation_df.transpose()             # transpose, making time the index
        observation_df_T.index = pd.to_datetime(
            observation_df_T.index, format = "%Y-%m-%d_%H:%M:%S"  # index variable as type datetime
        )
        
        # specify resampling frequency 
        frequency = str(int(frequency_secs/60))+"min"    

        # interpolate and resample frequency
        buffer_df = observation_df_T.resample(frequency).asfreq()
        with Parallel(n_jobs=cpu_pool) as parallel:
            
            jobs = []
            interp_chunks = ()
            step = 200
            for a, i in enumerate(range(0, len(observation_df_T.columns), step)):
                
                start = i
                if (i+step-1) < buffer_df.shape[1]:
                    stop = i+(step)
                else:
                    stop = buffer_df.shape[1]
                    
                jobs.append(
                    delayed(_interpolate_one)(observation_df_T.iloc[:,start:stop], interpolation_limit, frequency)
                )
                
            interp_chunks = parallel(jobs)

        observation_df_T = pd.DataFrame(
            data = np.concatenate(interp_chunks, axis = 1), 
            columns = buffer_df.columns, 
            index = buffer_df.index
        )
        
        # re-transpose, making link the index
        observation_df_new = observation_df_T.transpose()
    
    else:
        observation_df_new = pd.DataFrame()

    return observation_df_new

def _interpolate_one(df, interpolation_limit, frequency):
    
    interp_out = (df.resample('min').
                        interpolate(
                            limit = interpolation_limit, 
                            limit_direction = 'both'
                        ).
                        resample(frequency).
                        asfreq().
                        to_numpy()
                       )
    return interp_out

def _read_timeseries_files(filepath, timeseries_dates, t0, final_persist_datetime):
    # Search for most recent RFC timseries file based on offset hours and lookback window
    # for each location.
    files = glob.glob(filepath + '/*')
    # create temporary dataframe with file names, split up by location and datetime
    df = pd.DataFrame([f.split('/')[-1].split('.') for f in files], columns=['Datetime','dt','ID','rfc','ext'])
    df = df[df['Datetime'].isin(timeseries_dates)][['ID','Datetime']]
    # For each location, find the most recent timeseries file (within timeseries window calculated a priori)
    df['Datetime'] = df['Datetime'].apply(lambda _: datetime.strptime(_, '%Y-%m-%d_%H'))
    df = df.groupby('ID').max().reset_index()
    df['Datetime'] = df['Datetime'].dt.strftime('%Y-%m-%d_%H')

    # Loop through list of timeseries files and store relevent information in dataframe.
    file_list = (df['Datetime'] + '.60min.' + df['ID'] + '.RFCTimeSeries.ncdf').tolist()
    rfc_df = pd.DataFrame()
    for f in file_list:
        ds = xr.open_dataset(filepath + '/' + f)
        sliceStartTime = datetime.strptime(ds.attrs.get('sliceStartTimeUTC'), '%Y-%m-%d_%H:%M:%S')
        sliceTimeResolutionMinutes = ds.attrs.get('sliceTimeResolutionMinutes')
        df = ds.to_dataframe().reset_index().sort_values('forecastInd')[['stationId','discharges','synthetic_values','totalCounts','timeSteps']]
        df['Datetime'] = pd.date_range(sliceStartTime, periods=df.shape[0], freq=sliceTimeResolutionMinutes+'T')
        # Filter out forecasts that go beyond the rfc_persist_days parameter. This isn't necessary, but removes
        # excess data, keeping the dataframe of observations as small as possible.
        df = df[df['Datetime']<final_persist_datetime]
        # Locate where t0 is in the timeseries
        df['timeseries_idx'] = df.index[df.Datetime == t0][0]
        df['file'] = f

        # Validate data to determine whether or not it will be used.
        use_rfc = _validate_RFC_data(
            df['stationId'][0],
            df.discharges,
            df.synthetic_values,
            filepath,
            f,
            300, #NOTE: this is t-route's default timestep. This will need to be verifiied again within t-route...
            False
        )
        df['use_rfc'] = use_rfc
        df['da_timestep'] = int(sliceTimeResolutionMinutes)*60

        rfc_df = pd.concat([rfc_df, df])
    rfc_df['stationId'] = rfc_df['stationId'].str.decode('utf-8').str.strip()
    return rfc_df

def _read_lastobs_file(
        lastobsfile,
        station_id = "stationId",
        ref_t_attr_id = "modelTimeAtOutput",
        time_idx_id = "timeInd",
        obs_discharge_id = "discharge",
        discharge_nan = -9999.0,
        time_shift = 0,
        ):
    with xr.open_dataset(lastobsfile) as ds:
        gages = ds[station_id].values
        
        ref_time = datetime.strptime(ds.attrs[ref_t_attr_id], "%Y-%m-%d_%H:%M:%S")
        
        last_ts = ds[time_idx_id].values[-1]
        
        df_discharge = (
            ds[obs_discharge_id].to_dataframe().                 # discharge to MultiIndex DF
            replace(to_replace = discharge_nan, value = np.nan). # replace null values with nan
            unstack(level = 0)                                   # unstack to single Index (timeInd)    
        )
        
        last_obs_index = (
            df_discharge.
            apply(pd.Series.last_valid_index).                   # index of last non-nan value, each gage
            to_numpy()                                           # to numpy array
        )
        last_obs_index = np.nan_to_num(last_obs_index, nan = last_ts).astype(int)
                        
        last_observations = []
        lastobs_times     = []
        for i, idx in enumerate(last_obs_index):
            last_observations.append(df_discharge.iloc[idx,i])
            lastobs_times.append(ds.time.values[i, idx].decode('utf-8'))
            
        last_observations = np.array(last_observations)
        lastobs_times     = pd.to_datetime(
            np.array(lastobs_times), 
            format="%Y-%m-%d_%H:%M:%S", 
            errors = 'coerce'
        )

        lastobs_times = (lastobs_times - ref_time).total_seconds()
        lastobs_times = lastobs_times - time_shift

    data_var_dict = {
        'gages'               : gages,
        'time_since_lastobs'  : lastobs_times,
        'lastobs_discharge'   : last_observations
    }

    lastobs_df = pd.DataFrame(data = data_var_dict)
    lastobs_df['gages'] = lastobs_df['gages'].str.decode('utf-8')

    return lastobs_df

def _read_lite_restart(file):
    '''
    Open lite restart pickle files. Can open either waterbody_restart or channel_restart
    
    Arguments
    -----------
        file (string): File path to lite restart file
        
    Returns
    ----------
        df (DataFrame): restart states
        t0 (datetime): restart datetime
    '''
    # open pickle file to pandas DataFrame
    df = pd.read_pickle(pathlib.Path(file))
    
    # extract restart time as datetime object
    t0 = df['time'].iloc[0].to_pydatetime()
    
    return df.drop(columns = 'time') , t0

def _write_lite_restart(
    values,
    t0,
    lite_restart
):
    '''
    Save initial conditions dataframes as pickle files
    
    Arguments
    -----------
        q0 (DataFrame):
        waterbodies_df (DataFrame):
        t0 (datetime.datetime):
        restart_parameters (string):
        
    Returns
    -----------
        
    '''
    
    output_directory = lite_restart.get('lite_restart_output_directory', None)
    if output_directory:
        
        # retrieve/reconstruct t-route produced output dataframes
        q0 = values['q0']
        q0_ids = values['q0_ids']
        q0 = pd.DataFrame(data=q0.reshape(len(q0_ids), -1), 
                          index=q0_ids,
                          columns=['qu0','qd0','h0'])

        waterbodies_df = values['waterbody_df']
        waterbodies_ids = values['waterbody_df_ids']
        
        waterbodies_df = pd.DataFrame(data=waterbodies_df.reshape(len(waterbodies_ids), -1), 
                                      index=waterbodies_ids,
                                      columns=['ifd','LkArea','LkMxE','OrificeA','OrificeC','OrificeE',
                                               'WeirC','WeirE','WeirL','qd0','h0','index'])
        waterbodies_df.index.name = 'lake_id'

        timestamp = t0 + timedelta(seconds=values['t-route_model_time'])

        # create restart filenames
        timestamp_str = timestamp.strftime("%Y%m%d%H%M")
        channel_restart_filename = 'channel_restart_' + timestamp_str
        waterbody_restart_filename = 'waterbody_restart_' + timestamp_str
        
        q0_out = q0.copy()
        q0_out['time'] = timestamp
        q0_out.to_pickle(pathlib.Path.joinpath(output_directory, channel_restart_filename))

        if not waterbodies_df.empty:
            wbody_initial_states = waterbodies_df.loc[:,['qd0','h0']]
            wbody_initial_states['time'] = timestamp
            wbody_initial_states.to_pickle(pathlib.Path.joinpath(output_directory, waterbody_restart_filename))

def _write_lastobs(
    values,
    t0,
    lastobs_output_folder=False,
):

    # join gageIDs to lastobs_df
    lastobs_df = values['lastobs_df']
    lastobs_df_ids = values['lastobs_df_ids']
    lastobs_df = pd.DataFrame(data=lastobs_df.reshape(len(lastobs_df_ids), -1),
                              index=lastobs_df_ids,
                              columns=['time_since_lastobs','lastobs_discharge','gages'])

    # timestamp of last simulation timestep
    modelTimeAtOutput = t0 + timedelta(seconds = values['t-route_model_time'])
    modelTimeAtOutput_str = modelTimeAtOutput.strftime('%Y-%m-%d_%H:%M:%S')

    # timestamp of last observation
    var = [timedelta(seconds=d) for d in lastobs_df.time_since_lastobs.fillna(0)]
    # shorvath (10/18/23): This code was copied from nhd_io.py based on V3 operations.
    # I changed 'modelTimeAtOutput - d' to 'modelTimeAtOutput + d' here because the '-' was
    # creating timestamps that are ahead of the model time which can't be right. I'm not sure
    # if this method was thoroughly checked for V3, but I've only made the update here.
    #TODO: Determine if this should also be changed in nhd_io.py...
    lastobs_timestamp = [modelTimeAtOutput + d for d in var]
    lastobs_timestamp_str = [d.strftime('%Y-%m-%d_%H:%M:%S') for d in lastobs_timestamp]
    lastobs_timestamp_str_array = np.asarray(lastobs_timestamp_str,dtype = '|S19').reshape(len(lastobs_timestamp_str),1)
    
    # create xarray Dataset similarly structured to WRF-generated lastobs netcdf files
    ds = xr.Dataset(
        {
            "stationId": (["stationIdInd"], lastobs_df["gages"].to_numpy(dtype = '|S15')),
            "time": (["stationIdInd", "timeInd"], np.asarray(lastobs_timestamp_str_array,dtype = '|S19')),
            "discharge": (["stationIdInd", "timeInd"], lastobs_df["lastobs_discharge"].to_numpy().reshape(len(lastobs_df["lastobs_discharge"]),1)),
        }
    )
    ds.attrs["modelTimeAtOutput"] = modelTimeAtOutput_str

    # write-out LastObs file as netcdf
    if isinstance(lastobs_output_folder, pathlib.Path):
        lastobs_output_folder = str(lastobs_output_folder)
    output_path = pathlib.Path(lastobs_output_folder + "/nudgingLastObs." + modelTimeAtOutput_str + ".nc").resolve()
    ds.to_netcdf(str(output_path))