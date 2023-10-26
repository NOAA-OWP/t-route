import numpy as np
import pandas as pd
import yaml
from datetime import datetime, timedelta
from joblib import delayed, Parallel
import xarray as xr
import glob

from troute.routing.fast_reach.reservoir_RFC_da import _validate_RFC_data

from nwm_routing.log_level_set import log_level_set
from troute.config import Config

class DAforcing_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_data_assimilation_parameters', '_forcing_parameters', '_compute_parameters',
                     '_usgs_df', 'reservoir_usgs_df', 'reservoir_usace_df', '_rfc_timeseries_df', '_lastobs_df' ]

        if bmi_cfg_file:
            (compute_parameters,
             forcing_parameters,
             data_assimilation_parameters,
            ) = _read_config_file(bmi_cfg_file)
            
            self._compute_parameters = compute_parameters
            self._forcing_parameters = forcing_parameters
            self._data_assimilation_parameters = data_assimilation_parameters

            print(data_assimilation_parameters)

            #############################
            # Read DA files:
            #############################
            nudging = data_assimilation_parameters.get('streamflow_da', {}).get('streamflow_nudging', False)
            
            usgs_persistence = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_persistence_da', {}).get('reservoir_persistence_usgs', False)
            usace_persistence = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_persistence_da', {}).get('reservoir_persistence_usace', False)
            rfc = data_assimilation_parameters.get('reservoir_da', {}).get('reservoir_rfc_da', {}).get('reservoir_rfc_forecasts', False)

            #usace_persistence = False

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

            #import pdb; pdb.set_trace()

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

            # to do: flatten the following arrays:
            #   Confirm that there is no usace dataframe (only a reservoir usace dataframe)? 
            #   self._reservoir_usace_df

            #import pdb; pdb.set_trace()

            # save time reference: start_datetime
            self.dateNull = start_datetime

            # read in metadata for BMI compliant arrays:

            # USGS Observations
            if nudging or usgs_persistence:

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
                    = _time_stations_from_df(self._usgs_df,start_datetime)
                # save metadata in class instances
                self.datesSecondsArray_usgs = _datesSecondsArray_usgs
                self.nDates_usgs = _nDates_usgs
                self.stationArray_usgs = _stationArray_usgs
                self.stationStringLengthArray_usgs = _stationStringLengthArray_usgs
                self.nStations_usgs = _nStations_usgs
                # flatten the actual USGS datafrane into a numpy ndarray
                _usgsArray = _flatten_array(self._usgs_df)
                # ... and save it with the class instance
                self.usgsArray = _usgsArray
         
                #import pdb; pdb.set_trace()

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
                    = _time_stations_from_df(self._reservoir_usgs_df,start_datetime)
                # save metadata in class instances
                self.datesSecondsArray_reservoir_usgs = _datesSecondsArray_reservoir_usgs
                self.nDates_reservoir_usgs = _nDates_reservoir_usgs
                self.stationArray_reservoir_usgs = _stationArray_reservoir_usgs
                self.stationStringLengthArray_reservoir_usgs = _stationStringLengthArray_reservoir_usgs
                self.nStations_reservoir_usgs = _nStations_reservoir_usgs
                # flatten the actual USGS datafrane into a numpy ndarray
                _reservoirUsgsArray = _flatten_array(self._reservoir_usgs_df)
                # ... and save it with the class instance
                self.reservoirUsgsArray = _reservoirUsgsArray            

            # USACE Observations        
            if usace_persistence:

                # see detailed comments in USGS branch
                ( _datesSecondsArray_reservoir_usace, _nDates_reservoir_usace, \
                    _stationArray_reservoir_usace, _stationStringLengthArray_reservoir_usace, \
                    _nStations_reservoir_usace) \
                    = _time_stations_from_df(self._reservoir_usace_df,start_datetime)
                # save metadata in class instances
                self.datesSecondsArray_reservoir_usace = _datesSecondsArray_reservoir_usace
                self.nDates_reservoir_usace = _nDates_reservoir_usace
                self.stationArray_reservoir_usace = _stationArray_reservoir_usace
                self.stationStringLengthArray_reservoir_usace = _stationStringLengthArray_reservoir_usace
                self.nStations_reservoir_usace = _nStations_reservoir_usace
                # flatten the actual USACE datafrane into a numpy ndarray
                _reservoirUsaceArray = _flatten_array(self._reservoir_usace_df)
                # ... and save it with the class instance
                self.reservoirUsaceArray = _reservoirUsaceArray  

                #import pdb; pdb.set_trace()

            testRevert = True
            # The following code segments revert the array troute -> bmi-flattened
            # array conversions, and are placed here only temporarily for
            # verification purposes, to be moved to DataAssimilation eventually.
            # If "reverse" code is not removed, it can be temporarily disabled
            # by setting the testRevert flag to False
            if (testRevert):

                # USGS dataframe
                
                #import pdb; pdb.set_trace()
                # Unflatten the arrays
                df_raw_usgs = _unflatten_array(self.usgsArray,self.nDates_usgs,\
                                          self.nStations_usgs)
                #import pdb; pdb.set_trace()
                # Decode time/date axis
                timeAxisName = 'time'
                freqString = '5T'
                df_withDates_usgs = _time_retrieve_from_arrays(df_raw_usgs, self.dateNull, \
                            self.datesSecondsArray_usgs, timeAxisName, freqString)
                #import pdb; pdb.set_trace()
                # Decode station ID axis
                stationAxisName = 'stationId'
                df_withStationsAndDates_usgs = _stations_retrieve_from_arrays(df_withDates_usgs, \
                    self.stationArray_usgs, self.stationStringLengthArray_usgs, \
                    stationAxisName)

                # Reservoir USGS dataframe

                #import pdb; pdb.set_trace()
                # Unflatten the arrays
                df_raw_reservoirUsgs = _unflatten_array(self.reservoirUsgsArray,self.nDates_reservoir_usgs,\
                                          self.nStations_reservoir_usgs)
                #import pdb; pdb.set_trace()
                # Decode time/date axis
                timeAxisName = 'time'
                freqString = '15T'
                df_withDates_reservoirUsgs = _time_retrieve_from_arrays(df_raw_reservoirUsgs, self.dateNull, \
                            self.datesSecondsArray_reservoir_usgs, timeAxisName, freqString)
                #import pdb; pdb.set_trace()
                # Decode station ID axis
                stationAxisName = 'stationId'
                df_withStationsAndDates_reservoirUsgs = _stations_retrieve_from_arrays(df_withDates_reservoirUsgs, \
                    self.stationArray_reservoir_usgs, self.stationStringLengthArray_reservoir_usgs, \
                    stationAxisName)
                
                # Reservoir USACE dataframe

                #import pdb; pdb.set_trace()
                # Unflatten the arrays
                df_raw_reservoirUsace = _unflatten_array(self.reservoirUsaceArray,self.nDates_reservoir_usace,\
                                          self.nStations_reservoir_usace)
                #import pdb; pdb.set_trace()
                # Decode time/date axis
                timeAxisName = 'time'
                freqString = '15T'
                df_withDates_reservoirUsace = _time_retrieve_from_arrays(df_raw_reservoirUsace, self.dateNull, \
                            self.datesSecondsArray_reservoir_usace, timeAxisName, freqString)
                #import pdb; pdb.set_trace()
                # Decode station ID axis
                stationAxisName = 'stationId'
                df_withStationsAndDates_reservoirUsace = _stations_retrieve_from_arrays(df_withDates_reservoirUsace, \
                    self.stationArray_reservoir_usace, self.stationStringLengthArray_reservoir_usace, \
                    stationAxisName)                

                import pdb; pdb.set_trace()

        else:
            raise(RuntimeError("No config file provided."))


def _time_stations_from_df(dataFrame, timeBase):

    # get columns (date-time), extract as list at first
    datesList = (dataFrame.columns).tolist()
    # then subtract timebase
    datesListSeconds = [int((d-timeBase).total_seconds()) for d in datesList]
    nDates = len(datesListSeconds)
    datesSecondsArray = np.array(datesListSeconds)

    # get station IDs, as list
    stations = (dataFrame).index.tolist()
    nStations = len(stations)

    # build station array (one giant 1D ndarray with ASCII encodings
    # concatenated together), aided by an array of how many strings each
    # station ID consists of - completely USGS, USACE, etc agnostic 
    #
    # define array defining string length for each station ID
    stationStringLengthArray = np.empty(nStations, dtype=np.int64)
    stationStringLengthArray[0] = len(stations[0])

    # start big array with ID for first station
    stationArray = np.empty(stationStringLengthArray[0], dtype=np.int64)

    # init some aux. variables
    firstIter = True
    stationIndex = 0
    # iterate through all stations
    for station in stations:
        # capture string length of present station
        stationStringLengthArray[stationIndex] = len(station)
        stationIndex += 1        
        # convert string to ASCII code
        charArray = np.array(station, dtype='c')
        addArray = np.array([ord(ch) for ch in charArray])        
        # build up big array
        if (firstIter):
            stationArray = addArray
            firstIter = False
        else:
            stationArray = np.concatenate([stationArray, addArray])

    return datesSecondsArray, nDates, stationArray, stationStringLengthArray, nStations


def _flatten_array(dataFrame):

    # convert to numpy array first
    array1 = dataFrame.to_numpy(copy=True, dtype=np.float16)
    # flatten it
    array2 = array1.reshape(-1)

    return array2


def _unflatten_array(array_1D,nx,ny):

    # reshape numpy array
    array_2D = array_1D.reshape(ny,nx)

    # convert to data frame
    df_raw = pd.DataFrame(array_2D)

    return df_raw


def _time_retrieve_from_arrays(dataFrame, timeBase, datesSecondsArrays, \
                               timeAxisName, freqString):

    # convert array to list at first
    datesList = (datesSecondsArrays).tolist()

    # add timebase back in, and convert to datetime format
    datesList = [(timedelta(seconds=d)+timeBase) for d in datesList]

    # convert to DatetimeIndex
    dateTimeIndex = pd.DatetimeIndex(datesList, name=timeAxisName, freq=freqString)

    # transpose dataframe
    dataFrameTranspose = dataFrame.T
    # add time axis as index on transposed dataframe (time is in rows)
    dataFrameTranspose = dataFrameTranspose.set_index(dateTimeIndex)
    # revert to original (time is in columns)
    dataFrameWithTime = dataFrameTranspose.T

    return dataFrameWithTime


def _stations_retrieve_from_arrays(dataFrame, stationsArray, stationStringLengthArray, \
                               stationAxisName):

    # first create station outputs as list
    stationsList = []

    # mark place in global stationsArray
    i=0

    for stringLength in stationStringLengthArray:

        # extract small array with ASCII codes of one station ID
        stationID_asArray = stationsArray[i:i+stringLength]
        # increment global pointer along array
        i=i+stringLength
        # convert to list of characters
        firstRun = True
        for asciiCode in stationID_asArray:
            # decode Ascii back to character and stick it to the end of the list
            chrAdd = chr(asciiCode)
            if firstRun == True:
                stationID_asString = chrAdd
                firstRun = False
            else:
                stationID_asString += chrAdd

        #import pdb; pdb.set_trace()

        # done with one station ID
        stationsList.append(stationID_asString)

    # construct index    
    index = pd.Index(stationsList, name=stationAxisName, dtype=object)

    dataFrame.index = (index)

    #import pdb; pdb.set_trace()

    return dataFrame



        


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

    return (
        compute_parameters,
        forcing_parameters,
        data_assimilation_parameters,
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

        print('date ',d)
        print('file ',filepath + '/' + d + '*')

        if f:
            temp_df = xr.open_dataset(f[0])[['stationId','time','discharge','discharge_quality']].to_dataframe()
            observation_df = pd.concat([observation_df, temp_df])
    

    print('Dates ',dates)
    print('Observation station ID ',observation_df)

    #import pdb; pdb.set_trace()

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
    
    print('Enter interpolate and resample')

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