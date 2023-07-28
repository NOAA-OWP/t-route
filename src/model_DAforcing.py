import numpy as np
import netCDF4 as nc
import pandas as pd
import yaml
from array import array
import datetime
import pathlib
import logging
from joblib import delayed, Parallel
import netCDF4
import xarray as xr


from troute.routing.fast_reach.reservoir_hybrid_da import reservoir_hybrid_da
from troute.routing.fast_reach.reservoir_RFC_da import reservoir_RFC_da, preprocess_RFC_data
from troute.DataAssimilation import build_lastobs_df, _create_usgs_df, _reindex_link_to_lake_id
#import troute.nhd_network_utilities_v02 as nnu


class DAforcing_model():

    def __init__(self, bmi_cfg_file=None):
        """
        
        """
        __slots__ = ['_data_assimilation_parameters', '_streamflow_da_parameters',
                     '_lake_number', '_res_type', '_gage_lakeID_crosswalk_file', 
                     '_time', '_time_step', 
                     '_rfc_timeseries_idx', '_rfc_gage_id', '_rfc_timeseries_folder',
                     '_rfc_timeseries_file', '_rfc_timeseries_offset_hours', '_rfc_forecast_persist_days',
                     '_time_step_seconds', '_total_counts', '_rfc_timeseries_discharges', '_rfc_update_time', 
                     '_use_RFC',
                     '_timeslice_obs', '_timeslice_datetime', '_timeslice_stationId',
                     '_timeslice_discharge_quality', 
                     '_nudging', '_lastobs_stationId', '_time_since_lastobs', '_lastobs_discharge' ]

        if bmi_cfg_file:
            (#bmi_parameters, #TODO We might not need any bmi specific parameters
             log_parameters, #TODO Update all logging/warnings throughout standalone reservoirs...
             compute_parameters,
             restart_parameters,
             forcing_parameters,
             data_assimilation_parameters,
            ) = _read_config_file(bmi_cfg_file)

            self._compute_parameters = compute_parameters
            self._forcing_parameters = forcing_parameters
            self._data_assimilation_parameters = data_assimilation_parameters

            #if compute_parameters:
            self._time_step = forcing_parameters.get("dt", 300) #time step in sec
            self._nts = forcing_parameters.get("nts", 288) #the number of time steps
            #NOTE: For any each operation either of AnA or Forecast, t0 (model start date) is
            #provided from config file.
            self._t0 = restart_parameters.get("model_start_time", None)   
            rfc_parameters = data_assimilation_parameters.get("rfc_reservoir_da", None) 
            
            if not self._t0:
                raise(RuntimeError("No start_datetime provided in config file."))

            if rfc_parameters:
                #self._rfc_gage_id = rfc_parameters.get("reservoir_rfc_gage_id", None)
                
                #self._reservoir_rfc_forecast_lookback_hours= rfc_parameters.get("reservoir_rfc_forecast_lookback_hours",None)
                self._rfc_timeseries_folder = rfc_parameters.get("rfc_forecast_time_series_folder", None)
                self._gage_lakeID_crosswalk_file = rfc_parameters.get("gage_lakeID_crosswalk_file", None)
                self._rfc_forecast_peyrsist_days = rfc_parameters.get("reservoir_rfc_forecast_persist_days", 11)
                self._rfc_timeseries_offset_hours = rfc_parameters.get("reservoir_rfc_timeseries_offset_hours", 28)
        else:
            raise(RuntimeError("No config file provided."))
            
        self._time = 0.0  # initiall set zero but will be in sync with troute_model's self._time (sec)
    
    def preprocess_static_vars(self, values: dict):
        '''
        Create static data structures for reservoir DA objects and streamflow nudging 

        Parameters
        ----------

            
        Returns
        -------

        '''
        data_assimilation_parameters = self._data_assimilation_parameters
        streamflow_da_parameters     = data_assimilation_parameters.get('streamflow_da', None)

        if streamflow_da_parameters:
            self._nudging = streamflow_da_parameters.get('streamflow_nudging', False)
        
        self._lastobs_file = streamflow_da_parameters.get("wrf_hydro_lastobs_file", None)
        self._lastobs_crosswalk_file = streamflow_da_parameters.get("gage_segID_crosswalk_file", None)
        self._lastobs_start = streamflow_da_parameters.get("wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0)
        self._streamflow_da_parameters = streamflow_da_parameters 

        self._lake_number = values['waterbody__lake_number'] # NHDWaterbodyComID in RouteLink.nc
        self._res_type = values['waterbody__type_number']

        if self._res_type==4 or self._res_type==5:
        # prepare static data for RFC DA
            self.gage_lake_crosswalk_beta()

    def gage_lake_crosswalk_beta(self):
        '''
        Compute gage to lake crosswalk DataFrame using lake_gage_df.csv

        Parameters
        ----------
        values: dict
            The model state data structure.
            
        Returns
        -------
        rfc_gage_id
        '''
        if self._res_type==4 or self._res_type==5:
            #RFC reservoir DA
            file= self._gage_lakeID_crosswalk_file
            lake_gage_df = pd.read_csv(file)
            self._rfc_gage_id = lake_gage_df.loc[(lake_gage_df['lake_id']==int(self._lake_number[0]))&
                                            (lake_gage_df['reservoir_type']==self._res_type[0]), 
                                            'gage_id'].values[0]

    '''
    def gage_lake_crosswalk(self, values:dict,):

        Compute gage to lake crosswalk DataFrame using reservoir_index_AnA.nc

        Parameters
        ----------
        values: dict
            The model state data structure.
            
        Returns
        -------

      
        dataset = nc.Dataset(file)
        lake_id = pd.Series(dataset.variables['lake_id'][:])
        reservoir_type = pd.Series(dataset.variables['reservoir_type'][:])
        usgs_lake_id = pd.Series(dataset.variables['usgs_lake_id'][:])
        usace_lake_id = pd.Series(dataset.variables['usace_lake_id'][:])
        rfc_lake_id = pd.Series(dataset.variables['rfc_lake_id'][:])

        # Get rfc_gage_id correponding to rfc_lake_id
        rfc_gage_id_len = dataset.variables['rfc_gage_id'].shape[0]
        masked_array  = dataset.variables['rfc_gage_id'][:,:]
        rfc_gage_id =[]
        for idx in range(rfc_gage_id_len):
            arr = masked_array[idx]
            #import pdb; pdb.set_trace()
            strings = arr.data.astype(str)
            gage_id=""
            for string in strings:
                gage_id=gage_id + string
            rfc_gage_id.append(gage_id)
        rfc_gage_id = pd.Series(rfc_gage_id)        
        
        # Get usgs_gage_id corresponding to usgs_lake_id
        usgs_gage_id_len = dataset.variables['usgs_gage_id'].shape[0]
        masked_array     = dataset.variables['usgs_gage_id'][:,:]
        usgs_gage_id     = []
        for idx in range(usgs_gage_id_len):
            arr     = masked_array[idx]
            strings = arr.data.astype(str)
            gage_id = ""
            for string in strings:
                gage_id = gage_id + string
            usgs_gage_id.append(gage_id)
        usgs_gage_id = pd.Series(usgs_gage_id)  

        # Get usace_gage_id corresponding to usace_lake_id
        usace_gage_id_len = dataset.variables['usace_gage_id'].shape[0]
        masked_array      = dataset.variables['usace_gage_id'][:,:]
        usace_gage_id      = []
        for idx in range(usace_gage_id_len):
            arr     = masked_array[idx]
            strings = arr.data.astype(str)
            gage_id = ""
            for string in strings:
                gage_id = gage_id + string
            usace_gage_id.append(gage_id)
        usace_gage_id = pd.Series(usace_gage_id) 

        series_dict={
                    'lake_id': lake_id,
                    'reservoir_type': reservoir_type,
                    'usgs_lake_id': usgs_lake_id,
                    'usgs_gage_id': usgs_gage_id,
                    'usace_lake_id': usace_lake_id,
                    'usace_gage_id': usace_gage_id,
                    'rfc_lake_id': rfc_lake_id,
                    'rfc_gage_id': rfc_gage_id,
                    }        

        # 
    '''

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

         # compute data assimilation time series and related dynamic parameters
        if self._res_type==2 or self._res_type==3 or self._nudging:                

            #Create 1d numpy array of last_obs_df with with each element starting with stationID, 
            #followed by time_since_lastobs and lastobs_discharge 
            last_obs_df       = pd.DataFrame()            
            lastobs_file = self._streamflow_da_parameters.get("wrf_hydro_lastobs_file", None)

            if lastobs_file:                
                #lastobs_crosswalk_file = streamflow_da_parameters.get("gage_segID_crosswalk_file", None)
                lastobs_start = self._streamflow_da_parameters.get("wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0)
       
                lastobs_df = build_lastobs_df(
                                                lastobs_file,
                                                #self._lastobs_crosswalk_file,
                                                lastobs_start,
                                                )
                #1D array with each element starting with stationID, followed by time_since_lastobs and lastobs_discharge 
                lastobs_arr =  np.array(lastobs_df.to_records(index=True))
                self._lastobs_stationId  = lastobs_arr['gages'].tolist()
                self._time_since_lastobs =  lastobs_arr['time_since_lastobs'].tolist()                
                self._lastobs_discharge  =  lastobs_arr['lastobs_discharge'].tolist()
                import pdb; pdb.set_trace()
            # replace link ids with lake ids, for gages at waterbody outlets, 
            # otherwise, gage data will not be assimilated at waterbody outlet
            # segments.
            #if network.link_lake_crosswalk:
            #    self._last_obs_df = _reindex_link_to_lake_id(self._last_obs_df, network.link_lake_crosswalk)
       
            #Create 1d numpy arrays with each element having station id and discharge timeseries or 
            #stationID and discharge qual. timeseries
            timeslice_obs_df  = pd.DataFrame()
            timeslice_qual_df = pd.DataFrame()
            # Create da_sets: sets of TimeSlice files for each loop          
            dt = datetime.datetime.strptime(self._t0,"%Y-%m-%d_%H:%M:%S") 
            formatted_date = datetime.datetime(dt.year, dt.month, dt.day, dt.hour, dt.minute)
            hours = self._time_step*self._nts/(60*60)
            td = datetime.timedelta(hours=hours)
            final_timestamp = formatted_date + td         
 
            if self._data_assimilation_parameters:
                #da_sets = nnu.build_da_sets(data_assimilation_parameters, run_sets, network.t0)
                da_sets = build_da_sets(self._data_assimilation_parameters, final_timestamp, self._t0)
            da_run=da_sets[0]
                       
            run_parameters = {
                            'dt': self._forcing_parameters.get('dt'),
                            'nts': self._forcing_parameters.get('nts'),
                            'cpu_pool': self._compute_parameters.get('cpu_pool')
                            }
            #from_files=True,c
            #value_dict=None, 
            persistence_reservoir_da = self._data_assimilation_parameters.get("persistence_reservoir_da", None)
            streamflow_da_parameters = self._data_assimilation_parameters.get('streamflow_da', None)
            
            # Compared to _create_usgs_df used in DataAssimilation.py, this _create_timeslice_df doesn't
            # do crosswalk b/t gageID and linkID and QC, all of which will be performed within t-route.
            timeslice_obs_df = _create_timeslice_df(
                                            persistence_reservoir_da,
                                            #self._data_assimilation_parameters, 
                                            streamflow_da_parameters, 
                                            run_parameters, 
                                            self._t0,
                                            #network, 
                                            da_run,
                                            )
            # convert a dataframe to a 1D numpy array while including the index value (station id) as part of the array
            timeslice_obs_arr      = np.array(timeslice_obs_df.to_records(index=True)) # 1d array with each element having station id and discharge timeseries
            n_elements_each_record = len(timeslice_obs_arr.dtype) - 1
            stationId              = timeslice_obs_arr['link']
            stationId_repeat_list  = (np.repeat(stationId, n_elements_each_record)).tolist()
            # repeated timestamps by n_element_each_record for 
            timestamps             = timeslice_obs_df.columns.strftime('%Y-%m-%d %H:%M:%S').tolist()
            datetime_repeat_list   = timestamps*len(stationId)
            # observed time series for one stationId after another 
            obs_list=[]
            for idx in range(len(stationId)):
                obs_at_station = timeslice_obs_arr[idx]
                temp_list = obs_at_station.tolist()
                obs_list.extend(temp_list[1:])

            self._timeslice_obs               = obs_list   
            self._timeslice_datetime          = datetime_repeat_list
            self._timeslice_stationId         = stationId_repeat_list
            self._timeslice_discharge_quality = np.full(len(obs_list),100)
            import pdb; pdb.set_trace()
        if self._res_type==4 or self._res_type==5:
            (self._use_RFC, 
             self._rfc_timeseries_discharges, 
             self._rfc_timeseries_idx, 
             self._rfc_update_time, 
             self._da_time_step, 
             self._total_counts,
             self._rfc_timeseries_file) = preprocess_RFC_data( 
                                                        self._t0,
                                                        self._rfc_timeseries_offset_hours,
                                                        self._rfc_gage_id,
                                                        self._rfc_timeseries_folder,
                                                        int(self._lake_number[0]),
                                                        self._time_step
                                                        )
        import pdb; pdb.set_trace()

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

    #bmi_parameters = data.get("bmi_parameters", None)
    log_parameters = data.get("log_parameters", None)
    compute_parameters = data.get("compute_parameters", None)
    restart_parameters = compute_parameters.get("restart_parameters", None)
    forcing_parameters = compute_parameters.get("forcing_parameters", None)
    data_assimilation_parameters = compute_parameters.get("data_assimilation_parameters", None)
    rfc_parameters = data_assimilation_parameters.get("rfc_reservoir_da", None)

    return (
        #bmi_parameters,
        log_parameters,
        compute_parameters,
        restart_parameters,
        forcing_parameters,
        data_assimilation_parameters,
        )

def build_da_sets(da_params, final_timestamp, t0):
    """
    Create set lists of USGS and/or USACE TimeSlice files used for 
    streamflow and reservoir data assimilation
        
    Arguments
    --------
    - da_params (dict): user-input data assimilation parameters
    - run_sets (list) : forcing files for each run set in the simlation
    - t0 (datetime)   : model initialization time
        
    Returns
    -------
    - da_sets (list)  : lists of USGS and USACE TimeSlice files for each run set 
                         in the simulation
        
    Notes
    -----
        
    """

    persistence_reservoir_da = da_params.get('persistence_reservoir_da', False)
    # check for user-input usgs and usace timeslice directories
    usgs_timeslices_folder = persistence_reservoir_da.get(
        "usgs_timeslices_folder",
        None
    )
    usace_timeslices_folder = persistence_reservoir_da.get(
        "usace_timeslices_folder",
        None
    )

    # User-specified DA ON/OFF preferences
    usace_da = False
    usgs_da = False
    #reservoir_da = da_params.get('reservoir_da', False)
    if persistence_reservoir_da:
        usgs_da = persistence_reservoir_da.get('reservoir_persistence_usgs', False)
        usace_da = persistence_reservoir_da.get('reservoir_persistence_usace', False)
        
    nudging = False
    streamflow_da = da_params.get('streamflow_da', False)
    if streamflow_da:
        nudging = streamflow_da.get('streamflow_nudging', False)
            
    if not usgs_da and not usace_da and not nudging:
        # if all DA capabilities are OFF, return empty dictionary
        #da_sets = [{} for _ in run_sets]
        da_sets= {}

    # if no user-input timeslice folders, a list of empty dictionaries
    elif not usgs_timeslices_folder and not usace_timeslices_folder:
        # if no timeslice folders, return empty dictionary
        #da_sets = [{} for _ in run_sets]
        da_sets = {}

    # if user-input timeslice folders are present, build TimeSlice sets
    else:            
        # create Path objects for each TimeSlice directory
        if usgs_timeslices_folder:
            usgs_timeslices_folder = pathlib.Path(usgs_timeslices_folder)
        if usace_timeslices_folder:
            usace_timeslices_folder = pathlib.Path(usace_timeslices_folder)
            
        # the number of timeslice files appended to the front- and back-ends
        # of the TimeSlice file interpolation stack
        pad_hours = persistence_reservoir_da.get("timeslice_lookback_hours",1)
        timeslice_pad = pad_hours * 4 # number of 15 minute TimeSlices in the pad

        # timedelta of TimeSlice data - typically 15 minutes
        dt_timeslice = datetime.timedelta(minutes = 15)

        da_sets = [] # initialize list to store TimeSlice set lists

        # Loop through each run set and build lists of available TimeSlice files
        if 1==1:
        #for (i, set_dict) in enumerate(run_sets):
            i=0    
            # Append an empty dictionary to the loop, which be used to hold
            # lists of USGS and USACE TimeSlice files.
            da_sets.append({})

            # timestamps of TimeSlice files desired for run set i
            t0 = datetime.datetime.strptime(t0,"%Y-%m-%d_%H:%M:%S") 
            timestamps = pd.date_range(
                t0 - dt_timeslice * timeslice_pad,
                #run_sets[i]['final_timestamp'] + dt_timeslice * 4,
                final_timestamp + dt_timeslice * 4,
                freq=dt_timeslice
            )

            # identify available USGS TimeSlices in run set i
            if (usgs_timeslices_folder and nudging) or (usgs_timeslices_folder and usgs_da):
                filenames_usgs = (timestamps.strftime('%Y-%m-%d_%H:%M:%S') 
                            + '.15min.usgsTimeSlice.ncdf').to_list()
                    
                # identify available USGS TimeSlices
                filenames_usgs = _check_timeslice_exists(
                        filenames_usgs, 
                        usgs_timeslices_folder
                )
                    
                # Add available TimeSlices to da_sets list
                da_sets[i]['usgs_timeslice_files'] = filenames_usgs
                    
            # identify available USACE TimeSlices in run set i
            if usace_timeslices_folder and usace_da:
                filenames_usace = (timestamps.strftime('%Y-%m-%d_%H:%M:%S') 
                                + '.15min.usaceTimeSlice.ncdf').to_list()
                    
                # identify available USACE TimeSlices
                filenames_usace = _check_timeslice_exists(
                        filenames_usace, 
                        usace_timeslices_folder
                )
                    
                # Add available TimeSlices to da_sets list
                da_sets[i]['usace_timeslice_files'] = filenames_usace

            # reset initialization time for loop set i+1
            #t0 = run_sets[i]['final_timestamp']
            t0 = final_timestamp
                
    return da_sets

def _check_timeslice_exists(filenames, timeslices_folder):
    """
    Check that each TimeSlice file in a list of files exists. Return a list of 
    available files.
    
    Arguments
    ---------
    - filenames               (list of str): TimeSlice filenames
    - timeslices_folder (pathlib.PosixPath): TimeSlice directory
    
    Returns
    -------
    - filenames_existing (list of chr): Existing TimeSlice filenames in 
                                        TimeSlice directory
    
    Notes
    -----
    - This is a utility function used by build_da_sets
    
    """
    LOG = logging.getLogger('')
    
    # check that all USGS TimeSlice files in the set actually exist
    drop_list = []
    for f in filenames:
        try:
            J = pathlib.Path(timeslices_folder.joinpath(f))     
            assert J.is_file() == True
        except AssertionError:
            LOG.warning("Missing TimeSlice file %s", J)
            drop_list.append(f)

    # Assemble a list of existing TimeSlice files, only
    filenames_existing = [x for x in filenames if x not in drop_list]   
    
    return filenames_existing


def _create_timeslice_df(
                    persistence_reservoir_da,
                    #data_assimilation_parameters, 
                    streamflow_da_parameters,                    
                    run_parameters, 
                    t0,
                    #network, 
                    da_run):
    '''
    Function for reading USGS or USACE timeslice files and creating a dataframe
    of USGS or USACE gage observations. This dataframe is used for streamflow
    nudging and can be used for constructing USGS/USACE persistence reservoir DA.
    
    Arguments:
    ----------
    - persistence_reservoir_da     (dict): user input data re persistence reservoir da
    - run_parameters               (dict): user input data re subset of compute configuration
    - da_run                       (list): list of data assimilation files separated by for loop chunks
    
    Returns:
    --------
    - usgs_df (DataFrame): dataframe of USGS gage observations
    '''

    usgs_timeslices_folder  = persistence_reservoir_da.get("usgs_timeslices_folder", None)
    usace_timeslices_folder = persistence_reservoir_da.get("usace_timeslices_folder", None)
    #crosswalk_file         = streamflow_da_parameters.get("gage_segID_crosswalk_file", None)
    crosswalk_gage_field    = streamflow_da_parameters.get('crosswalk_gage_field','gages')
    crosswalk_segID_field   = streamflow_da_parameters.get('crosswalk_segID_field','link')
    #da_decay_coefficient   = data_assimilation_parameters.get("da_decay_coefficient",120)
    qc_threshold            = persistence_reservoir_da.get("qc_threshold",1)
    interpolation_limit     = persistence_reservoir_da.get("interpolation_limit_min",59)


    # User-specified DA ON/OFF preferences
    usace_da = False
    usgs_da = False
    usgs_da = persistence_reservoir_da.get('reservoir_persistence_usgs', False)
    usace_da = persistence_reservoir_da.get('reservoir_persistence_usace', False)
    usgs_timeslices_folder = persistence_reservoir_da.get('usgs_timeslices_folder', None)
    usace_timeslices_folder = persistence_reservoir_da.get('usace_timeslices_folder', None)    

    if usgs_timeslices_folder and usgs_da:
        usgs_timeslices_folder = pathlib.Path(usgs_timeslices_folder)
        timeslice_files = [usgs_timeslices_folder.joinpath(f) for f in 
                        da_run['usgs_timeslice_files']]

    if usace_timeslices_folder and usace_da:    
        usace_timeslices_folder = pathlib.Path(usace_timeslices_folder)
        timeslice_files = [usace_timeslices_folder.joinpath(f) for f in 
                  da_run['usace_timeslice_files']]
	
    #t0 = datetime.datetime.strptime(t0,"%Y-%m-%d_%H:%M:%S") 

    if timeslice_files:
        timeslice_obs_df = get_obs_from_timeslices(
                #network.link_gage_df,
                crosswalk_gage_field,
                crosswalk_segID_field,
                timeslice_files,
                qc_threshold,
                interpolation_limit,
                run_parameters.get("dt"),
                t0,
                #network.t0,
                run_parameters.get("cpu_pool", None)
                ) #.loc[network.link_gage_df.index]
        
		
		# replace link ids with lake ids, for gages at waterbody outlets, 
		# otherwise, gage data will not be assimilated at waterbody outlet
		# segments.
        #if network.link_lake_crosswalk:
        #    usgs_df = _reindex_link_to_lake_id(usgs_df, network.link_lake_crosswalk)
    
    else:
        timeslice_obs_df  = pd.DataFrame()
    
    return timeslice_obs_df

def get_obs_from_timeslices(
    #crosswalk_df,
    crosswalk_gage_field,
    crosswalk_dest_field,
    timeslice_files,
    qc_threshold,
    interpolation_limit,
    frequency_secs,
    t0,
    cpu_pool, 
):
    """
    Read observations from TimeSlice files, interpolate available observations
    and organize into a Pandas DataFrame
    
    Aguments
    --------                                          
    - crosswalk_gage_field         (str): fieldname of gage ID data in crosswalk dataframe
    
    - crosswalk_dest_field         (str): fieldname of destination data in link_gage_df. 
                                          For streamflow DA, this is the field
                                          containing segment IDs. For reservoir DA, 
                                          this is the field containing waterbody IDs.
                                          
    - timeslice_files (list of PosixPath): Full paths to existing TimeSlice files
    
    - qc_threshold                  (int): Numerical observation quality 
                                           threshold. Observations with quality
                                           flags less than this value will be 
                                           removed and replaced witn nan.
                                           
    - interpolation_limit           (int): Maximum gap duration (minuts) over 
                                           which observations may be interpolated 
                                           
    - t0                       (datetime): Initialization time of simulation set
    
    - cpu_pool                      (int): Number of CPUs used for parallel 
                                           TimeSlice reading and interolation
    
    Returns
    -------
    - observation_df_new (Pandas DataFrame): 
    
    Notes
    -----
    The max_fill is applied when the series is being considered at a 1 minute interval
    so 14 minutes ensures no over-interpolation with 15-minute gage records, but creates
    square-wave signals at gages reporting hourly...
    therefore, we advise a 59 minute gap filling tolerance.
    
    """
    # TODO: 
    # - Parallelize this reading for speedup using netCDF4
    # - Generecize the function to read USACE or USGS data
    # - only return gages that are in the model domain (consider mask application)

    # open TimeSlce files, organize data into dataframes
    with Parallel(n_jobs=cpu_pool) as parallel:
        jobs = []
        for f in timeslice_files:
            jobs.append(delayed(_read_timeslice_file)(f))
        timeslice_dataframes = parallel(jobs)

    # create lists of observations and obs quality dataframes returned 
    # from _read_timeslice_file
    timeslice_obs_frames = []
    timeslice_qual_frames = []
    for d in timeslice_dataframes:
        timeslice_obs_frames.append(d[0])  # TimeSlice gage observation
        timeslice_qual_frames.append(d[1]) # TimeSlice observation qual

    # concatenate dataframes
    timeslice_obs_df  = pd.concat(timeslice_obs_frames, axis = 1)
    timeslice_qual_df = pd.concat(timeslice_qual_frames, axis = 1) 

    # While it should be Link <> gage crosswalk data, we use here a fake crosswalk_df
    # that has gage <> gage instead.
    fake_crosswalk_df = (pd.DataFrame({
                                crosswalk_dest_field: timeslice_obs_df.index,
                                crosswalk_gage_field:  timeslice_obs_df.index    
                        }))

    df = fake_crosswalk_df
    df[crosswalk_gage_field] = np.asarray(df[crosswalk_gage_field]).astype('<U15')
    df = df.set_index(crosswalk_gage_field)
    
    # join crosswalk data with timeslice data, indexed on crosswalk destination field
    observation_df = (df.join(timeslice_obs_df).
               reset_index().
               set_index(crosswalk_dest_field).
               drop([crosswalk_gage_field], axis=1))

    observation_qual_df = (df.join(timeslice_qual_df).
               reset_index().
               set_index(crosswalk_dest_field).
               drop([crosswalk_gage_field], axis=1))

    # ---- Laugh testing ------
    # screen-out erroneous qc flags
    observation_qual_df = (observation_qual_df.
                           mask(observation_qual_df < 0, np.nan).
                           mask(observation_qual_df > 1, np.nan)
                          )

    # screen-out poor quality flow observations
    observation_df = (observation_df.
                      mask(observation_qual_df < qc_threshold, np.nan).
                      mask(observation_df <= 0, np.nan)
                     )

    # ---- Interpolate USGS observations to the input frequency (frequency_secs)
    # This step is especially essential as actual published time of flow for stations
    # can be different from timestamp included in the names of timeslice files.
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
    # TODO: Witine t-route or reservoir module, link value should be replaced by actual linkId 
    # correponding to stationId (link here is actually still stationId)
    observation_df_new = observation_df_T.transpose()    

    return observation_df_new
 

def _read_timeslice_file(f):

    with netCDF4.Dataset(
        filename = f,
        mode = 'r',
        format = "NETCDF4"
    ) as ds:
        
        discharge = ds.variables['discharge'][:].filled(fill_value = np.nan)
        stns      = ds.variables['stationId'][:].filled(fill_value = np.nan)
        t         = ds.variables['time'][:].filled(fill_value = np.nan)
        qual      = ds.variables['discharge_quality'][:].filled(fill_value = np.nan)
        
    stationId = np.apply_along_axis(''.join, 1, stns.astype(np.str))
    time_str = np.apply_along_axis(''.join, 1, t.astype(np.str))
    stationId = np.char.strip(stationId)

    timeslice_observations = (pd.DataFrame({
                                'stationId' : stationId,
                                'datetime'  : time_str,
                                'discharge' : discharge
                            }).
                             set_index(['stationId', 'datetime']).
                             unstack(1, fill_value = np.nan)['discharge'])
    
    observation_quality = (pd.DataFrame({
                                'stationId' : stationId,
                                'datetime'  : time_str,
                                'quality'   : qual/100
                            }).
                             set_index(['stationId', 'datetime']).
                             unstack(1, fill_value = np.nan)['quality'])
    
    return timeslice_observations, observation_quality

def build_lastobs_df(
        lastobsfile,
        #crosswalk_file,
        time_shift           = 0,
        crosswalk_gage_field = "gages",
        crosswalk_link_field = "link",
        obs_discharge_id     = "discharge",
        time_idx_id          = "timeInd",
        station_id           = "stationId",
        station_idx_id       = "stationIdInd",
        time_id              = "time",
        discharge_nan        = -9999.0,
        ref_t_attr_id        = "modelTimeAtOutput",
        route_link_idx       = "feature_id",
    ):
    '''
    Constructs a DataFame of "lastobs" data used in streamflow DA routine
    "lastobs" information is just like it sounds. It is the magnitude and
    timing of the last valid observation at each gage in the model domain. 
    We use this information to jump start initialize the DA process, both 
    for forecast and AnA simulations. 
    
    Arguments
    ---------
    
    Returns
    -------
    
    Notes
    -----
    
    '''
    # open crosswalking file and construct dataframe relating gageID to segmentID
    #with xr.open_dataset(crosswalk_file) as ds:
    #    gage_list = list(map(bytes.strip, ds[crosswalk_gage_field].values))
    #    gage_mask = list(map(bytes.isalnum, gage_list))
    #    gage_da   = list(map(bytes.strip, ds[crosswalk_gage_field][gage_mask].values))
    #    data_var_dict = {
    #        crosswalk_gage_field: gage_da,
    #        crosswalk_link_field: ds[crosswalk_link_field].values[gage_mask],
    #    }
    #    gage_link_df = pd.DataFrame(data = data_var_dict).set_index([crosswalk_gage_field])
      
    with xr.open_dataset(lastobsfile) as ds:        
        gages    = np.char.strip(ds[station_id].values)
        
        ref_time = datetime.datetime.strptime(ds.attrs[ref_t_attr_id], "%Y-%m-%d_%H:%M:%S")
        
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

    lastobs_df = (
        pd.DataFrame(data = data_var_dict).
        set_index('gages')
        #join(gage_link_df, how = 'inner').
        #reset_index().
        #set_index(crosswalk_link_field)
    )
    

    #lastobs_df = lastobs_df[
    #    [
    #        'gages',
    #        'time_since_lastobs',
    #        'lastobs_discharge',
    #    ]
    #]
    
    return lastobs_df

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