from abc import ABC, abstractmethod
import troute.nhd_io as nhd_io
import pandas as pd
import numpy as np
import pathlib
import time
import xarray as xr
from collections import defaultdict

#FIXME parameterize into construciton
showtiming = True
verbose = True

def build_data_assimilation_csv(data_assimilation_parameters):

    return nhd_io.get_usgs_df_from_csv(
        data_assimilation_parameters["data_assimilation_csv"],
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
    )

def build_data_assimilation_folder(data_assimilation_parameters, run_parameters):

    usgs_timeslices_folder = pathlib.Path(
        data_assimilation_parameters["data_assimilation_timeslices_folder"],
    ).resolve()
    if "data_assimilation_filter" in data_assimilation_parameters:
        da_glob_filter = data_assimilation_parameters["data_assimilation_filter"]
        usgs_files = sorted(usgs_timeslices_folder.glob(da_glob_filter))
    elif "usgs_timeslice_files" in data_assimilation_parameters:
        usgs_files = data_assimilation_parameters.get("usgs_timeslice_files", None)
        usgs_files = [usgs_timeslices_folder.joinpath(f) for f in usgs_files]
    else:
        print("No Files Found for DA")
        # TODO: Handle this with a real exception

    return nhd_io.get_usgs_from_time_slices_folder(
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
        usgs_files,
        data_assimilation_parameters.get("qc_threshold", 1),
        data_assimilation_parameters.get("data_assimilation_interpolation_limit", 59),
        #FIXME/TODO collapse these into da parameters
        run_parameters["dt"],
        run_parameters["t0"],
    )

def _reindex_link_to_lake_id(target_df, crosswalk):
    '''
    Utility function for replacing link ID index values
    with lake ID values in a dataframe. This is used to 
    reinedex dataframes used for streamflow DA such that 
    data from data from gages located at waterbody outlets
    can be assimilated. 
    
    Arguments:
    ----------
    - target_df (DataFrame): Data frame to be reinexed
    - crosswalk      (dict): Relates lake ids to outlet link ids
    
    Returns:
    --------
    - target_df (DataFrame): Re-indexed with lake ids replacing 
                             link ids
    '''

    # evaluate intersection of link ids and target_df index values
    # i.e. what are the index positions of link ids that need replacing?
    linkids = np.fromiter(crosswalk.values(), dtype = int)
    gageidxs = target_df.index.to_numpy()
    lake_index_intersect = np.intersect1d(
        gageidxs, 
        linkids, 
        return_indices = True
    )

    # replace link ids with lake IDs in the target_df index array
    lakeids = np.fromiter(crosswalk.keys(), dtype = int)
    gageidxs[lake_index_intersect[1]] = lakeids[lake_index_intersect[2]]

    # (re) set the target_df index
    target_df.set_index(gageidxs, inplace = True)
    
    return target_df

def _create_usgs_df(data_assimilation_parameters, streamflow_da_parameters, run_parameters, network, da_run):
    '''
    Function for reading USGS timeslice files and creating a dataframe
    of USGS gage observations. This dataframe is used for streamflow
    nudging and can be used for constructing USGS reservoir dataframes.
    
    Arguments:
    ----------
    - data_assimilation_parameters (dict): user input data re data assimilation
    - streamflow_da_parameters     (dict): user input data re streamflow nudging
    - run_parameters               (dict): user input data re subset of compute configuration
    - network                    (Object): network object created from abstract class
    - da_run                       (list): list of data assimilation files separated by for loop chunks
    
    Returns:
    --------
    - usgs_df (DataFrame): dataframe of USGS gage observations
    '''
    usgs_timeslices_folder = data_assimilation_parameters.get("usgs_timeslices_folder", None)
    lastobs_file           = streamflow_da_parameters.get("wrf_hydro_lastobs_file", None)
    lastobs_start          = data_assimilation_parameters.get("wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time",0)
    lastobs_type           = data_assimilation_parameters.get("wrf_lastobs_type", "error-based")
    crosswalk_file         = streamflow_da_parameters.get("gage_segID_crosswalk_file", None)
    crosswalk_gage_field   = streamflow_da_parameters.get('crosswalk_gage_field','gages')
    crosswalk_segID_field  = streamflow_da_parameters.get('crosswalk_segID_field','link')
    da_decay_coefficient   = data_assimilation_parameters.get("da_decay_coefficient",120)
    qc_threshold           = data_assimilation_parameters.get("qc_threshold",1)
    interpolation_limit    = data_assimilation_parameters.get("interpolation_limit_min",59)
    
    # TODO: join timeslice folder and files into complete path upstream
    usgs_timeslices_folder = pathlib.Path(usgs_timeslices_folder)
    usgs_files = [usgs_timeslices_folder.joinpath(f) for f in 
                  da_run['usgs_timeslice_files']]
	
    if usgs_files:
        usgs_df = (
            nhd_io.get_obs_from_timeslices(
                network.link_gage_df,
                crosswalk_gage_field,
                crosswalk_segID_field,
                usgs_files,
                qc_threshold,
                interpolation_limit,
                run_parameters.get("dt"),
                network.t0,
                run_parameters.get("cpu_pool", None)
            ).
            loc[network.link_gage_df.index]
        )
		
		# replace link ids with lake ids, for gages at waterbody outlets, 
		# otherwise, gage data will not be assimilated at waterbody outlet
		# segments.
        if network.link_lake_crosswalk:
            usgs_df = _reindex_link_to_lake_id(usgs_df, network.link_lake_crosswalk)
    
    else:
        usgs_df = pd.DataFrame()
    
    return usgs_df

def _create_reservoir_df(data_assimilation_parameters, reservoir_da_parameters, streamflow_da_parameters, run_parameters, network, da_run, lake_gage_crosswalk, res_source):
    '''
    Function for reading USGS/USACE timeslice files and creating a dataframe
    of reservoir observations and initial parameters. 
    These dataframes are used for reservoir DA.
    
    Arguments:
    ----------
    - data_assimilation_parameters (dict): user input data re data assimilation
    - reservoir_da_parameters      (dict): user input data re reservoir data assimilation
    - streamflow_da_parameters     (dict): user input data re streamflow nudging
    - run_parameters               (dict): user input data re subset of compute configuration
    - network                    (Object): network object created from abstract class
    - da_run                       (list): list of data assimilation files separated
                                           by for loop chunks
    - lake_gage_crosswalk          (dict): usgs/usace gage ids and corresponding segment ids at
                                           which they are located
    - res_source                    (str): either 'usgs' or 'usace', specifiying which type of
                                           reservoir dataframe to create (must match lake_gage_crosswalk
    
    Returns:
    --------
    - reservoir_usgs/usace_df       (DataFrame): USGS/USACE reservoir observations
    - reservoir_usgs/usace_param_df (DataFrame): USGS/USACE reservoir hybrid DA initial parameters
    '''
    res_timeslices_folder  = data_assimilation_parameters.get(res_source + "_timeslices_folder",None)
    crosswalk_file         = reservoir_da_parameters.get("gage_lakeID_crosswalk_file", None)
    crosswalk_gage_field   = streamflow_da_parameters.get('crosswalk_' + res_source + '_gage_field',res_source + '_gage_id')
    crosswalk_lakeID_field = streamflow_da_parameters.get('crosswalk_' + res_source + '_lakeID_field',res_source + '_lake_id')
    qc_threshold           = data_assimilation_parameters.get("qc_threshold",1)
    interpolation_limit    = data_assimilation_parameters.get("interpolation_limit_min",59)
	
    # TODO: join timeslice folder and files into complete path upstream in workflow
    res_timeslices_folder = pathlib.Path(res_timeslices_folder)
    res_files = [res_timeslices_folder.joinpath(f) for f in
                 da_run[res_source + '_timeslice_files']]
			
    if res_files:
		
        reservoir_df = nhd_io.get_obs_from_timeslices(
            lake_gage_crosswalk,
            crosswalk_gage_field,
            crosswalk_lakeID_field,
            res_files,
            qc_threshold,
            interpolation_limit,
            900,                      # 15 minutes, as secs
            network.t0,
            run_parameters.get("cpu_pool", None)
        )
		
    else:
        reservoir_df = pd.DataFrame() 
	
    # create reservoir hybrid DA initial parameters dataframe    
    if reservoir_df.empty == False:
        reservoir_param_df = pd.DataFrame(
            data = 0, 
            index = reservoir_df.index ,
            columns = ['update_time']
        )
        reservoir_param_df['prev_persisted_outflow'] = np.nan
        reservoir_param_df['persistence_update_time'] = 0
        reservoir_param_df['persistence_index'] = 0
    else:
        reservoir_param_df = pd.DataFrame()
        
    return reservoir_df, reservoir_param_df
    
def _set_reservoir_da_params(run_results):
    '''
    Update persistence reservoir DA parameters for subsequent loops
    Arguments:
    ----------
    - run_results (list): output from the compute kernel sequence, organized
        (because that is how it comes out of the kernel) by network.
        For each item in the result, there are seven elements, the
        fifth (usgs) and sixth (usace) of which are lists of five elements 
        containing: 1) a list of the segments ids where data assimilation 
        was performed (if any) in that network; 2) a list of the lupdate time; 
        3) a list of the previously persisted outflow; 4) a list of the 
        persistence index; 5) a list of hte persistence update time.
    
    Returns:
    --------
    - reservoir_usgs_param_df (DataFrame): USGS reservoir DA parameters
    - reservoir_usace_param_df (DataFrame): USACE reservoir DA parameters
    '''
    
    reservoir_usgs_param_df = pd.DataFrame(data = [], 
                                           index = [], 
                                           columns = [
                                               'update_time', 'prev_persisted_outflow', 
                                               'persistence_update_time', 'persistence_index'
                                           ]
                                          )
    reservoir_usace_param_df = pd.DataFrame(data = [], 
                                           index = [], 
                                           columns = [
                                               'update_time', 'prev_persisted_outflow', 
                                               'persistence_update_time', 'persistence_index'
                                           ]
                                          )
    
    for r in run_results:
        
        if len(r[4][0]) > 0:
            tmp_usgs = pd.DataFrame(data = r[4][1], index = r[4][0], columns = ['update_time'])
            tmp_usgs['prev_persisted_outflow'] = r[4][2]
            tmp_usgs['persistence_update_time'] = r[4][4]
            tmp_usgs['persistence_index'] = r[4][3]
            reservoir_usgs_param_df = pd.concat([reservoir_usgs_param_df, tmp_usgs])
        
        if len(r[5][0]) > 0:
            tmp_usace = pd.DataFrame(data = r[5][1], index = r[5][0], columns = ['update_time'])
            tmp_usace['prev_persisted_outflow'] = r[5][2]
            tmp_usace['persistence_update_time'] = r[5][4]
            tmp_usace['persistence_index'] = r[5][3]
            reservoir_usace_param_df = pd.concat([reservoir_usace_param_df, tmp_usace])
    
    return reservoir_usgs_param_df, reservoir_usace_param_df

def build_lastobs_df(
        lastobsfile,
        crosswalk_file,
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
    with xr.open_dataset(crosswalk_file) as ds:
        gage_list = list(map(bytes.strip, ds[crosswalk_gage_field].values))
        gage_mask = list(map(bytes.isalnum, gage_list))
        gage_da   = list(map(bytes.strip, ds[crosswalk_gage_field][gage_mask].values))
        data_var_dict = {
            crosswalk_gage_field: gage_da,
            crosswalk_link_field: ds[crosswalk_link_field].values[gage_mask],
        }
        gage_link_df = pd.DataFrame(data = data_var_dict).set_index([crosswalk_gage_field])
            
    with xr.open_dataset(lastobsfile) as ds:
        
        gages    = np.char.strip(ds[station_id].values)
        
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

    lastobs_df = (
        pd.DataFrame(data = data_var_dict).
        set_index('gages').
        join(gage_link_df, how = 'inner').
        reset_index().
        set_index(crosswalk_link_field)
    )
    lastobs_df = lastobs_df[
        [
            'gages',
            'time_since_lastobs',
            'lastobs_discharge',
        ]
    ]
    
    return lastobs_df

def new_lastobs(run_results, time_increment):
    """
    Creates new "lastobs" dataframe for the next simulation chunk.

    Arguments:
    ----------
    - run_results (list): output from the compute kernel sequence, organized
        (because that is how it comes out of the kernel) by network.
        For each item in the result, there are seven elements, the
        fourth of which is a tuple containing: 1) a list of the
        segments ids where data assimilation was performed (if any)
        in that network; 2) a list of the last valid observation
        applied at that segment; 3) a list of the time in seconds
        from the beginning of the last simulation that the
        observation was applied.
    - time_increment (int): length of the prior simulation. To prepare the
        next lastobs state, we have to convert the time since the prior
        simulation start to a time since the new simulation start.
        If the most recent observation was right at the end of the
        prior loop, then the value in the incoming run_result will
        be equal to the time_increment and the output value will be
        zero. If observations were not present at the last timestep,
        the last obs time will be calculated to a negative value --
        the number of seconds ago that the last valid observation
        was used for assimilation.
        
    Returns:
    --------
    - lastobs_df (DataFrame): Last gage observations data for DA
    """
    df = pd.concat(
        [
            pd.DataFrame(
                # TODO: Add time_increment (or subtract?) from time_since_lastobs
                np.array([rr[3][1],rr[3][2]]).T,
                index=rr[3][0],
                columns=["time_since_lastobs", "lastobs_discharge"]
            )
            for rr in run_results
            if not rr[3][0].size == 0
        ],
        copy=False,
    )
    df["time_since_lastobs"] = df["time_since_lastobs"] - time_increment

    return df

def read_reservoir_parameter_file(
    reservoir_parameter_file, 
    usgs_hybrid,
    usace_hybrid,
    rfc_forecast,
    lake_index_field="lake_id", 
    usgs_gage_id_field = "usgs_gage_id",
    usgs_lake_id_field = "usgs_lake_id",
    usace_gage_id_field = "usace_gage_id",
    usace_lake_id_field = "usace_lake_id",
    lake_id_mask=None,
):

    """
    Reads reservoir parameter file, which is separate from the LAKEPARM file.
    Extracts reservoir "type" codes and returns in a DataFrame
    type 1: Levelool
    type 2: USGS Hybrid Persistence
    type 3: USACE Hybrid Persistence
    type 4: RFC
    This function is only called if Hybrid Persistence or RFC type reservoirs
    are active.
    
    Arguments
    ---------
    - reservoir_parameter_file (str): full file path of the reservoir parameter
                                      file
    
    - usgs_hybrid          (boolean): If True, then USGS Hybrid DA will be coded
    
    - usace_hybrid         (boolean): If True, then USACE Hybrid DA will be coded
    
    - rfc_forecast         (boolean): If True, then RFC Forecast DA will be coded
    
    - lake_index_field         (str): field containing lake IDs in reservoir 
                                      parameter file
    
    - lake_id_mask     (dict_values): Waterbody IDs in the model domain 
    
    Returns
    -------
    - df1 (Pandas DataFrame): Reservoir type codes, indexed by lake_id
    
    Notes
    -----
    
    """
    with xr.open_dataset(reservoir_parameter_file) as ds:
        ds = ds.swap_dims({"feature_id": lake_index_field})
        ds_new = ds["reservoir_type"]
        df1 = ds_new.sel({lake_index_field: list(lake_id_mask)}).to_dataframe()
        
        ds_vars = [i for i in ds.data_vars] 
        
        if (usgs_gage_id_field in ds_vars) and (usgs_lake_id_field in ds_vars):
            usgs_crosswalk = pd.DataFrame(
                data = ds[usgs_gage_id_field].to_numpy(), 
                index = ds[usgs_lake_id_field].to_numpy(), 
                columns = [usgs_gage_id_field]
            )
            usgs_crosswalk.index.name = usgs_lake_id_field
        else:
            usgs_crosswalk = None
        
        if (usace_gage_id_field in ds_vars) and (usace_lake_id_field in ds_vars):
            usace_crosswalk = pd.DataFrame(
                data = ds[usace_gage_id_field].to_numpy(), 
                index = ds[usace_lake_id_field].to_numpy(), 
                columns = [usace_gage_id_field]
            )
            usace_crosswalk.index.name = usace_lake_id_field
        else:
            usace_crosswalk = None
        
    # drop duplicate indices
    df1 = (df1.reset_index()
           .drop_duplicates(subset="lake_id")
           .set_index("lake_id")
          )
    
    # recode to levelpool (1) for reservoir DA types set to false
    if usgs_hybrid == False:
        df1[df1['reservoir_type'] == 2] = 1
    if usace_hybrid == False:
        df1[df1['reservoir_type'] == 3] = 1
    if rfc_forecast == False:
        df1[df1['reservoir_type'] == 4] = 1
    
    return df1, usgs_crosswalk, usace_crosswalk

class DataAssimilation(ABC):
    """
    
    """

    @property
    @abstractmethod
    def usgs_df(self):
        pass

class NudgingDA(DataAssimilation):
    """
    
    """
    slots = ["_usgs_df", "_lastobs_df", "_da_params"]
    def __init__(self, data_assimilation_parameters, run_parameters):
        data_assimilation_csv = data_assimilation_parameters.get(
            "data_assimilation_csv", None
        )
        data_assimilation_folder = data_assimilation_parameters.get(
            "data_assimilation_timeslices_folder", None
        )
        # TODO: Copy comments from nhd_network_utilitys `build_data_assimilation_lastobs`
        lastobs_file = data_assimilation_parameters.get("wrf_hydro_lastobs_file", None)
        lastobs_start = data_assimilation_parameters.get(
            "wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0
        )
        lastobs_type = data_assimilation_parameters.get("wrf_lastobs_type", "error-based")
        lastobs_crosswalk_file = data_assimilation_parameters.get(
            "wrf_hydro_da_channel_ID_crosswalk_file", None
        )

        self._da_params = {}

        if data_assimilation_csv or data_assimilation_folder or lastobs_file:
            self._da_params["da_decay_coefficient"] = data_assimilation_parameters.get("da_decay_coefficient", 120)
            # TODO: Add parameters here for interpolation length (14/59), QC threshold (1.0)

            if showtiming:
                start_time = time.time()
            if verbose:
                print("creating usgs time_slice data array ...")
            self._last_obs_df = nhd_io.build_lastobs_df(
                lastobs_file,
                lastobs_crosswalk_file,
                lastobs_type,  # TODO: Confirm that we are using this; delete it if not.
                lastobs_start,
            )
            if data_assimilation_csv:
                self._usgs_df = build_data_assimilation_csv(data_assimilation_parameters)
            elif data_assimilation_folder:
                self._usgs_df = build_data_assimilation_folder(data_assimilation_parameters, run_parameters)

            if not self._last_obs_df.index.empty:
                if not self._usgs_df.empty and not self._usgs_df.index.equals(self._last_obs_df.index):
                    print("USGS Dataframe Index Does Not Match Last Observations Dataframe Index")
                    self._usgs_df = self._usgs_df.loc[self._last_obs_df.index]
            if verbose:
                print("usgs array complete")
            if showtiming:
                print("... in %s seconds." % (time.time() - start_time))
        else:
            self._last_obs_df = pd.DataFrame()
            self._usgs_df = pd.DataFrame()

    @property
    def asssimilation_parameters(self):
        return self._da_params
    
    @property
    def last_obs(self):
        return self._last_obs_df

    @property
    def usgs_df(self):
        return self._usgs_df

    
class AllDA(DataAssimilation):
    """
    Object containing all data assimilation information. Created from user input
    and a network object that defines the network and forcing parameters. TODO: this
    could be broken up into separate objects specific to each type of DA (streamflow
    nudging, reservoir persistence, etc.). But for now it is all done here.
    
    Arguments:
    ----------
    - data_assimilation_parameters (dict): user input data re data assimilation
    - run_parameters               (dict): user input data re subset of compute configuration
    - waterbody_parameters         (dict): user input data re waterbody parameters
    - network                    (Object): network object created from abstract class
    - da_run                       (list): list of data assimilation files separated
                                           by for loop chunks
    
    Returns:
    --------
    - data_assimilation               (Object): Object containing all data assimilation information
        - usgs_df                  (DataFrame): dataframe of USGS gage observations
        - lastobs_df               (DataFrame): Last gage observations data for DA
        - reservoir_usgs_df        (DataFrame): USGS reservoir observations
        - reservoir_usgs_param_df  (DataFrame): USGS reservoir DA parameters
        - reservoir_usace_df       (DataFrame): USACE reservoir observations
        - reservoir_usace_param_df (DataFrame): USACE reservoir DA parameters
        - da_parameter_dict             (dict): user input data re data assimilation
    """
    __slots__ = ["_usgs_df", "_last_obs_df", "_reservoir_usgs_df", "_reservoir_usgs_param_df", "_reservoir_usace_df", "_reservoir_usace_param_df", "_da_parameter_dict"]
    def __init__(self, data_assimilation_parameters, run_parameters, waterbody_parameters, network, da_run):
                
        #---------------------------------------------------------------------------
        # Determine which DA features to turn on
        #---------------------------------------------------------------------------

        # isolate user-input parameters for streamflow data assimilation
        streamflow_da_parameters = data_assimilation_parameters.get('streamflow_da', None)

        # determine if user explictly requests streamflow DA
        nudging = False
        if streamflow_da_parameters:
            nudging = streamflow_da_parameters.get('streamflow_nudging', False)

        usgs_timeslices_folder = data_assimilation_parameters.get('usgs_timeslices_folder', None)
        
        # isolate user-input parameters for reservoir data assimilation
        reservoir_da_parameters = data_assimilation_parameters.get('reservoir_da', None)

        # check if user explictly requests USGS and/or USACE reservoir DA
        usgs_persistence  = False
        usace_persistence = False
        if reservoir_da_parameters:
            usgs_persistence  = reservoir_da_parameters.get('reservoir_persistence_usgs', False)
            usace_persistence = reservoir_da_parameters.get('reservoir_persistence_usace', False)
            param_file   = reservoir_da_parameters.get('gage_lakeID_crosswalk_file', None)
        
        # check if RFC-type reservoirs are set to true
        rfc_params = waterbody_parameters.get('rfc')
        if rfc_params:
            rfc_forecast = rfc_params.get('reservoir_rfc_forecasts',False)
            param_file = rfc_params.get('reservoir_parameter_file',None)
        else:
            rfc_forecast = False
        
        level_pool_params = waterbody_parameters.get('level_pool', defaultdict(list))
                
        
        #--------------------------------------------------------------------------------
        # Produce lastobs_df if streamflow nudging is turned on
        #--------------------------------------------------------------------------------
        lastobs_df = pd.DataFrame()
        da_parameter_dict = {}
        
        if nudging:
            lastobs_file = streamflow_da_parameters.get("wrf_hydro_lastobs_file", None)
            lastobs_crosswalk_file = streamflow_da_parameters.get("gage_segID_crosswalk_file", None)
            lastobs_start = streamflow_da_parameters.get("wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0)
            
            if lastobs_file:
                lastobs_df = build_lastobs_df(
                    lastobs_file,
                    lastobs_crosswalk_file,
                    lastobs_start,
                )
            
            da_parameter_dict["da_decay_coefficient"] = data_assimilation_parameters.get("da_decay_coefficient", 120)
            da_parameter_dict["diffusive_streamflow_nudging"] = streamflow_da_parameters.get("diffusive_streamflow_nudging", False)

            # replace link ids with lake ids, for gages at waterbody outlets, 
            # otherwise, gage data will not be assimilated at waterbody outlet
            # segments.
            if network.link_lake_crosswalk:
                lastobs_df = _reindex_link_to_lake_id(lastobs_df, network.link_lake_crosswalk)
            
        self._last_obs_df = lastobs_df
        self._da_parameter_dict = da_parameter_dict
        # TODO: Add parameters here for interpolation length (14/59), QC threshold (1.0)
        
        
        #--------------------------------------------------------------------------------
        # Assemble USGS observation data for Streamflow DA or USGS Reservoir Persistence
        #--------------------------------------------------------------------------------

        # if user requested nudging or usgs_persistence and specified a USGS TimeSlice directory, 
        # then build and return USGS dataframe
        if (nudging or usgs_persistence) and usgs_timeslices_folder:

            self._usgs_df = _create_usgs_df(data_assimilation_parameters, streamflow_da_parameters, run_parameters, network, da_run)
                
        #--------------------------------------------------------------------------------
        # Assemble Reservoir dataframes
        #--------------------------------------------------------------------------------
        
        if usgs_persistence:
            # if usgs_df is already created, make reservoir_usgs_df from that rather than reading in data again
            if self._usgs_df.empty == False: 
                
                gage_lake_df = (
                    network.usgs_lake_gage_crosswalk.
                    reset_index().
                    set_index(['usgs_gage_id']) # <- TODO use input parameter for this
                )
                
                # build dataframe that crosswalks gageIDs to segmentIDs
                gage_link_df = (
                    network.link_gage_df['gages'].
                    reset_index().
                    set_index(['gages'])
                )
                
                # build dataframe that crosswalks segmentIDs to lakeIDs
                link_lake_df = (
                    gage_lake_df.
                    join(gage_link_df, how = 'inner').
                    reset_index().set_index('link').
                    drop(['index'], axis = 1)
                )

                # resample `usgs_df` to 15 minute intervals
                usgs_df_15min = (
                    self._usgs_df.
                    transpose().
                    resample('15min').asfreq().
                    transpose()
                )

                # subset and re-index `usgs_df`, using the segID <> lakeID crosswalk
                reservoir_usgs_df = (
                    usgs_df_15min.join(link_lake_df, how = 'inner').
                    reset_index().
                    set_index('usgs_lake_id').
                    drop(['index'], axis = 1)
                )
                
                # create reservoir hybrid DA initial parameters dataframe    
                if reservoir_usgs_df.empty == False:
                    reservoir_usgs_param_df = pd.DataFrame(
                        data = 0, 
                        index = reservoir_usgs_df.index ,
                        columns = ['update_time']
                    )
                    reservoir_usgs_param_df['prev_persisted_outflow'] = np.nan
                    reservoir_usgs_param_df['persistence_update_time'] = 0
                    reservoir_usgs_param_df['persistence_index'] = 0
                else:
                    reservoir_usgs_param_df = pd.DataFrame()
                
            else:
                (
                    reservoir_usgs_df,
                    reservoir_usgs_param_df
                ) = _create_reservoir_df(
                    data_assimilation_parameters,
                    reservoir_da_parameters,
                    streamflow_da_parameters,
                    run_parameters,
                    network,
                    da_run,
                    lake_gage_crosswalk = network.usgs_lake_gage_crosswalk,
                    res_source = 'usgs')
        else:
            reservoir_usgs_df = pd.DataFrame()
            reservoir_usgs_param_df = pd.DataFrame()
            
        if usace_persistence:
            (
                reservoir_usace_df,
                reservoir_usace_param_df
            ) = _create_reservoir_df(
                data_assimilation_parameters,
                reservoir_da_parameters,
                streamflow_da_parameters,
                run_parameters,
                network,
                da_run,
                lake_gage_crosswalk = network.usace_lake_gage_crosswalk,
                res_source = 'usace')
        else:
            reservoir_usace_df = pd.DataFrame()
            reservoir_usace_param_df = pd.DataFrame()
        
        self._reservoir_usgs_df = reservoir_usgs_df
        self._reservoir_usgs_param_df = reservoir_usgs_param_df
        self._reservoir_usace_df = reservoir_usace_df
        self._reservoir_usace_param_df = reservoir_usace_param_df
        
        # Trim the time-extent of the streamflow_da usgs_df
        # what happens if there are timeslice files missing on the front-end? 
        # if the first column is some timestamp greater than t0, then this will throw
        # an error. Need to think through this more. 
        if not self._usgs_df.empty:
            self._usgs_df = self._usgs_df.loc[:,network.t0:]


    def update(self, run_results, data_assimilation_parameters, run_parameters, network, da_run):
        '''
        Function to update data assimilation object for the next loop iteration.
        
        Arguments:
        ----------
        - run_results                  (list): output from the compute kernel sequence,
                                               organized (because that is how it comes 
                                               out of the kernel) by network.
                                               For each item in the result, there are 
                                               seven elements, the fifth (usgs) and sixth 
                                               (usace) of which are lists of five elements 
                                               containing: 1) a list of the segments ids 
                                               where data assimilation was performed (if any) 
                                               in that network; 2) a list of the lupdate time; 
                                               3) a list of the previously persisted outflow; 
                                               4) a list of the persistence index; 5) a list 
                                               of the persistence update time.
        - data_assimilation_parameters (dict): user input data re data assimilation
        - run_parameters               (dict): user input data re subset of compute configuration
        - network                    (Object): network object created from abstract class
        - da_run                       (list): list of data assimilation files separated
                                               by for loop chunks
        
        Returns:
        --------
        - data_assimilation               (Object): Object containing all data assimilation information
            - usgs_df                  (DataFrame): dataframe of USGS gage observations
            - lastobs_df               (DataFrame): Last gage observations data for DA
            - reservoir_usgs_df        (DataFrame): USGS reservoir observations
            - reservoir_usgs_param_df  (DataFrame): USGS reservoir DA parameters
            - reservoir_usace_df       (DataFrame): USACE reservoir observations
            - reservoir_usace_param_df (DataFrame): USACE reservoir DA parameters
        '''
        # get reservoir DA initial parameters for next loop itteration
        self._reservoir_usgs_param_df, self._reservoir_usace_param_df = _set_reservoir_da_params(run_results)
        
        # update usgs_df if it is not empty
        streamflow_da_parameters = data_assimilation_parameters.get('streamflow_da', None)
        reservoir_da_parameters = data_assimilation_parameters.get('reservoir_da', None)
        
        if streamflow_da_parameters.get('streamflow_nudging', False):
            self._usgs_df = _create_usgs_df(data_assimilation_parameters, streamflow_da_parameters, run_parameters, network, da_run)
            
            if reservoir_da_parameters.get('reservoir_persistence_usgs', False):
                
                gage_lake_df = (
                    network.usgs_lake_gage_crosswalk.
                    reset_index().
                    set_index(['usgs_gage_id']) # <- TODO use input parameter for this
                )
                
                # build dataframe that crosswalks gageIDs to segmentIDs
                gage_link_df = (
                    network.link_gage_df['gages'].
                    reset_index().
                    set_index(['gages'])
                )
                
                # build dataframe that crosswalks segmentIDs to lakeIDs
                link_lake_df = (
                    gage_lake_df.
                    join(gage_link_df, how = 'inner').
                    reset_index().set_index('link').
                    drop(['index'], axis = 1)
                )
                
                # resample `usgs_df` to 15 minute intervals
                usgs_df_15min = (
                    self._usgs_df.
                    transpose().
                    resample('15min').asfreq().
                    transpose()
                )
                
                # subset and re-index `usgs_df`, using the segID <> lakeID crosswalk
                self._reservoir_usgs_df = (
                    usgs_df_15min.join(link_lake_df, how = 'inner').
                    reset_index().
                    set_index('usgs_lake_id').
                    drop(['index'], axis = 1)
                )
        
        elif reservoir_da_parameters.get('reservoir_persistence_usgs', False):
            (
                self._reservoir_usgs_df,
                _,
            ) = _create_reservoir_df(
                data_assimilation_parameters,
                reservoir_da_parameters,
                streamflow_da_parameters,
                run_parameters,
                network,
                da_run,
                lake_gage_crosswalk = network.usgs_lake_gage_crosswalk,
                res_source = 'usgs')
        
        # USACE
        if reservoir_da_parameters.get('reservoir_persistence_usace', False):
            
            (
                self._reservoir_usace_df,
                _,
            ) = _create_reservoir_df(
                data_assimilation_parameters,
                reservoir_da_parameters,
                streamflow_da_parameters,
                run_parameters,
                network,
                da_run,
                lake_gage_crosswalk = network.usace_lake_gage_crosswalk,
                res_source = 'usace')
        
        # if there are no TimeSlice files available for hybrid reservoir DA in the next loop, 
        # but there are DA parameters from the previous loop, then create a
        # dummy observations df. This allows the reservoir persistence to continue across loops.
        # USGS Reservoirs
        if not network.waterbody_types_df.empty:
            if 2 in network.waterbody_types_df['reservoir_type'].unique():
                if self._reservoir_usgs_df.empty and len(self._reservoir_usgs_param_df.index) > 0:
                    self._reservoir_usgs_df = pd.DataFrame(
                        data    = np.nan, 
                        index   = self._reservoir_usgs_param_df.index, 
                        columns = [network.t0]
                    )

            # USACE Reservoirs   
            if 3 in network.waterbody_types_df['reservoir_type'].unique():
                if self._reservoir_usace_df.empty and len(self._reservoir_usace_param_df.index) > 0:
                    self._reservoir_usace_df = pd.DataFrame(
                        data    = np.nan, 
                        index   = self._reservoir_usace_param_df.index, 
                        columns = [network.t0]
                    )

            '''
            # update RFC lookback hours if there are RFC-type reservoirs in the simulation domain
            if 4 in network.waterbody_types_df['reservoir_type'].unique():
                waterbody_parameters = update_lookback_hours(run_parameters.get("dt"), run_parameters.get("nts"), waterbody_parameters)
            '''
        
        if streamflow_da_parameters:
            if streamflow_da_parameters.get('streamflow_nudging', False):
                self._last_obs_df = new_lastobs(run_results, run_parameters.get("dt") * run_parameters.get("nts"))
        
        # Trim the time-extent of the streamflow_da usgs_df
        # what happens if there are timeslice files missing on the front-end? 
        # if the first column is some timestamp greater than t0, then this will throw
        # an error. Need to think through this more. 
        if not self._usgs_df.empty:
            self._usgs_df = self._usgs_df.loc[:,network.t0:]
    
    @property
    def assimilation_parameters(self):
        return self._da_parameter_dict
    
    @property
    def lastobs_df(self):
        return self._last_obs_df

    @property
    def usgs_df(self):
        return self._usgs_df
    
    @property
    def reservoir_usgs_df(self):
        return self._reservoir_usgs_df
    
    @property
    def reservoir_usgs_param_df(self):
        return self._reservoir_usgs_param_df
    
    @property
    def reservoir_usace_df(self):
        return self._reservoir_usace_df
    
    @property
    def reservoir_usace_param_df(self):
        return self._reservoir_usace_param_df
