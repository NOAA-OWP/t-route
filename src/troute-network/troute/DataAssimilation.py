from abc import ABC, abstractmethod
import troute.nhd_io as nhd_io
import troute.nhd_network_utilities_v02 as nnu
import pandas as pd
import numpy as np
import pathlib
import time

#FIXME parameterize into construciton
showtiming = True
verbose = True

def build_data_assimilation_csv(data_assimilation_parameters):

    return nhd_io.get_usgs_df_from_csv(
        data_assimilation_parameters["data_assimilation_csv"],
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
    )

def build_data_assimilation_folder_v02(data_assimilation_parameters, run_parameters, network):

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
    
    return nhd_io.get_obs_from_timeslices(
        network.link_gage_df,
        # next two lines assume we've already verified that streamflow_da is true and contains the necessary dictionary items.
        #TODO: reorganize how DA objects are called/created in __main__.py to ensure we have all necessary info...
        data_assimilation_parameters.get('streamflow_da',None).get('crosswalk_gage_field','gages'),
        data_assimilation_parameters.get('streamflow_da',None).get('crosswalk_segID_field','link'),
        usgs_files,
        data_assimilation_parameters.get("qc_threshold", 1),
        data_assimilation_parameters.get("data_assimilation_interpolation_limit", 59),
        900, # 15 minutes, as secs
        run_parameters["t0"],
        network.cpu_pool, #verify
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
            loc[link_gage_df.index]
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
    
    """
    slots = ["_usgs_df", "_last_obs_df", "_reservoir_usgs_df", "_reservoir_usgs_param_df", "_reservoir_usace_df", "_reservoir_usace_param_df", "_da_parameter_dict", "_waterbody_types_df", "_usgs_lake_gage_crosswalk", "_usace_lake_gage_crosswalk"]
    def __init__(self, data_assimilation_parameters, run_parameters, waterbody_parameters, network, da_run):
        lastobs_df, da_parameter_dict = nnu.build_data_assimilation_lastobs(
            data_assimilation_parameters
        )
    
        # replace link ids with lake ids, for gages at waterbody outlets, 
        # otherwise, gage data will not be assimilated at waterbody outlet
        # segments.
        if network.link_lake_crosswalk:
            lastobs_df = _reindex_link_to_lake_id(lastobs_df, network.link_lake_crosswalk)
            
        self._last_obs_df = lastobs_df
        self._da_parameter_dict = da_parameter_dict
        # TODO: Add parameters here for interpolation length (14/59), QC threshold (1.0)
        
        #---------------------------------------------------------------------------
        # Determine which DA features to turn on
        #---------------------------------------------------------------------------

        # isolate user-input parameters for streamflow data assimilation
        streamflow_da_parameters = data_assimilation_parameters.get('streamflow_da', None)

        # determine if user explictly requests streamflow DA
        nudging = False
        if streamflow_da_parameters:
            nudging = streamflow_da_parameters.get('streamflow_nudging', False)

        usgs_timeslices_folder = data_assimilation_parameters.get('usgs_timeslices_folder',None)
        
        # isolate user-input parameters for reservoir data assimilation
        reservoir_da_parameters = data_assimilation_parameters.get('reservoir_da', None)

        # check if user explictly requests USGS and/or USACE reservoir DA
        usgs_persistence  = False
        usace_persistence = False
        if reservoir_da_parameters:
            usgs_persistence  = reservoir_da_parameters.get('reservoir_persistence_usgs', False)
            usace_persistence = reservoir_da_parameters.get('reservoir_persistence_usace', False)
            param_file   = reservoir_da.get('gage_lakeID_crosswalk_file',None)
        
        # check if RFC-type reservoirs are set to true
        rfc_params = waterbody_parameters.get('rfc')
        if rfc_params:
            rfc_forecast = rfc_params.get('reservoir_rfc_forecasts',False)
            param_file = rfc_params.get('reservoir_parameter_file',None)
        else:
            rfc_forecast = False
        
        level_pool_params = waterbody_parameters.get('level_pool', defaultdict(list))
                
                
        #--------------------------------------------------------------------------------
        # Assemble USGS observation data for Streamflow DA or USGS Reservoir Persistence
        #--------------------------------------------------------------------------------

        # if user requested nudging or usgs_persistence and a specified a USGS TimeSlice directory, 
        # then build and return USGS dataframe
        if (nudging or usgs_persistence) and usgs_timeslices_folder:

            start_time = time.time()
            LOG.info(
                "Creating a DataFrame of USGS gage observations ..."
            )
            
            self._usgs_df = _create_usgs_df(data_assimilation_parameters, streamflow_da_parameters, run_parameters, network, da_run)
            
            LOG.debug(
                "USGS observation DataFrame creation complete in %s seconds." \
                % (time.time() - start_time)
            )
            
        ##### is the following snippet of code needed?
        if not self._last_obs_df.index.empty:
            if not self._usgs_df.empty and not self._usgs_df.index.equals(self._last_obs_df.index):
                print("USGS Dataframe Index Does Not Match Last Observations Dataframe Index")
                self._usgs_df = self._usgs_df.loc[self._last_obs_df.index]
        #####
                
        #--------------------------------------------------------------------------------
        # Assemble Reservoir dataframes
        #--------------------------------------------------------------------------------
        
        # if any reservoir DA is turned on, load info from reservoir parameter file:
        break_network_at_waterbodies = waterbody_parameters.get("break_network_at_waterbodies", False)
        if not network.wbody_conn: 
            # Turn off any further reservoir processing if the network contains no waterbodies
            break_network_at_waterbodies = False
        
        if break_network_at_waterbodies:
            if (param_file and (usgs_persistence or usace_persistence)) or (param_file and rfc_forecast):
                waterbody_type_specified = True
                (
                    waterbody_types_df, 
                    usgs_lake_gage_crosswalk, 
                    usace_lake_gage_crosswalk
                ) = nhd_io.read_reservoir_parameter_file(
                    param_file,
                    usgs_persistence,
                    usace_persistence,
                    rfc_forecast,
                    level_pool_params.get("level_pool_waterbody_id", 'lake_id'),
                    reservoir_da_parameters.get('crosswalk_usgs_gage_field', 'usgs_gage_id'),
                    reservoir_da_parameters.get('crosswalk_usgs_lakeID_field', 'usgs_lake_id'),
                    reservoir_da_parameters.get('crosswalk_usace_gage_field', 'usace_gage_id'),
                    reservoir_da_parameters.get('crosswalk_usace_lakeID_field', 'usace_lake_id'),
                    network.wbody_conn.values(),
                )
            else:
                waterbody_type_specified = True
                waterbody_types_df = pd.DataFrame(data = 1, index = waterbodies_df.index, columns = ['reservoir_type'])
                usgs_lake_gage_crosswalk = None
                usace_lake_gage_crosswalk = None
        
        else:
            # Declare empty dataframes
            waterbody_types_df = pd.DataFrame()
            usgs_lake_gage_crosswalk = None
            usace_lake_gage_crosswalk = None
        
        self._waterbodies_df = waterbodies_df
        self._usgs_lake_gage_crosswalk = usgs_lake_gage_crosswalk
        self._usace_lake_gage_crosswalk = usace_lake_gage_crosswalk
        

        if usgs_persistence:
            if usgs_df.emtpy == False: # if usgs_df is already created, make reservoir_usgs_df from that rather than reading in data again
                
                gage_lake_df = (
                    usgs_lake_gage_crosswalk.
                    reset_index().
                    set_index(['usgs_gage_id']) # <- TODO use input parameter for this
                )
                
                # build dataframe that crosswalks gageIDs to segmentIDs
                gage_link_df = (
                    link_gage_df['gages'].
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
                    usgs_df.
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
                    lake_gage_crosswalk = usgs_lake_gage_crosswalk,
                    res_source = 'usgs')
        else:
            reservoir_usgs_df = pd.DataFrame()
            reservoir_usgs_param_df = pd.DataFrame()
            
        if usace_persistence:
            (
                reservoir_usace_df,
                reservoir_ucace_param_df
            ) = _create_reservoir_df(
                data_assimilation_parameters,
                reservoir_da_parameters,
                streamflow_da_parameters,
                run_parameters,
                network,
                da_run,
                lake_gage_crosswalk = usace_lake_gage_crosswalk,
                res_source = 'usace')
        else:
            reservoir_usace_df = pd.DataFrame()
            reservoir_usace_param_df = pd.DataFrame()
        
        self._reservoir_usgs_df = reservoir_usgs_df
        self._reservoir_usgs_param_df = reservoir_usgs_param_df
        self._reservoir_usace_df = reservoir_usace_df
        self._reservoir_usace_param_df = reservoir_usace_param_df


    @property
    def assimilation_parameters(self):
        return self._da_parameter_dict
    
    @property
    def last_obs(self):
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
    
    @property
    def waterbody_types_df(self):
        return self._waterbody_types_df
    
    @property
    def usgs_lake_gage_crosswalk(self):
        return self._usgs_lake_gage_crosswalk
    
    @property
    def usace_lake_gage_crosswalk(self):
        return self._usace_lake_gage_crosswalk