import time
import pandas as pd
import pathlib
import numpy as np
import pathlib
import xarray as xr
from datetime import datetime
from collections import defaultdict
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io
import logging

LOG = logging.getLogger('')


def nwm_network_preprocess(
    supernetwork_parameters,
    waterbody_parameters,
    preprocessing_parameters,
    compute_parameters,
    data_assimilation_parameters,
):

    hybrid_params = compute_parameters.get("hybrid_parameters", False)
    if hybrid_params:
        
        domain_file    = hybrid_params.get("diffusive_domain",   None)
        run_hybrid     = hybrid_params.get('run_hybrid_routing', False)
        topobathy_file = hybrid_params.get("topobathy_domain",   None)
        use_topobathy  = hybrid_params.get('use_natl_xsections', False)
        
        if domain_file and run_hybrid:
            
            # read diffusive domain dictionary from yaml or json
            diffusive_domain = nhd_io.read_diffusive_domain(domain_file)
            
            if topobathy_file and use_topobathy:
                
                # read topobathy domain netcdf file, set index to 'comid'
                # TODO: replace 'comid' with a user-specified indexing variable name.
                # ... if for whatever reason there is not a `comid` variable in the 
                # ... dataframe returned from read_netcdf, then the code would break here.
                topobathy_data = (nhd_io.read_netcdf(topobathy_file).set_index('comid'))
                
                # TODO: Request GID make comID variable an integer in their product, so
                # we do not need to change variable types, here.
                topobathy_data.index = topobathy_data.index.astype(int)
                
            else:
                topobathy_data = pd.DataFrame()
                LOG.debug('No natural cross section topobathy data provided.')
             
            # initialize a dictionary to hold network data for each of the diffusive domains
            diffusive_network_data = {}
        
        else:
            diffusive_domain = None
            diffusive_network_data = None
            topobathy_data = pd.DataFrame()
            LOG.debug('No diffusive domain file spefified in configuration file.')
    else:
        diffusive_domain = None
        diffusive_network_data = None
        topobathy_data = pd.DataFrame()

    LOG.info("creating supernetwork connections set")

    start_time = time.time()

    # STEP 1: Build basic network connections graph,
    # read network parameters, identify waterbodies and gages, if any.
    connections, param_df, wbody_conn, gages = nnu.build_connections(
        supernetwork_parameters,
    )
    
    break_network_at_waterbodies = waterbody_parameters.get(
        "break_network_at_waterbodies", False
    )
    break_network_at_gages = supernetwork_parameters.get(
        "break_network_at_gages", False
    )

    if (
        not wbody_conn
    ):  # Turn off any further reservoir processing if the network contains no waterbodies
        break_network_at_waterbodies = False

    if break_network_at_waterbodies:
        connections = nhd_network.replace_waterbodies_connections(
            connections, wbody_conn
        )

    
    LOG.debug("supernetwork connections set complete in %s seconds." % (time.time() - start_time))

    ################################
    ## STEP 3a: Read waterbody parameter file
    # waterbodies_values = supernetwork_values[12]
    # waterbodies_segments = supernetwork_values[13]
    # connections_tailwaters = supernetwork_values[4]

    waterbody_type_specified = False
    if break_network_at_waterbodies:
        # Read waterbody parameters
        waterbodies_df = nhd_io.read_waterbody_df(
            waterbody_parameters, {"level_pool": wbody_conn.values()}
        )

        # Remove duplicate lake_ids and rows
        waterbodies_df = (
            waterbodies_df.reset_index()
            .drop_duplicates(subset="lake_id")
            .set_index("lake_id")
        )

        # Declare empty dataframe
        waterbody_types_df = pd.DataFrame()

        wb_params_level_pool = waterbody_parameters.get(
            'level_pool', defaultdict(list)
        )
        
        # Check if hybrid-usgs or hybrid-usace reservoir DA is set to True
        reservoir_da = data_assimilation_parameters.get(
            'reservoir_da', 
            None
        )
        if reservoir_da:
            usgs_hybrid  = reservoir_da.get(
                'reservoir_persistence_usgs', 
                False
            )
            usace_hybrid = reservoir_da.get(
                'reservoir_persistence_usace', 
                False
            )
            param_file   = reservoir_da.get(
                'gage_lakeID_crosswalk_file',
                None
            )
        else:
            param_file = None
            usace_hybrid = False
            usgs_hybrid = False
            
        # check if RFC-type reservoirs are set to true
        rfc_params = waterbody_parameters.get('rfc')
        if rfc_params:
            rfc_forecast = rfc_params.get(
                'reservoir_rfc_forecasts',
                False
            )
            param_file = rfc_params.get('reservoir_parameter_file',None)
        else:
            rfc_forecast = False

        if param_file and reservoir_da:
            waterbody_type_specified = True
            waterbody_types_df = nhd_io.read_reservoir_parameter_file(
                param_file,
                usgs_hybrid,
                usace_hybrid,
                rfc_forecast,
                wb_params_level_pool.get("level_pool_waterbody_id", 'lake_id'),
                wbody_conn.values(),
            )
        else:
            waterbody_type_specified = True
            waterbody_types_df = pd.DataFrame(data = 1, index = waterbodies_df.index, columns = ['reservoir_type'])

    else:
        # Declare empty dataframes
        waterbody_types_df = pd.DataFrame()
        waterbodies_df = pd.DataFrame()

    # build diffusive domain data and edit MC domain data for hybrid simulation
    if diffusive_domain:
        
        rconn = nhd_network.reverse_network(connections)

        for tw in diffusive_domain:
            
            # ===== build diffusive network data objects ==== 
            diffusive_network_data[tw] = {}
            
            # add diffusive domain segments
            diffusive_network_data[tw]['mainstem_segs'] = diffusive_domain[tw]

            # diffusive domain tributary segments
            trib_segs = []
            for s in diffusive_domain[tw]:
                us_list = rconn[s]
                for u in us_list:
                    if u not in diffusive_domain[tw]:
                        trib_segs.append(u)
            diffusive_network_data[tw]['tributary_segments'] = trib_segs
            
            # diffusive domain connections object
            diffusive_network_data[tw]['connections'] = {k: connections[k] for k in (diffusive_domain[tw] + trib_segs)}
            
            # diffusive domain reaches and upstream connections. 
            # break network at tributary segments
            _, reaches, rconn_diff = nnu.organize_independent_networks(
                diffusive_network_data[tw]['connections'],
                set(trib_segs),
            )
            diffusive_network_data[tw]['rconn'] = rconn_diff
            diffusive_network_data[tw]['reaches'] = reaches[tw]
            
            # RouteLink parameters
            diffusive_network_data[tw]['param_df'] = param_df.filter(
                (diffusive_domain[tw] + trib_segs),
                axis = 0,
            )
        
            # ==== remove diffusive domain segs from MC domain ====
        
            # drop indices from param_df
            param_df = param_df.drop(diffusive_domain[tw])
        
            # remove keys from connections dictionary
            for s in diffusive_domain[tw]:
                connections.pop(s)
            
            # update downstream connections of trib segs
            for us in trib_segs:
                connections[us] = []
    
    # STEP 2: Identify Independent Networks and Reaches by Network
    
    start_time = time.time()

    LOG.info("organizing connections into reaches ...")

    network_break_segments = set()
    if break_network_at_waterbodies:
        network_break_segments = network_break_segments.union(wbody_conn.values())
    if break_network_at_gages:
        network_break_segments = network_break_segments.union(gages['gages'].keys())
        
    independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections,
        network_break_segments,
    )
    
    LOG.debug("reach organization complete in %s seconds." % (time.time() - start_time))
    
    if preprocessing_parameters.get('preprocess_only', False):

        LOG.debug("saving preprocessed network data to disk for future use")
        # todo: consider a better default than None
        destination_folder = preprocessing_parameters.get('preprocess_output_folder', None)
        if destination_folder:

            output_filename = preprocessing_parameters.get(
                'preprocess_output_filename', 
                'preprocess_output'
            )

            outputs = {}
            outputs.update(
                {'connections': connections,
                 'param_df': param_df,
                 'wbody_conn': wbody_conn,
                 'waterbodies_df': waterbodies_df,
                 'waterbody_types_df': waterbody_types_df,
                 'break_network_at_waterbodies': break_network_at_waterbodies,
                 'waterbody_type_specified': waterbody_type_specified,
                 'independent_networks': independent_networks,
                 'reaches_bytw': reaches_bytw,
                 'rconn': rconn,
                 'link_gage_df': pd.DataFrame.from_dict(gages),
                 'diffusive_network_data': diffusive_network_data,
                 'topobathy_data': topobathy_data,
                }
            )
            try:
                np.save(
                    pathlib.Path(destination_folder).joinpath(output_filename),
                    outputs
                )
            except:
                LOG.critical('Canonot find %s. Aborting preprocessing routine' % pathlib.Path(destination_folder))
                quit()
                
            LOG.debug(
                "writing preprocessed network data to %s"\
                % pathlib.Path(destination_folder).joinpath(output_filename + '.npy'))
            LOG.critical(
                "Preprocessed network data written to %s aborting preprocessing sequence" \
                % pathlib.Path(destination_folder).joinpath(output_filename + '.npy'))
            quit()

        else:
            LOG.critical(
                "No destination folder specified for preprocessing. Please specify preprocess_output_folder in configuration file. Aborting preprocessing routine"
            )
            quit()

    return (
        connections,
        param_df,
        wbody_conn,
        waterbodies_df,
        waterbody_types_df,
        break_network_at_waterbodies,  # Could this be inferred from the wbody_conn or waterbodies_df  # Could this be inferred from the wbody_conn or waterbodies_df? Consider making this name less about the network and more about the reservoir simulation.
        waterbody_type_specified,  # Seems like this could be inferred from waterbody_types_df...
        independent_networks,
        reaches_bytw,
        rconn,
        pd.DataFrame.from_dict(gages),
        diffusive_network_data,
        topobathy_data,
    )

def unpack_nwm_preprocess_data(preprocessing_parameters):
    
    preprocess_filepath = preprocessing_parameters.get('preprocess_source_file',None)
    if preprocess_filepath:
        try:
            inputs = np.load(pathlib.Path(preprocess_filepath),allow_pickle='TRUE').item()
        except:
            LOG.critical('Canonot find %s' % pathlib.Path(preprocess_filepath))
            quit()
              
        connections = inputs.get('connections',None)            
        param_df = inputs.get('param_df',None)
        wbody_conn = inputs.get('wbody_conn',None)
        waterbodies_df = inputs.get('waterbodies_df',None)
        waterbody_types_df = inputs.get('waterbody_types_df',None)
        break_network_at_waterbodies = inputs.get('break_network_at_waterbodies',None)
        waterbody_type_specified = inputs.get('waterbody_type_specified',None)
        independent_networks = inputs.get('independent_networks',None)
        reaches_bytw = inputs.get('reaches_bytw',None)
        rconn = inputs.get('rconn',None)
        gages = inputs.get('link_gage_df',None)
        diffusive_network_data = inputs.get('diffusive_network_data',None)
        topobathy_data = inputs.get('topobathy_data',None)
                
    else:
        LOG.critical("use_preprocessed_data = True, but no preprocess_source_file is specified. Aborting the simulation.")
        quit()
                         
    return (
        connections,
        param_df,
        wbody_conn,
        waterbodies_df,
        waterbody_types_df,
        break_network_at_waterbodies,
        waterbody_type_specified,
        independent_networks,
        reaches_bytw,
        rconn,
        gages,
        diffusive_network_data,
        topobathy_data,
    )


def nwm_initial_warmstate_preprocess(
    break_network_at_waterbodies,
    restart_parameters,
    data_assimilation_parameters,
    segment_index,
    waterbodies_df,
):

    '''
    Assemble model initial condition data:
        - waterbody inital states (outflow and pool elevation)
        - channel initial states (flow and depth)
        - initial time
        - most recent gage observations
    Additionally, a dictionary of data assimilation parameters is assembled.
    
    Arguments
    ---------
    - break_network_at_waterbodies (bool): If True, waterbody initial states will
                                           be appended to the waterbody parameter
                                           dataframe. If False, waterbodies will
                                           not be simulated and the waterbody
                                           parameter datataframe wil not be changed
    
    - restart_parameters           (dict): User-input simulation restart 
                                           parameters
    
    - data_assimilation_parameters (dict): User-input data assimilation 
                                           parameters
    
    - segment_index        (Pandas Index): All segment IDs in the simulation 
                                           doamin
    
    - waterbodies_df   (Pandas DataFrame): Waterbody parameters   
    
    Returns
    -------
    - waterbodies_df (Pandas DataFrame): Waterbody parameters with initial
                                         states (outflow and pool elevation)
    
    - q0             (Pandas DataFrame): Initial flow and depth states for each
                                         segment in the model domain
    
    - t0                     (datetime): Datetime of the model initialization
    
    - lastobs_df     (Pandas DataFrame): Last gage observations data for DA
    
    - da_parameter_dict          (dict): Data assimilation parameters
    
    Notes
    -----
    - I don't know if it makes sense to create the da_parameter_dict, here. 
      Consider moving the creation to another location..
    '''

    #----------------------------------------------------------------------------
    # Assemble waterbody initial states (outflow and pool elevation
    #----------------------------------------------------------------------------
    
    if break_network_at_waterbodies:

        start_time = time.time()
        LOG.info("setting waterbody initial states ...")

        # if a lite restart file is provided, read initial states from it.
        if restart_parameters.get("lite_waterbody_restart_file", None):
            
            waterbodies_initial_states_df, _ = nhd_io.read_lite_restart(
                restart_parameters['lite_waterbody_restart_file']
            )
            
        # read waterbody initial states from WRF-Hydro type restart file
        elif restart_parameters.get("wrf_hydro_waterbody_restart_file", None):
            waterbodies_initial_states_df = nhd_io.get_reservoir_restart_from_wrf_hydro(
                restart_parameters["wrf_hydro_waterbody_restart_file"],
                restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file"],
                restart_parameters.get("wrf_hydro_waterbody_ID_crosswalk_file_field_name", 'lake_id'),
                restart_parameters["wrf_hydro_waterbody_crosswalk_filter_file"],
                restart_parameters.get(
                    "wrf_hydro_waterbody_crosswalk_filter_file_field_name",
                    'NHDWaterbodyComID'
                ),
            )
        
        # if no restart file is provided, default initial states
        else:
            # TODO: Consider adding option to read cold state from route-link file
            waterbodies_initial_ds_flow_const = 0.0
            waterbodies_initial_depth_const = -1e9
            # Set initial states from cold-state
            waterbodies_initial_states_df = pd.DataFrame(
                0,
                index=waterbodies_df.index,
                columns=[
                    "qd0",
                    "h0",
                ],
                dtype="float32",
            )
            # TODO: This assignment could probably by done in the above call
            waterbodies_initial_states_df["qd0"] = waterbodies_initial_ds_flow_const
            waterbodies_initial_states_df["h0"] = waterbodies_initial_depth_const
            waterbodies_initial_states_df["index"] = range(
                len(waterbodies_initial_states_df)
            )

        waterbodies_df = pd.merge(
            waterbodies_df, waterbodies_initial_states_df, on="lake_id"
        )

        LOG.debug(
            "waterbody initial states complete in %s seconds."\
            % (time.time() - start_time))
        start_time = time.time()

    #----------------------------------------------------------------------------
    # Assemble channel initial states (flow and depth)
    # also establish simulation initialization timestamp
    #----------------------------------------------------------------------------
    
    start_time = time.time()
    LOG.info("setting channel initial states ...")

    # if lite restart file is provided, the read channel initial states from it
    if restart_parameters.get("lite_channel_restart_file", None):
        
        q0, t0 = nhd_io.read_lite_restart(
            restart_parameters['lite_channel_restart_file']
        )
    
    # build initial states from user-provided restart parameters
    else:
        q0 = nnu.build_channel_initial_state(restart_parameters, segment_index)

        # get initialization time from restart file
        if restart_parameters.get("wrf_hydro_channel_restart_file", None):
            channel_initial_states_file = restart_parameters[
                "wrf_hydro_channel_restart_file"
            ]
            t0_str = nhd_io.get_param_str(
                channel_initial_states_file, 
                "Restart_Time"
            )
        else:
            t0_str = "2015-08-16_00:00:00"

        # convert timestamp from string to datetime
        t0 = datetime.strptime(t0_str, "%Y-%m-%d_%H:%M:%S")
        
    # get initial time from user inputs
    if restart_parameters.get("start_datetime", None):
        t0_str = restart_parameters.get("start_datetime")
        
        def _try_parsing_date(text):
            for fmt in (
                "%Y-%m-%d_%H:%M", 
                "%Y-%m-%d_%H:%M:%S", 
                "%Y-%m-%d %H:%M", 
                "%Y-%m-%d %H:%M:%S", 
                "%Y/%m/%d %H:%M", 
                "%Y/%m/%d %H:%M:%S"
            ):
                try:
                    return datetime.strptime(text, fmt)
                except ValueError:
                    pass
            LOG.error('No valid date format found for start_datetime input. Please use format YYYY-MM-DD_HH:MM')
            quit()
            
        t0 = _try_parsing_date(t0_str)
    else:
        if t0_str == "2015-08-16_00:00:00":
            LOG.info('No user-input start_datetime and no restart file, start time arbitrarily 2015-08-16_00:00:00')
        else:
            LOG.info('No user-specified start_datetime, continuing with start time from restart file: %s', t0_str)
        
    #----------------------------------------------------------------------------
    # Assemble streamflow DA lastobs data
    #----------------------------------------------------------------------------
    
    lastobs_df, da_parameter_dict = nnu.build_data_assimilation_lastobs(
        data_assimilation_parameters
    )

    LOG.debug(
        "channel initial states complete in %s seconds."\
        % (time.time() - start_time)
    )
    start_time = time.time()

    return waterbodies_df, q0, t0, lastobs_df, da_parameter_dict
    # TODO: This returns a full dataframe (waterbodies_df) with the
    # merged initial states for waterbodies, but only the
    # initial state values (q0; not merged with the channel properties)
    # for the channels --
    # That is because that is how they are used downstream. Need to
    # trace that back and decide if there is one of those two ways
    # that is optimal and make both returns that way.


def nwm_forcing_preprocess(
    run,
    forcing_parameters,
    da_run,
    data_assimilation_parameters,
    break_network_at_waterbodies,
    segment_index,
    link_gage_df,
    lastobs_index,
    cpu_pool,
    t0,
):
    """
    Assemble model forcings. Forcings include hydrological lateral inflows (qlats)
    and observations for streamflow and/or reservoir data assimilation schemes
    
    Aguments
    --------
    - run                          (dict): List of forcing files pertaining
                                           to a single run-set
    
    - forcing_parameters           (dict): User-input simulation forcing 
                                           parameters
    
    - da_run                       (dict): Lists of TimeSlice files for a
                                           single run-set
    
    - data_assimilation_parameters (dict): User-input parameters controlling
                                           data assimilation routines
    
    - break_network_at_waterbodies (bool): ????
    
    - segment_index                    (): ????
    
    - link_gage_df     (Pandas DataFrame): Crosswalking between segment ID and
                                           USGS gage IDs in the model domain
    
    - lastobs_index        (Pandas Index): ????
    
    - cpu_pool                      (int): Number of CPUs in the process-parall
                                           pool
    
    - t0                       (datetime): Timestamp of the simualtion initial 
                                           condition

    Returns
    -------
    - qlats_df           (Pandas DataFrame): Lateral inflow data, indexed by 
                                             segment ID
                                             
    - usgs_df            (Pandas DataFrame): Streamflow DA observations, indexed by
                                             segment ID
                                             
    - reservoir_usgs_df  (Pandas DataFrame): Reservoir DA observations 
                                             (USGS persistence reservoirs), indexed 
                                             by waterbody ID
                                             
    - reservoir_usace_df (Pandas DataFrame): Reservoir DA observations 
                                             (USACE persistence reservoirs), indexed 
                                             by waterbody ID
    
    Notes
    -----
    - This function adds new keys to the `run` input dictionary
    
    """

    # Unpack user-specified forcing parameters
    dt                           = forcing_parameters.get("dt", None)
    qts_subdivisions             = forcing_parameters.get("qts_subdivisions", None)
    qlat_input_folder            = forcing_parameters.get("qlat_input_folder", None)
    qlat_file_index_col          = forcing_parameters.get("qlat_file_index_col", "feature_id")
    qlat_file_value_col          = forcing_parameters.get("qlat_file_value_col", "q_lateral")
    qlat_file_gw_bucket_flux_col = forcing_parameters.get("qlat_file_gw_bucket_flux_col", "qBucket")
    qlat_file_terrain_runoff_col = forcing_parameters.get("qlat_file_terrain_runoff_col", "qSfcLatRunoff")

    # TODO: find a better way to deal with these defaults and overrides.
    run["t0"]                           = run.get("t0", t0)
    run["nts"]                          = run.get("nts")
    run["dt"]                           = run.get("dt", dt)
    run["qts_subdivisions"]             = run.get("qts_subdivisions", qts_subdivisions)
    run["qlat_input_folder"]            = run.get("qlat_input_folder", qlat_input_folder)
    run["qlat_file_index_col"]          = run.get("qlat_file_index_col", qlat_file_index_col)
    run["qlat_file_value_col"]          = run.get("qlat_file_value_col", qlat_file_value_col)
    run["qlat_file_gw_bucket_flux_col"] = run.get("qlat_file_gw_bucket_flux_col", qlat_file_gw_bucket_flux_col)
    run["qlat_file_terrain_runoff_col"] = run.get("qlat_file_terrain_runoff_col", qlat_file_terrain_runoff_col)

    #---------------------------------------------------------------------------
    # Assemble lateral inflow data
    #---------------------------------------------------------------------------
    
    start_time = time.time()
    LOG.info("Creating a DataFrame of lateral inflow forcings ...")

    qlats_df = nnu.build_qlateral_array(
        run,
        cpu_pool,
        segment_index,
    )

    LOG.debug(
        "lateral inflow DataFrame creation complete in %s seconds." \
        % (time.time() - start_time)
    )
    
    #---------------------------------------------------------------------------
    # Assemble streamflow DA observation data
    #---------------------------------------------------------------------------
    
    # isolate user-input parameters for streamflow data assimilation
    streamflow_da_parameters = data_assimilation_parameters.get(
        'streamflow_da', 
        None
    )
    
    # determine if user explictly requests streamflow DA
    nudging = False
    if streamflow_da_parameters:
        nudging = streamflow_da_parameters.get('streamflow_nudging', False)
        
    # if user requested nudging and a specified a USGS TimeSlice directory, 
    # then build and return USGS dataframe
    if nudging and da_run['usgs_timeslice_files']:
        
        start_time = time.time()
        LOG.info(
            "Creating a DataFrame of USGS gage observations for streamflow DA ..."
        )
        
        usgs_timeslices_folder = data_assimilation_parameters.get(
                                    "usgs_timeslices_folder", 
                                    None
                                )
        lastobs_file           = streamflow_da_parameters.get(
                                    "wrf_hydro_lastobs_file", 
                                    None
                                )
        lastobs_start          = data_assimilation_parameters.get(
            "wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time",
                                    0
                                )
        lastobs_type           = data_assimilation_parameters.get(
                                    "wrf_lastobs_type", 
                                    "error-based"
                                )
        crosswalk_file         = streamflow_da_parameters.get(
                                    "gage_segID_crosswalk_file", 
                                    None
                                )
        crosswalk_gage_field   = streamflow_da_parameters.get(
                                    'crosswalk_gage_field',
                                    'gages'
                                )
        crosswalk_segID_field  = streamflow_da_parameters.get(
                                    'crosswalk_segID_field',
                                    'link'
                                )
        da_decay_coefficient   = data_assimilation_parameters.get(
                                    "da_decay_coefficient",
                                    120
                                )
        qc_threshold           = data_assimilation_parameters.get(
                                    "qc_threshold",
                                    1
                                )
        interpolation_limit    = data_assimilation_parameters.get(
                                    "interpolation_limit_min",
                                    59
                                )
        
        # TODO: join timeslice folder and files into complete path upstream
        usgs_timeslices_folder = pathlib.Path(usgs_timeslices_folder)
        usgs_files = [usgs_timeslices_folder.joinpath(f) for f in 
                      da_run['usgs_timeslice_files']]
        
        if usgs_files:
            usgs_df = (
                nhd_io.get_obs_from_timeslices(
                    crosswalk_file,
                    crosswalk_gage_field,
                    crosswalk_segID_field,
                    usgs_files,
                    qc_threshold,
                    interpolation_limit,
                    run["dt"],
                    run["t0"],
                    cpu_pool
                ).
                loc[link_gage_df.index]
            )
            
        else:
            usgs_df = pd.DataFrame()
        
        LOG.debug(
            "Streamflow DA USGS observation DataFrame creation complete in %s seconds." \
            % (time.time() - start_time)
        )
                 
    else: 
        usgs_df = pd.DataFrame()
    
    #---------------------------------------------------------------------------
    # Assemble reservoir data assimilation observation data
    #---------------------------------------------------------------------------
    
    # isolate user-input parameters for reservoir data assimilation
    reservoir_da_parameters = data_assimilation_parameters.get(
        'reservoir_da', 
        None
    )
    
    # check if user explictly requests USGS and/or USACE reservoir DA
    usgs_persistence  = False
    usace_persistence = False
    if reservoir_da_parameters:
        usgs_persistence  = reservoir_da_parameters.get(
            'reservoir_persistence_usgs', 
            False
        )
        usace_persistence = reservoir_da_parameters.get(
            'reservoir_persistence_usace', 
            False
        )
        
    #---------------------------------------------
    # observations for USGS reservoir DA
    #---------------------------------------------
    
    # create crosswalking dataframes and identify USGS lakeIDs in the model domain
    if usgs_persistence:
        
        crosswalk_file = reservoir_da_parameters["gage_lakeID_crosswalk_file"]
        crosswalk_file = pathlib.Path(crosswalk_file)
        
        with xr.open_dataset(crosswalk_file) as ds:
            gage_lake_dict = {
                'gage': ds.usgs_gage_id.values,
                'usgs_lake_id': ds.usgs_lake_id.values
            }
            gage_lake_df = (
                pd.DataFrame(data = gage_lake_dict).
                set_index('gage')
            )
        
        # build dataframe that crosswalks gageIDs to segmentIDs
        gage_link_df = (
            link_gage_df['gages'].
            reset_index().
            set_index(['gages'])
        )
        
        # extract the USGS lakeIDs within the model domain
        usgs_lakes_in_domain = (
            gage_lake_df.join(gage_link_df, how = 'inner')['usgs_lake_id'].
            to_numpy()
        )
    
    # if USGS TimeSlices have already been opened and assembled for 
    # streamflow DA, then simply take TimeSlice observations from `usgs_df`,
    # no need to open TimeSlice files again.
    if usgs_persistence and nudging and usgs_df.empty == False:
        
        start_time = time.time()
        LOG.info(
            "Creating a DataFrame of USGS gage observations for Reservoir DA ..."
        )
        
        # build dataframe that crosswalks segmentIDs to lakeIDs
        link_lake_df = (
            gage_lake_df.
            join(gage_link_df, how = 'inner').
            reset_index().set_index('index').
            drop(['level_0'], axis = 1)
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
        
        LOG.debug(
            "Reservoir DA USGS observation DataFrame creation complete in %s seconds." \
            % (time.time() - start_time)
        )
                    
    elif usgs_persistence and not nudging:
        
        start_time = time.time()
        LOG.info("Creating a DataFrame of USGS gage observations for Reservoir DA ...")
        
        usgs_timeslices_folder = data_assimilation_parameters.get(
                                   "usgs_timeslices_folder",
                                   None
                                )
        crosswalk_file         = reservoir_da_parameters.get(
                                    "gage_lakeID_crosswalk_file", 
                                    None
                                )
        crosswalk_gage_field   = streamflow_da_parameters.get(
                                    'crosswalk_usgs_gage_field',
                                    'usgs_gage_id'
                                )
        crosswalk_lakeID_field  = streamflow_da_parameters.get(
                                    'crosswalk_usgs_lakeID_field',
                                    'usgs_lake_id'
                                )
        qc_threshold            = data_assimilation_parameters.get(
                                    "qc_threshold",
                                    1
                                )
        interpolation_limit     = data_assimilation_parameters.get(
                                    "interpolation_limit_min",
                                    59
                                )
        
        # TODO: join timeslice folder and files into complete path upstream in workflow
        usgs_timeslices_folder = pathlib.Path(usgs_timeslices_folder)
        usgs_files = [usgs_timeslices_folder.joinpath(f) for f in
                      da_run['usgs_timeslice_files']]
                
        if usgs_files:
            
            reservoir_usgs_df = nhd_io.get_obs_from_timeslices(
                crosswalk_file,
                crosswalk_gage_field,
                crosswalk_lakeID_field,
                usgs_files,
                qc_threshold,
                interpolation_limit,
                900,                      # 15 minutes, as secs
                run["t0"],
                cpu_pool
            )
            
        else:
            reservoir_usgs_df = pd.DataFrame()
            
        LOG.debug(
            "Reservoir DA USGS observation DataFrame creation complete in %s seconds." \
            % (time.time() - start_time)
        )  
        
    else:
        reservoir_usgs_df = pd.DataFrame()
        
    #---------------------------------------------
    # observations for USACE reservoir DA
    #---------------------------------------------  
    
    if usace_persistence:
        start_time = time.time()
        LOG.info("Creating a DataFrame of USACE gage observations for Reservoir DA ...")
        
        usace_timeslices_folder = data_assimilation_parameters.get(
                                   "usace_timeslices_folder",
                                   None
                                )
        crosswalk_file          = reservoir_da_parameters.get(
                                    "gage_lakeID_crosswalk_file", 
                                    None
                                )
        crosswalk_gage_field    = streamflow_da_parameters.get(
                                    'crosswalk_usace_gage_field',
                                    'usace_gage_id'
                                )
        crosswalk_lakeID_field   = streamflow_da_parameters.get(
                                    'crosswalk_usace_lakeID_field',
                                    'usace_lake_id'
                                )
        qc_threshold            = data_assimilation_parameters.get(
                                    "qc_threshold",
                                    1
                                )
        interpolation_limit     = data_assimilation_parameters.get(
                                    "interpolation_limit_min",
                                    59
                                )
        
        # TODO: join timeslice folder and files into complete path upstream in workflow
        usace_timeslices_folder = pathlib.Path(usace_timeslices_folder)
        usace_files = [usace_timeslices_folder.joinpath(f) for f in da_run['usace_timeslice_files']]
                
        if usace_files:
            reservoir_usace_df = nhd_io.get_obs_from_timeslices(
                crosswalk_file,
                crosswalk_gage_field,
                crosswalk_lakeID_field,
                usace_files,
                qc_threshold,
                interpolation_limit,
                900,                      # 15 minutes, as secs
                run["t0"],
                cpu_pool
            )
        else:
            reservoir_usace_df = pd.DataFrame()
            
        LOG.debug(
            "Reservoir DA USACE observation DataFrame creation complete in %s seconds." \
            % (time.time() - start_time)
        )
        
    else:
        reservoir_usace_df = pd.DataFrame()
    
    #---------------------------------------------------------------------------
    # Assemble coastal coupling data [WIP]
    
    coastal_boundary_elev = forcing_parameters.get("coastal_boundary_elev_data", None)
    coastal_ncdf = forcing_parameters.get("coastal_ncdf", None)

    if coastal_boundary_elev:
        LOG.info("creating coastal dataframe ...")
        coastal_df = nhd_io.build_coastal_dataframe(coastal_boundary_elev)

    if coastal_ncdf:
        LOG.info("creating coastal ncdf dataframe ...")
        coastal_ncdf_df = nhd_io.build_coastal_ncdf_dataframe(coastal_ncdf)

    #---------------------------------------------------------------------------
    # Trim the time-extent of the streamflow_da usgs_df
    # what happens if there are timeslice files missing on the front-end? 
    # if the first column is some timestamp greater than t0, then this will throw
    # an error. Need to think through this more. 
    if not usgs_df.empty:
        usgs_df = usgs_df.loc[:,t0:]
    
    return qlats_df, usgs_df, reservoir_usgs_df, reservoir_usace_df
