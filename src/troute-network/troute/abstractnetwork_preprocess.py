import json
import pathlib
from functools import partial, reduce
from itertools import chain
from datetime import datetime, timedelta
from collections import defaultdict, deque
import logging
import os
import time

import pandas as pd
import numpy as np
import netCDF4
from joblib import delayed, Parallel
import pyarrow as pa
import pyarrow.parquet as pq

import troute.nhd_io as nhd_io
import troute.nhd_network as nhd_network
import troute.nhd_network_utilities_v02 as nnu


LOG = logging.getLogger('')

def build_diffusive_domain(
    compute_parameters,
    param_df,
    connections,
):
    
    hybrid_params = compute_parameters.get("hybrid_parameters", False)
    if hybrid_params:
        # switch parameters
        # if run_hybrid = False, run MC only
        # if run_hybrid = True, if use_topobathy = False, run MC+diffusive on RouteLink.nc
        #    "      "      "  , if use_topobathy = True,  if run_refactored_network = False, run MC+diffusive on original hydrofabric
        #    "      "      "  , if use_topobathy = True,  if run_refactored_network = True,  run MC+diffusive on refactored hydrofabric
        run_hybrid             = hybrid_params.get('run_hybrid_routing', False)
        use_topobathy          = hybrid_params.get('use_natl_xsections', False) 
        run_refactored         = hybrid_params.get('run_refactored_network', False)
        
        # file path parameters of non-refactored hydrofabric defined by RouteLink.nc
        domain_file            = hybrid_params.get("diffusive_domain",   None)
        topobathy_file         = hybrid_params.get("topobathy_domain",   None)
        
        # file path parameters of refactored hydrofabric for diffusive wave channel routing 
        refactored_domain_file    = hybrid_params.get("refactored_domain", None)
        refactored_topobathy_file = hybrid_params.get("refactored_topobathy_domain", None)
        #-------------------------------------------------------------------------
        # for non-refactored hydofabric defined by RouteLink.nc
        # TODO: By default, make diffusive available for both non-refactored and refactored hydrofabric for now. Place a switch in the future. 
        if run_hybrid and domain_file:
            
            LOG.info('reading diffusive domain extent for MC/Diffusive hybrid simulation')
            
            # read diffusive domain dictionary from yaml or json
            diffusive_domain = nhd_io.read_diffusive_domain(domain_file)           
            
            if use_topobathy and topobathy_file:
                
                LOG.debug('Natural cross section data on original hydrofabric are provided.')
                
                # read topobathy domain netcdf file, set index to 'comid'
                # TODO: replace 'link' with a user-specified indexing variable name.
                # ... if for whatever reason there is not a `link` variable in the 
                # ... dataframe returned from read_netcdf, then the code would break here.
                topobathy_df = (nhd_io.read_netcdf(topobathy_file).set_index('link'))
                
                # TODO: Request GID make comID variable an integer in their product, so
                # we do not need to change variable types, here.
                topobathy_df.index = topobathy_df.index.astype(int)
                
            else:
                topobathy_df = pd.DataFrame()
                LOG.debug('No natural cross section topobathy data provided. Hybrid simualtion will run on compound trapezoidal geometry.')
             
            # initialize a dictionary to hold network data for each of the diffusive domains
            diffusive_network_data = {}
        
        else:
            diffusive_domain       = None
            diffusive_network_data = None
            topobathy_df           = pd.DataFrame()
            LOG.info('No diffusive domain file specified in configuration file. This is an MC-only simulation')
        unrefactored_topobathy_df  = pd.DataFrame()    
        #-------------------------------------------------------------------------
        # for refactored hydofabric 
        if run_hybrid and run_refactored and refactored_domain_file:
            
            LOG.info('reading refactored diffusive domain extent for MC/Diffusive hybrid simulation')
            
            # read diffusive domain dictionary from yaml or json
            refactored_diffusive_domain = nhd_io.read_diffusive_domain(refactored_domain_file)           
            
            if use_topobathy and refactored_topobathy_file:
                
                LOG.debug('Natural cross section data of refactored hydrofabric are provided.')
                
                # read topobathy domain netcdf file, set index to 'comid'
                # TODO: replace 'link' with a user-specified indexing variable name.
                # ... if for whatever reason there is not a `link` variable in the 
                # ... dataframe returned from read_netcdf, then the code would break here.
                topobathy_df = (nhd_io.read_netcdf(refactored_topobathy_file).set_index('link'))

                # unrefactored_topobaty_data is passed to diffusive kernel to provide thalweg elevation of unrefactored topobathy 
                # for crosswalking water elevations between non-refactored and refactored hydrofabrics. 
                unrefactored_topobathy_df       = (nhd_io.read_netcdf(topobathy_file).set_index('link'))
                unrefactored_topobathy_df.index = unrefactored_topobathy_df.index.astype(int)
                
            else:
                topobathy_df               = pd.DataFrame()
                LOG.debug('No natural cross section topobathy data of refactored hydrofabric provided. Hybrid simualtion will run on compound trapezoidal geometry.')
             
            # initialize a dictionary to hold network data for each of the diffusive domains
            refactored_diffusive_network_data = {}
        
        else:
            refactored_diffusive_domain       = None
            refactored_diffusive_network_data = None
            refactored_reaches                = {}
            LOG.info('No refactored diffusive domain file specified in configuration file. This is an MC-only simulation')     
   
    else:
        diffusive_domain                  = None
        diffusive_network_data            = None
        topobathy_df                      = pd.DataFrame()
        unrefactored_topobathy_df         = pd.DataFrame() 
        refactored_diffusive_domain       = None
        refactored_diffusive_network_data = None   
        refactored_reaches                = {}
        LOG.info('No hybrid parameters specified in configuration file. This is an MC-only simulation')
    
 #============================================================================
    # build diffusive domain data and edit MC domain data for hybrid simulation
    
    #
    if diffusive_domain:
        rconn_diff0 = nhd_network.reverse_network(connections)
        refactored_reaches = {}
        
        for tw in diffusive_domain:
            mainstem_segs = diffusive_domain[tw]['links']
            # we want mainstem_segs start at a mainstem link right after the upstream boundary mainstem link, which is
            # in turn not under any waterbody. This boundary mainstem link should be turned into a tributary segment.
            upstream_boundary_mainstem_link = diffusive_domain[tw]['upstream_boundary_link_mainstem']         
            if upstream_boundary_mainstem_link[0] in mainstem_segs:
                mainstem_segs.remove(upstream_boundary_mainstem_link[0])
            
            # ===== build diffusive network data objects ==== 
            diffusive_network_data[tw] = {}

            # add diffusive domain segments
            diffusive_network_data[tw]['mainstem_segs'] =  mainstem_segs

            # diffusive domain tributary segments
            trib_segs = []
            
            for seg in mainstem_segs:
                us_list = rconn_diff0[seg]
                for u in us_list:
                    if u not in mainstem_segs:
                        trib_segs.append(u) 

            diffusive_network_data[tw]['tributary_segments'] = trib_segs
            # diffusive domain connections object
            diffusive_network_data[tw]['connections'] = {k: connections[k] for k in (mainstem_segs + trib_segs)}       

            # diffusive domain reaches and upstream connections. 
            # break network at tributary segments
            _, reaches, rconn_diff = nnu.organize_independent_networks(
                diffusive_network_data[tw]['connections'],
                set(trib_segs),
                set(),
            )
            
            diffusive_network_data[tw]['rconn'] = rconn_diff
            diffusive_network_data[tw]['reaches'] = reaches[tw]

            # RouteLink parameters
            diffusive_network_data[tw]['param_df'] = param_df.filter(
                (mainstem_segs + trib_segs),
                axis = 0,
            )
            diffusive_network_data[tw]['upstream_boundary_link'] = upstream_boundary_mainstem_link

            if refactored_diffusive_domain: 
                diffusive_parameters = {'geo_file_path': refactored_topobathy_file}
                refactored_connections = nnu.build_refac_connections(diffusive_parameters)

                # list of stream segments of a single refactored diffusive domain 
                refac_tw = refactored_diffusive_domain[tw]['refac_tw']
                rlinks_tw = refactored_diffusive_domain[tw]['rlinks']
                refactored_connections_tw = {}   

                # Subset a connection dictionary (upstream segment as key : downstream segments as values) from refactored_connections
                # for a single refactored diffusive domain defined by a current tw. 
                for k in rlinks_tw:
                    if k in refactored_connections.keys() and k != refac_tw:
                        refactored_connections_tw[k] = refactored_connections[k]

                refactored_diffusive_network_data[refac_tw] = {}                
                refactored_diffusive_network_data[refac_tw]['tributary_segments'] = trib_segs
                refactored_diffusive_network_data[refac_tw]['connections'] = refactored_connections_tw                 

                for k in trib_segs:
                    refactored_diffusive_network_data[refac_tw]['connections'][k]= [refactored_diffusive_domain[tw]['incoming_tribs'][k]]

                # diffusive domain reaches and upstream connections. 
                # break network at tributary segments
                _, refactored_reaches_batch, refactored_conn_diff = nnu.organize_independent_networks(
                                                            refactored_diffusive_network_data[refac_tw]['connections'],
                                                            set(trib_segs),
                                                            set(),
                                                            )

                refactored_reaches[refac_tw] = refactored_reaches_batch[refac_tw]
                refactored_diffusive_network_data[refac_tw]['mainstem_segs'] = refactored_diffusive_domain[tw]['rlinks']
                refactored_diffusive_network_data[refac_tw]['upstream_boundary_link'] = diffusive_network_data[tw]['upstream_boundary_link'] 
            else:
                refactored_reaches={}

            # ==== remove diffusive domain segs from MC domain ====        
            # drop indices from param_df
            param_df = param_df.drop(mainstem_segs)

            # remove keys from connections dictionary
            for s in mainstem_segs:
                connections.pop(s)

            # update downstream connections of trib segs
            for us in trib_segs:
                connections[us] = []
    
    return (
        param_df,
        connections,
        diffusive_network_data,
        topobathy_df,
        refactored_diffusive_domain,
        refactored_reaches,
        unrefactored_topobathy_df
    )

def create_independent_networks(
    waterbody_parameters, 
    connections, 
    wbody_conn, 
    gages = pd.DataFrame() #FIXME update default value when we update 'break_network_at_gages',
    ):

    LOG.info("organizing connections into reaches ...")
    start_time = time.time() 
    gage_break_segments = set()
    wbody_break_segments = set()
    
    break_network_at_waterbodies = waterbody_parameters.get(
        "break_network_at_waterbodies", False
    )
    
    # if streamflow DA, then break network at gages
    #TODO update to work with HYFeatures, need to determine how we'll do DA...
    break_network_at_gages = False
    
    if break_network_at_waterbodies:
        wbody_break_segments = wbody_break_segments.union(wbody_conn.values())
        
    if break_network_at_gages:
        gage_break_segments = gage_break_segments.union(gages['gages'].keys())
 
    independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections,
        wbody_break_segments,
        gage_break_segments,
    )
    
    LOG.debug("reach organization complete in %s seconds." % (time.time() - start_time))

    return independent_networks, reaches_bytw, rconn

def initial_warmstate_preprocess(
    break_network_at_waterbodies,
    restart_parameters,
    segment_index,
    waterbodies_df,
    ):

    '''
    Assemble model initial condition data:
        - waterbody inital states (outflow and pool elevation)
        - channel initial states (flow and depth)
        - initial time
    
    Arguments
    ---------
    - break_network_at_waterbodies (bool): If True, waterbody initial states will
                                           be appended to the waterbody parameter
                                           dataframe. If False, waterbodies will
                                           not be simulated and the waterbody
                                           parameter datataframe wil not be changed
    - restart_parameters           (dict): User-input simulation restart 
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
    
    Notes
    -----
    '''

    # generalize waterbody ID's to be used with any network
    index_id = waterbodies_df.index.names[0]

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
                restart_parameters.get("wrf_hydro_waterbody_ID_crosswalk_file_field_name", index_id),
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
            waterbodies_df, waterbodies_initial_states_df, on=index_id
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
        # FIXME: Change it for hyfeature!
        q0, t0 = nhd_io.read_lite_restart(
            restart_parameters['lite_channel_restart_file']
        )
        t0_str = None

    # when a restart file for hyfeature is provied, then read initial states from it.
    elif restart_parameters.get("hyfeature_channel_restart_file", None):        
        q0 = nnu.build_channel_initial_state(restart_parameters, segment_index)        
        channel_initial_states_file = restart_parameters["hyfeature_channel_restart_file"]
        df     = pd.read_csv(channel_initial_states_file)
        t0_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
        t0     = datetime.strptime(t0_str,"%Y-%m-%d_%H:%M:%S")

    # build initial states from user-provided restart parameters
    else:
        # FIXME: Change it for hyfeature!
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

    LOG.debug(
        "channel initial states complete in %s seconds."\
        % (time.time() - start_time)
    )
    start_time = time.time()

    return (
            waterbodies_df,
            q0, 
            t0,
            )
    # TODO: This returns a full dataframe (waterbodies_df) with the
    # merged initial states for waterbodies, but only the
    # initial state values (q0; not merged with the channel properties)
    # for the channels --
    # That is because that is how they are used downstream. Need to
    # trace that back and decide if there is one of those two ways
    # that is optimal and make both returns that way.

def build_forcing_sets(
    supernetwork_parameters,
    forcing_parameters,
    t0,
    ):

    run_sets           = forcing_parameters.get("qlat_forcing_sets", None)
    qlat_input_folder  = forcing_parameters.get("qlat_input_folder", None)
    nts                = forcing_parameters.get("nts", None)
    max_loop_size      = forcing_parameters.get("max_loop_size", 12)
    dt                 = forcing_parameters.get("dt", None)

    geo_file_type      = supernetwork_parameters.get('geo_file_type')

    try:
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        assert qlat_input_folder.is_dir() == True
    except TypeError:
        raise TypeError("Aborting simulation because no qlat_input_folder is specified in the forcing_parameters section of the .yaml control file.") from None
    except AssertionError:
        raise AssertionError("Aborting simulation because the qlat_input_folder:", qlat_input_folder,"does not exist. Please check the the nexus_input_folder variable is correctly entered in the .yaml control file") from None

    forcing_glob_filter = forcing_parameters.get("qlat_file_pattern_filter", "*.NEXOUT")

    if forcing_glob_filter=="nex-*":
        print("Reformating qlat nexus files as hourly binary files...")
        binary_folder = forcing_parameters.get('binary_nexus_file_folder', None)
        qlat_files = qlat_input_folder.glob(forcing_glob_filter)

        #Check that directory/files specified will work
        if not binary_folder:
            raise(RuntimeError("No output binary qlat folder supplied in config"))
        elif not os.path.exists(binary_folder):
            raise(RuntimeError("Output binary qlat folder supplied in config does not exist"))
        elif len(list(pathlib.Path(binary_folder).glob('*.parquet'))) != 0:
            raise(RuntimeError("Output binary qlat folder supplied in config is not empty (already contains '.parquet' files)"))

        #Add tnx for backwards compatability
        qlat_files_list = list(qlat_files) + list(qlat_input_folder.glob('tnx*.csv'))
        #Convert files to binary hourly files, reset nexus input information
        qlat_input_folder, forcing_glob_filter = nex_files_to_binary(qlat_files_list, binary_folder)
        forcing_parameters["qlat_input_folder"] = qlat_input_folder
        forcing_parameters["qlat_file_pattern_filter"] = forcing_glob_filter
        
    # TODO: Throw errors if insufficient input data are available
    if run_sets:        
        #FIXME: Change it for hyfeature
        '''
        # append final_timestamp variable to each set_list
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        for (s, _) in enumerate(run_sets):
            final_chrtout = qlat_input_folder.joinpath(run_sets[s]['qlat_files'
                    ][-1])
            final_timestamp_str = nhd_io.get_param_str(final_chrtout,
                    'model_output_valid_time')
            run_sets[s]['final_timestamp'] = \
                datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')
        '''  
    elif qlat_input_folder:        
        # Construct run_set dictionary from user-specified parameters

        # get the first and seconded files from an ordered list of all forcing files
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        all_files          = sorted(qlat_input_folder.glob(forcing_glob_filter))
        first_file         = all_files[0]
        second_file        = all_files[1]

        # Deduce the timeinterval of the forcing data from the output timestamps of the first
        # two ordered CHRTOUT files
        if geo_file_type=='HYFeaturesNetwork':
            df     = read_file(first_file)
            t1_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
            t1     = datetime.strptime(t1_str,"%Y-%m-%d_%H:%M:%S")
            df     = read_file(second_file)
            t2_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
            t2     = datetime.strptime(t2_str,"%Y-%m-%d_%H:%M:%S")
        elif geo_file_type=='NHDNetwork':
            t1 = nhd_io.get_param_str(first_file, "model_output_valid_time")
            t1 = datetime.strptime(t1, "%Y-%m-%d_%H:%M:%S")
            t2 = nhd_io.get_param_str(second_file, "model_output_valid_time")
            t2 = datetime.strptime(t2, "%Y-%m-%d_%H:%M:%S")
        
        dt_qlat_timedelta = t2 - t1
        dt_qlat = dt_qlat_timedelta.seconds

        # determine qts_subdivisions
        qts_subdivisions = dt_qlat / dt
        if dt_qlat % dt == 0:
            qts_subdivisions = int(dt_qlat / dt)
        # make sure that qts_subdivisions = dt_qlat / dt
        forcing_parameters['qts_subdivisions']= qts_subdivisions

        # the number of files required for the simulation
        nfiles = int(np.ceil(nts / qts_subdivisions))
        
        # list of forcing file datetimes
        #datetime_list = [t0 + dt_qlat_timedelta * (n + 1) for n in
        #                 range(nfiles)]
        # ** Correction ** Because qlat file at time t is constantly applied throughout [t, t+1],
        #               ** n + 1 should be replaced by n
        datetime_list = [t0 + dt_qlat_timedelta * (n) for n in
                         range(nfiles)]        
        datetime_list_str = [datetime.strftime(d, '%Y%m%d%H%M') for d in
                             datetime_list]

        # list of forcing files
        forcing_filename_list = [d_str + forcing_glob_filter[1:] for d_str in
                                 datetime_list_str]
        
        # check that all forcing files exist
        for f in forcing_filename_list:
            try:
                J = pathlib.Path(qlat_input_folder.joinpath(f))     
                assert J.is_file() == True
            except AssertionError:
                raise AssertionError("Aborting simulation because forcing file", J, "cannot be not found.") from None
                
        # build run sets list
        run_sets = []
        k = 0
        j = 0
        nts_accum = 0
        nts_last = 0
        while k < len(forcing_filename_list):
            run_sets.append({})

            if k + max_loop_size < len(forcing_filename_list):
                run_sets[j]['qlat_files'] = forcing_filename_list[k:k
                    + max_loop_size]
            else:
                run_sets[j]['qlat_files'] = forcing_filename_list[k:]

            nts_accum += len(run_sets[j]['qlat_files']) * qts_subdivisions
            if nts_accum <= nts:
                run_sets[j]['nts'] = int(len(run_sets[j]['qlat_files'])
                                         * qts_subdivisions)
            else:
                run_sets[j]['nts'] = int(nts - nts_last)

            final_qlat = qlat_input_folder.joinpath(run_sets[j]['qlat_files'][-1]) 
            if geo_file_type=='NHDNetwork':           
                final_timestamp_str = nhd_io.get_param_str(final_qlat,'model_output_valid_time')
            elif geo_file_type=='HYFeaturesNetwork':
                df = read_file(final_qlat)
                final_timestamp_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")           
            
            run_sets[j]['final_timestamp'] = \
                datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')

            nts_last = nts_accum
            k += max_loop_size
            j += 1

    return run_sets

def build_qlateral_array(
    run,
    cpu_pool,
    nexus_to_upstream_flowpath_dict,
    supernetwork_parameters, 
    segment_index=pd.Index([]),
):

    start_time = time.time()
    LOG.info("Creating a DataFrame of lateral inflow forcings ...")
    
    # TODO: set default/optional arguments
    qts_subdivisions = run.get("qts_subdivisions", 1)
    nts = run.get("nts", 1)
    qlat_input_folder = run.get("qlat_input_folder", None)
    qlat_input_file = run.get("qlat_input_file", None)

    geo_file_type = supernetwork_parameters.get('geo_file_type')

    if qlat_input_folder:
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        if "qlat_files" in run:
            qlat_files = run.get("qlat_files")
            qlat_files = [qlat_input_folder.joinpath(f) for f in qlat_files]
        elif "qlat_file_pattern_filter" in run:
            qlat_file_pattern_filter = run.get(
                "qlat_file_pattern_filter", "*CHRT_OUT*"
            )
            qlat_files = sorted(qlat_input_folder.glob(qlat_file_pattern_filter))

        qlat_file_index_col = run.get(
            "qlat_file_index_col", "feature_id"
        )
        qlat_file_value_col = run.get("qlat_file_value_col", "q_lateral")
        gw_bucket_col = run.get("qlat_file_gw_bucket_flux_col","qBucket")
        terrain_ro_col = run.get("qlat_file_terrain_runoff_col","qSfcLatRunoff")

        if geo_file_type=='NHDNetwork':
            # Parallel reading of qlateral data from CHRTOUT
            with Parallel(n_jobs=cpu_pool) as parallel:
                jobs = []
                for f in qlat_files:
                    jobs.append(
                        delayed(nhd_io.get_ql_from_chrtout)
                        #(f, qlat_file_value_col, gw_bucket_col, terrain_ro_col)
                        #delayed(nhd_io.get_ql_from_csv)
                        (f)                    
                    )
                ql_list = parallel(jobs)

            # get feature_id from a single CHRTOUT file
            with netCDF4.Dataset(qlat_files[0]) as ds:
                idx = ds.variables[qlat_file_index_col][:].filled()

            # package data into a DataFrame
            qlats_df = pd.DataFrame(
                np.stack(ql_list).T,
                index = idx,
                columns = range(len(qlat_files))
            )

            qlats_df = qlats_df[qlats_df.index.isin(segment_index)]
        elif geo_file_type=='HYFeaturesNetwork':
            dfs=[]
            for f in qlat_files:
                df = read_file(f).set_index(['feature_id']) 
                dfs.append(df)
            
            # lateral flows [m^3/s] are stored at NEXUS points with NEXUS ids
            nexuses_lateralflows_df = pd.concat(dfs, axis=1)  
            
            # Take flowpath ids entering NEXUS and replace NEXUS ids by the upstream flowpath ids
            qlats_df = pd.concat( (nexuses_lateralflows_df.loc[int(k)].rename(v)
                                for k,v in nexus_to_upstream_flowpath_dict.items() ), axis=1
                                ).T 
            qlats_df.columns=range(len(qlat_files))
            qlats_df = qlats_df[qlats_df.index.isin(segment_index)]

            # The segment_index has the full network set of segments/flowpaths. 
            # Whereas the set of flowpaths that are downstream of nexuses is a 
            # subset of the segment_index. Therefore, all of the segments/flowpaths
            # that are not accounted for in the set of flowpaths downstream of
            # nexuses need to be added to the qlateral dataframe and padded with
            # zeros.
            all_df = pd.DataFrame( np.zeros( (len(segment_index), len(qlats_df.columns)) ), index=segment_index,
                columns=qlats_df.columns )
            all_df.loc[ qlats_df.index ] = qlats_df
            qlats_df = all_df.sort_index()

    elif qlat_input_file:
        qlats_df = nhd_io.get_ql_from_csv(qlat_input_file)
    else:
        qlat_const = run.get("qlat_const", 0)
        qlats_df = pd.DataFrame(
            qlat_const,
            index=segment_index,
            columns=range(nts // qts_subdivisions),
            dtype="float32",
        )

    # TODO: Make a more sophisticated date-based filter
    max_col = 1 + nts // qts_subdivisions
    if len(qlats_df.columns) > max_col:
        qlats_df.drop(qlats_df.columns[max_col:], axis=1, inplace=True)

    if not segment_index.empty:
        qlats_df = qlats_df[qlats_df.index.isin(segment_index)]
    
    LOG.debug(
        "lateral inflow DataFrame creation complete in %s seconds." \
            % (time.time() - start_time)
            )

    return qlats_df

def nex_files_to_binary(nexus_files, binary_folder):
    for f in nexus_files:
        # read the csv file
        df = pd.read_csv(f, usecols=[1,2], names=['Datetime','qlat'])
        
        # convert and reformat datetime column
        df['Datetime']= pd.to_datetime(df['Datetime']).dt.strftime("%Y%m%d%H%M")

        # reformat the dataframe
        df['feature_id'] = get_id_from_filename(f)
        df = df.pivot(index="feature_id", columns="Datetime", values="qlat")
        df.columns.name = None

        for col in df.columns:
            table_new = pa.Table.from_pandas(df.loc[:, [col]])
            
            if not os.path.exists(f'{binary_folder}/{col}NEXOUT.parquet'):
                pq.write_table(table_new, f'{binary_folder}/{col}NEXOUT.parquet')
            
            else:
                table_old = pq.read_table(f'{binary_folder}/{col}NEXOUT.parquet')
                table = pa.concat_tables([table_old,table_new])
                pq.write_table(table, f'{binary_folder}/{col}NEXOUT.parquet')
    
    nexus_input_folder = binary_folder
    forcing_glob_filter = '*NEXOUT.parquet'

    return nexus_input_folder, forcing_glob_filter

def get_id_from_filename(file_name):
    id = os.path.splitext(file_name)[0].split('-')[1].split('_')[0]
    return int(id)

def read_file(file_name):
    extension = file_name.suffix
    if extension=='.csv':
        df = pd.read_csv(file_name)
    elif extension=='.parquet':
        df = pq.read_table(file_name).to_pandas().reset_index()
        df.index.name = None
    
    return df