import time
import pathlib
import logging
from datetime import datetime
from collections import defaultdict
from pathlib import Path
import os

import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd

import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io
from troute.nhd_network import reverse_dict
import troute.hyfeature_network_utilities as hnu

LOG = logging.getLogger('')

def read_geo_file(
    supernetwork_parameters,
    waterbody_parameters,
):
    
    geo_file_path = supernetwork_parameters["geo_file_path"]
        
    file_type = Path(geo_file_path).suffix
    if(  file_type == '.gpkg' ):        
        dataframe = read_geopkg(geo_file_path)
    elif( file_type == '.json') :
        edge_list = supernetwork_parameters['flowpath_edge_list']
        dataframe = read_json(geo_file_path, edge_list) 
    else:
        raise RuntimeError("Unsupported file type: {}".format(file_type))

    # Don't need the string prefix anymore, drop it
    mask = ~ dataframe['toid'].str.startswith("tnex") 
    dataframe = dataframe.apply(numeric_id, axis=1)
        
    # make the flowpath linkage, ignore the terminal nexus
    flowpath_dict = dict(zip(dataframe.loc[mask].toid, dataframe.loc[mask].id))
    
    # **********  need to be included in flowpath_attributes  *************
    dataframe['alt'] = 1.0 #FIXME get the right value for this... 

    cols = supernetwork_parameters.get('columns',None)
    
    if cols:
        dataframe = dataframe[list(cols.values())]
        # Rename parameter columns to standard names: from route-link names
        #        key: "link"
        #        downstream: "to"
        #        dx: "Length"
        #        n: "n"  # TODO: rename to `manningn`
        #        ncc: "nCC"  # TODO: rename to `mannningncc`
        #        s0: "So"  # TODO: rename to `bedslope`
        #        bw: "BtmWdth"  # TODO: rename to `bottomwidth`
        #        waterbody: "NHDWaterbodyComID"
        #        gages: "gages"
        #        tw: "TopWdth"  # TODO: rename to `topwidth`
        #        twcc: "TopWdthCC"  # TODO: rename to `topwidthcc`
        #        alt: "alt"
        #        musk: "MusK"
        #        musx: "MusX"
        #        cs: "ChSlp"  # TODO: rename to `sideslope`
        dataframe = dataframe.rename(columns=reverse_dict(cols))
        dataframe.set_index("key", inplace=True)
        dataframe.sort_index()

    # numeric code used to indicate network terminal segments
    terminal_code = supernetwork_parameters.get("terminal_code", 0)

    # There can be an externally determined terminal code -- that's this first value
    terminal_codes = set()
    terminal_codes.add(terminal_code)
    # ... but there may also be off-domain nodes that are not explicitly identified
    # but which are terminal (i.e., off-domain) as a result of a mask or some other
    # an interior domain truncation that results in a
    # otherwise valid node value being pointed to, but which is masked out or
    # being intentionally separated into another domain.
    terminal_codes = terminal_codes | set(
        dataframe[~dataframe["downstream"].isin(dataframe.index)]["downstream"].values
    )

    # build connections dictionary
    connections = nhd_network.extract_connections(
        dataframe, "downstream", terminal_codes=terminal_codes
    )

    #Load waterbody/reservoir info
    if waterbody_parameters:
        levelpool_params = waterbody_parameters.get('level_pool', None)
        if not levelpool_params:
            # FIXME should not be a hard requirement
            raise(RuntimeError("No supplied levelpool parameters in routing config"))
            
        lake_id = levelpool_params.get("level_pool_waterbody_id", "wb-id")
        waterbody_df = read_ngen_waterbody_df(
                    levelpool_params["level_pool_waterbody_parameter_file_path"],
                    lake_id,
                    )
            
        # Remove duplicate lake_ids and rows
        waterbody_df = (
                        waterbody_df.reset_index()
                        .drop_duplicates(subset=lake_id)
                        .set_index(lake_id)
                        )

        try:
            waterbody_types_df = read_ngen_waterbody_type_df(
                                    levelpool_params["reservoir_parameter_file"],
                                    lake_id,
                                    #self.waterbody_connections.values(),
                                    )
            # Remove duplicate lake_ids and rows
            waterbody_types_df =(
                                 waterbody_types_df.reset_index()
                                .drop_duplicates(subset=lake_id)
                                .set_index(lake_id)
                                )

        except ValueError:
            #FIXME any reservoir operations requires some type
            #So make this default to 1 (levelpool)
            waterbody_types_df = pd.DataFrame(index=waterbody_df.index)
            waterbody_types_df['reservoir_type'] = 1    
              
    return dataframe, flowpath_dict, connections, waterbody_df, waterbody_types_df, terminal_codes

def build_hyfeature_network(supernetwork_parameters,
                            waterbody_parameters,
):
    
    geo_file_path = supernetwork_parameters["geo_file_path"]
    cols          = supernetwork_parameters["columns"]
    terminal_code = supernetwork_parameters.get("terminal_code", 0)
    
    break_network_at_waterbodies = supernetwork_parameters.get("break_network_at_waterbodies", False)        
    break_network_at_gages       = supernetwork_parameters.get("break_network_at_gages", False)       
    break_points                 = {"break_network_at_waterbodies": break_network_at_waterbodies,
                                     "break_network_at_gages": break_network_at_gages}
        
    file_type = Path(geo_file_path).suffix
    if(  file_type == '.gpkg' ):        
        dataframe = hyf_network.read_geopkg(geo_file_path)
    elif( file_type == '.json') :
        edge_list = supernetwork_parameters['flowpath_edge_list']
        dataframe = hyf_network.read_json(geo_file_path, edge_list) 
    else:
        raise RuntimeError("Unsupported file type: {}".format(file_type))

    # Don't need the string prefix anymore, drop it
    mask = ~ dataframe['toid'].str.startswith("tnex") 
    dataframe = dataframe.apply(hyf_network.numeric_id, axis=1)
        
    # make the flowpath linkage, ignore the terminal nexus
    flowpath_dict = dict(zip(dataframe.loc[mask].toid, dataframe.loc[mask].id))
    waterbody_types_df = pd.DataFrame()
    waterbody_df = pd.DataFrame()
    waterbody_type_specified = False

 # **********  need to be included in flowpath_attributes  *************   
    # FIXME once again, order here can hurt....to hack `alt` in, either need to
    # put it as a column in the config, or do this AFTER the super constructor
    # otherwise the alt column gets sliced out...
    dataframe['alt'] = 1.0 #FIXME get the right value for this...
    
    #Load waterbody/reservoir info
    #For ngen HYFeatures, the reservoirs to be simulated
    #are determined by the lake.json file
    #we limit waterbody_connections to only the flowpaths
    #that coincide with a lake listed in this file
    #see `waterbody_connections`
    if waterbody_parameters:
        # FIXME later, DO ALL LAKE PARAMS BETTER
        levelpool_params = waterbody_parameters.get('level_pool', None)
        if not levelpool_params:
            # FIXME should not be a hard requirement
            raise(RuntimeError("No supplied levelpool parameters in routing config"))
            
        lake_id = levelpool_params.get("level_pool_waterbody_id", "wb-id")
        waterbody_df = read_ngen_waterbody_df(
                    levelpool_params["level_pool_waterbody_parameter_file_path"],
                    lake_id,
                    #self.waterbody_connections.values()
                    )
            
        # Remove duplicate lake_ids and rows
        waterbody_df = (
                        waterbody_df.reset_index()
                        .drop_duplicates(subset=lake_id)
                        .set_index(lake_id)
                        )
        waterbody_df["qd0"] = 0.0
        waterbody_df["h0"] = -1e9

        try:
            waterbody_types_df = read_ngen_waterbody_type_df(
                                    levelpool_params["reservoir_parameter_file"],
                                    lake_id,
                                    #self.waterbody_connections.values(),
                                    )
            # Remove duplicate lake_ids and rows
            waterbody_types_df =(
                                 waterbody_types_df.reset_index()
                                .drop_duplicates(subset=lake_id)
                                .set_index(lake_id)
                                )

        except ValueError:
            #FIXME any reservoir operations requires some type
            #So make this default to 1 (levelpool)
            waterbody_types_df = pd.DataFrame(index=waterbody_df.index)
            waterbody_types_df['reservoir_type'] = 1        
              
    return (dataframe,           
            flowpath_dict,  
            waterbody_types_df, 
            waterbody_df, 
            waterbody_type_specified,
            cols,
            terminal_code,
            break_points,
           )

def hyfeature_hybrid_routing_preprocess(
    connections,
    param_df,
    wbody_conn,
    gages,
    preprocessing_parameters,
    compute_parameters,
    waterbody_parameters, 
):
    '''
    Creation of routing network data objects. Logical ordering of lower-level
    function calls that build individual network data objects.
    
    Arguments
    ---------
    supernetwork_parameters      (dict): user input data re network extent
    waterbody_parameters         (dict): user input data re waterbodies
    preprocessing_parameters     (dict): user input data re preprocessing
    compute_parameters           (dict): user input data re compute configuration
    data_assimilation_parameters (dict): user input data re data assimilation
    
    Returns
    -------
    connections                 (dict of int: [int]): {segment id: [downsteram adjacent segment ids]}
    param_df                             (DataFrame): Hydraulic geometry and roughness parameters, by segment
    wbody_conn                    (dict of int: int): {segment id: associated lake id}
    waterbodies_df                       (DataFrame): Waterbody (reservoir) parameters
    waterbody_types_df                   (DataFrame): Waterbody type codes (1 - levelpool, 2 - USGS, 3 - USACE, 4 - RFC)
    break_network_at_waterbodies              (bool): If True, waterbodies occpy reaches of their own
    waterbody_type_specified                  (bool): If True, more than just levelpool waterbodies exist
    link_lake_crosswalk           (dict of int: int): {lake id: outlet segment id}
    independent_networks (dict of int: {int: [int]}): {tailwater id: {segment id: [upstream adjacent segment ids]}}
    reaches_bytw              (dict of int: [[int]]): {tailwater id: list or reach lists}
    rconn                       (dict of int: [int]): {segment id: [upstream adjacent segment ids]}
    pd.DataFrame.from_dict(gages)        (DataFrame): Gage ids and corresponding segment ids at which they are located
    diffusive_network_data            (dict or None): Network data objects for diffusive domain
    topobathy_df                         (DataFrame): Natural cross section data for diffusive domain
    
    Notes
    -----
    - waterbody_type_specified is likely an excessive return and can be removed and inferred from the 
      contents of waterbody_types_df
    - The values of the link_lake_crosswalk dictionary are the downstream-most segments within 
      the waterbody extent to which waterbody data are written. They are NOT the first segments 
      downsteram of the waterbody 
    '''

    #============================================================================
    # Establish diffusive domain for MC/diffusive hybrid simulations
    
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

    #============================================================================
    # Identify Independent Networks and Reaches by Network
    LOG.info("organizing connections into reaches ...")
    start_time = time.time() 
    gage_break_segments = set()
    wbody_break_segments = set()
    
    break_network_at_waterbodies = waterbody_parameters.get(
        "break_network_at_waterbodies", False
    )
    
    # if streamflow DA, then break network at gages
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
    # FIXME: Make this commented out alive
    '''
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
                 'link_lake_crosswalk': link_lake_crosswalk,
                 'independent_networks': independent_networks,
                 'reaches_bytw': reaches_bytw,
                 'rconn': rconn,
                 'link_gage_df': link_gage_df,
                 'usgs_lake_gage_crosswalk': usgs_lake_gage_crosswalk, 
                 'usace_lake_gage_crosswalk': usace_lake_gage_crosswalk,
                 'diffusive_network_data': diffusive_network_data,
                 'topobathy_data': topobathy_df,
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
        '''
    return(independent_networks,
           reaches_bytw,
           rconn,
           diffusive_network_data,
           topobathy_df, 
           refactored_diffusive_domain,
           refactored_reaches,
           unrefactored_topobathy_df,   
            )
    
def hyfeature_initial_warmstate_preprocess(
    # break_network_at_waterbodies,
    restart_parameters,
    # data_assimilation_parameters,
    segment_index,
    # waterbodies_df,
    # link_lake_crosswalk,
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
    - data_assimilation_parameters (dict): User-input data assimilation 
                                           parameters
    - segment_index        (Pandas Index): All segment IDs in the simulation 
                                           doamin
    - waterbodies_df   (Pandas DataFrame): Waterbody parameters
    - link_lake_crosswalk          (dict): Crosswalking between lake ids and the link
                                           id of the lake outlet segment
    
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

    #----------------------------------------------------------------------------
    # Assemble waterbody initial states (outflow and pool elevation
    #----------------------------------------------------------------------------
    '''
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
    '''

    #----------------------------------------------------------------------------
    # Assemble channel initial states (flow and depth)
    # also establish simulation initialization timestamp
    #----------------------------------------------------------------------------    
    start_time = time.time()
    LOG.info("setting channel initial states ...")

    # if lite restart file is provided, the read channel initial states from it
    if restart_parameters.get("lite_channel_restart_file", None):
        # FIXME: Change it for hyfeature!
        '''
        q0, t0 = nhd_io.read_lite_restart(
            restart_parameters['lite_channel_restart_file']
        )
        t0_str = None
        ''' 
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
        '''
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
        '''
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
            #waterbodies_df,
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

def hyfeature_forcing(
    run,
    forcing_parameters,
    hybrid_parameters,
    nexus_to_upstream_flowpath_dict,
    segment_index,
    cpu_pool,
    t0,
    coastal_boundary_depth_df,
):
    """
    Assemble model forcings. Forcings include hydrological lateral inflows (qlats)
    and coastal boundary depths for hybrid runs
    
    Aguments
    --------
    - run                (dict): List of forcing files pertaining to a 
                                 single run-set
    - forcing_parameters (dict): User-input simulation forcing parameters
    - hybrid_parameters  (dict): User-input simulation hybrid parameters
    - segment_index     (Int64): Reach segment ids
    - cpu_pool            (int): Number of CPUs in the process-parallel pool
    Returns
    -------
    - qlats_df                 (Pandas DataFrame): Lateral inflow data, indexed by 
                                                   segment ID
    - coastal_bounary_depth_df (Pandas DataFrame): Coastal boundary water depths,
                                                   indexed by segment ID
    
    Notes
    -----
    
    """
   
    # Unpack user-specified forcing parameters
    dt                           = forcing_parameters.get("dt", None)
    qts_subdivisions             = forcing_parameters.get("qts_subdivisions", None)
    nexus_input_folder           = forcing_parameters.get("nexus_input_folder", None)
    qlat_file_index_col          = forcing_parameters.get("qlat_file_index_col", "feature_id")
    qlat_file_value_col          = forcing_parameters.get("qlat_file_value_col", "q_lateral")
    qlat_file_gw_bucket_flux_col = forcing_parameters.get("qlat_file_gw_bucket_flux_col", "qBucket")
    qlat_file_terrain_runoff_col = forcing_parameters.get("qlat_file_terrain_runoff_col", "qSfcLatRunoff")

  
    # TODO: find a better way to deal with these defaults and overrides.
    run["t0"]                           = run.get("t0", t0)
    run["nts"]                          = run.get("nts")
    run["dt"]                           = run.get("dt", dt)
    run["qts_subdivisions"]             = run.get("qts_subdivisions", qts_subdivisions)
    run["nexus_input_folder"]           = run.get("nexus_input_folder", nexus_input_folder)
    run["qlat_file_index_col"]          = run.get("qlat_file_index_col", qlat_file_index_col)
    run["qlat_file_value_col"]          = run.get("qlat_file_value_col", qlat_file_value_col)
    run["qlat_file_gw_bucket_flux_col"] = run.get("qlat_file_gw_bucket_flux_col", qlat_file_gw_bucket_flux_col)
    run["qlat_file_terrain_runoff_col"] = run.get("qlat_file_terrain_runoff_col", qlat_file_terrain_runoff_col)
    
    #---------------------------------------------------------------------------
    # Assemble lateral inflow data
    #---------------------------------------------------------------------------
    
    start_time = time.time()
    LOG.info("Creating a DataFrame of lateral inflow forcings ...")

    # Place holder, if reading qlats from a file use this.
    # TODO: add an option for reading qlat data from BMI/model engine
    from_file = True
    if from_file:
        qlats_df = hnu.build_qlateral_array(
            run,
            cpu_pool,
            nexus_to_upstream_flowpath_dict,
            segment_index,
        )

    LOG.debug(
        "lateral inflow DataFrame creation complete in %s seconds." \
        % (time.time() - start_time)
    )

    #---------------------------------------------------------------------
    # Assemble coastal coupling data [WIP]
    #---------------------------------------------------------------------
    # Run if coastal_boundary_depth_df has not already been created:
    if coastal_boundary_depth_df.empty:
        coastal_boundary_elev_files = forcing_parameters.get('coastal_boundary_input_file', None) 
        coastal_boundary_domain_files = hybrid_parameters.get('coastal_boundary_domain', None)    
        
        if coastal_boundary_elev_files:
            start_time = time.time()    
            LOG.info("creating coastal dataframe ...")
            
            coastal_boundary_domain   = nhd_io.read_coastal_boundary_domain(coastal_boundary_domain_files)          
            coastal_boundary_depth_df = nhd_io.build_coastal_ncdf_dataframe(
                coastal_boundary_elev_files,
                coastal_boundary_domain,
            )
                
            LOG.debug(
                "coastal boundary elevation observation DataFrame creation complete in %s seconds." \
                % (time.time() - start_time)
            )            

    return qlats_df, coastal_boundary_depth_df

def read_ngen_waterbody_df(parm_file, lake_index_field="wb-id", lake_id_mask=None):
    """
    Reads .gpkg or lake.json file and prepares a dataframe, filtered
    to the relevant reservoirs, to provide the parameters
    for level-pool reservoir computation.
    """
    def node_key_func(x):
        return int(x[3:])
    if os.path.splitext(parm_file)[1]=='.gpkg':
        df = gpd.read_file(parm_file, layer="lake_attributes").set_index('id')
    elif os.path.splitext(parm_file)[1]=='.json':
        df = pd.read_json(parm_file, orient="index")

    df.index = df.index.map(node_key_func)
    df.index.name = lake_index_field

    if lake_id_mask:
        df = df.loc[lake_id_mask]
    return df

def read_ngen_waterbody_type_df(parm_file, lake_index_field="wb-id", lake_id_mask=None):
    """
    """
    #FIXME: this function is likely not correct. Unclear how we will get 
    # reservoir type from the gpkg files. Information should be in 'crosswalk'
    # layer, but as of now (Nov 22, 2022) there doesn't seem to be a differentiation
    # between USGS reservoirs, USACE reservoirs, or RFC reservoirs...
    def node_key_func(x):
        return int(x[3:])
    
    if os.path.splitext(parm_file)[1]=='.gpkg':
        df = gpd.read_file(parm_file, layer="crosswalk").set_index('id')
    elif os.path.splitext(parm_file)[1]=='.json':
        df = pd.read_json(parm_file, orient="index")

    df.index = df.index.map(node_key_func)
    df.index.name = lake_index_field
    if lake_id_mask:
        df = df.loc[lake_id_mask]
        
    return df

def read_geopkg(file_path):
    flowpaths = gpd.read_file(file_path, layer="flowpaths")
    attributes = gpd.read_file(file_path, layer="flowpath_attributes").drop('geometry', axis=1)
    #merge all relevant data into a single dataframe
    flowpaths = pd.merge(flowpaths, attributes, on='id')

    return flowpaths

def read_json(file_path, edge_list):
    dfs = []
    with open(edge_list) as edge_file:
        edge_data = json.load(edge_file)
        edge_map = {}
        for id_dict in edge_data:
            edge_map[ id_dict['id'] ] = id_dict['toid']
        with open(file_path) as data_file:
            json_data = json.load(data_file)  
            for key_wb, value_params in json_data.items():
                df = pd.json_normalize(value_params)
                df['id'] = key_wb
                df['toid'] = edge_map[key_wb]
                dfs.append(df)
        df_main = pd.concat(dfs, ignore_index=True)

    return df_main

def numeric_id(flowpath):
    id = flowpath['id'].split('-')[-1]
    toid = flowpath['toid'].split('-')[-1]
    flowpath['id'] = int(id)
    flowpath['toid'] = int(toid)

    return flowpath