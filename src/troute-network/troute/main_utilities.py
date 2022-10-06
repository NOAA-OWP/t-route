import argparse
import time
import logging
import pandas as pd
import numpy as np
from datetime import timedelta
#from .input import _input_handler_v03
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_io as nhd_io

#from .output import nwm_output_generator
from troute.NHDNetwork import NHDNetwork
from troute.HYFeaturesNetwork import HYFeaturesNetwork
from troute.DataAssimilation import AllDA
from troute.routing.compute import compute_nhd_routing_v02, compute_diffusive_routing

LOG = logging.getLogger('')

def initialize_network(argv):
    
    args = _handle_args_v03(argv)
    
    # unpack user inputs
    (
        log_parameters,
        preprocessing_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        hybrid_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    ) = _input_handler_v03(args)

    run_parameters = {
        'dt': forcing_parameters.get('dt'),
        'nts': forcing_parameters.get('nts'),
        'cpu_pool': compute_parameters.get('cpu_pool')
    }
    
#    showtiming = log_parameters.get("showtiming", None)
#    if showtiming:
#        task_times = {}
#        task_times['forcing_time'] = 0
#        task_times['route_time'] = 0
#        task_times['output_time'] = 0
#        main_start_time = time.time()
    
    cpu_pool = compute_parameters.get("cpu_pool", None)
 
    # Build routing network data objects. Network data objects specify river 
    # network connectivity, channel geometry, and waterbody parameters. Also
    # perform initial warmstate preprocess.
#    if showtiming:
#        network_start_time = time.time()

    if "ngen_nexus_file" in supernetwork_parameters:
         network = HYFeaturesNetwork(supernetwork_parameters,
                                    waterbody_parameters=waterbody_parameters,
                                    restart_parameters=restart_parameters,
                                    forcing_parameters=forcing_parameters,
                                    verbose=verbose, showtiming=showtiming) 
    else:
        network = NHDNetwork(supernetwork_parameters=supernetwork_parameters,
                             waterbody_parameters=waterbody_parameters,
                             restart_parameters=restart_parameters,
                             forcing_parameters=forcing_parameters,
                             compute_parameters=compute_parameters,
                             data_assimilation_parameters=data_assimilation_parameters,
                             preprocessing_parameters=preprocessing_parameters,
                             verbose=True,
                             showtiming=False #showtiming,          
                            )
    
#    if showtiming:
#        network_end_time = time.time()
#        task_times['network_time'] = network_end_time - network_start_time
    
    all_parameters = [log_parameters, 
                      preprocessing_parameters,
                      supernetwork_parameters,
                      waterbody_parameters,
                      compute_parameters,
                      forcing_parameters,
                      restart_parameters,
                      hybrid_parameters,
                      output_parameters,
                      parity_parameters,
                      data_assimilation_parameters,
                      run_parameters
                     ]
    
    return network, log_parameters, preprocessing_parameters, supernetwork_parameters, waterbody_parameters, compute_parameters, forcing_parameters, restart_parameters, hybrid_parameters, output_parameters, parity_parameters, data_assimilation_parameters, run_parameters

def build_run_sets(network, forcing_parameters, compute_parameters, data_assimilation_parameters,
                   output_parameters, parity_parameters):
    
    # Create run_sets: sets of forcing files for each loop
    run_sets = nnu.build_forcing_sets(forcing_parameters, network.t0)

    # Create da_sets: sets of TimeSlice files for each loop
    if "data_assimilation_parameters" in compute_parameters:
        da_sets = nnu.build_da_sets(data_assimilation_parameters, run_sets, network.t0)
        
    # Create parity_sets: sets of CHRTOUT files against which to compare t-route flows
    if "wrf_hydro_parity_check" in output_parameters:
        parity_sets = nnu.build_parity_sets(parity_parameters, run_sets)
    else:
        parity_sets = []
    
    return run_sets, da_sets, parity_sets

def build_forcings(network, run, forcing_parameters, hybrid_parameters, compute_parameters):
    
    cpu_pool = compute_parameters.get('cpu_pool', None)
    # Create forcing data within network object for first loop iteration
    network.assemble_forcings(run, forcing_parameters, hybrid_parameters, cpu_pool)
    
    return network

def build_data_assimilation(network, data_assimilation_parameters, waterbody_parameters, da_run, forcing_parameters, compute_parameters):
    
    #FIXME: hack to get run_parameters. This is done in input_handler_v2. Probably need
    # to find a better way to do this here though...
    if not 'run_parameters' in locals():
        run_parameters = {'dt': forcing_parameters.get('dt'),
                          'nts': forcing_parameters.get('nts'),
                          'cpu_pool': compute_parameters.get('cpu_pool', None)}
    
    # Create data assimilation object from da_sets for first loop iteration
    data_assimilation = AllDA(data_assimilation_parameters,
                              run_parameters,
                              waterbody_parameters,
                              network,
                              da_run)
    
#    if showtiming:
#        forcing_end_time = time.time()
#        task_times['forcing_time'] += forcing_end_time - network_end_time
    
    return data_assimilation

def run_routing(network, data_assimilation, run_sets, da_sets, compute_parameters, forcing_parameters, waterbody_parameters,
                output_parameters, hybrid_parameters, data_assimilation_parameters, run_parameters, parity_sets):
    '''
    
    '''
    parallel_compute_method = compute_parameters.get("parallel_compute_method", None)
    subnetwork_target_size = compute_parameters.get("subnetwork_target_size", 1)
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", 1)
    compute_kernel = compute_parameters.get("compute_kernel", "V02-caching")
    assume_short_ts = compute_parameters.get("assume_short_ts", False)
    return_courant = compute_parameters.get("return_courant", False)
    cpu_pool = compute_parameters.get("cpu_pool", 1)
        
    # Pass empty subnetwork list to nwm_route. These objects will be calculated/populated
    # on first iteration of for loop only. For additional loops this will be passed
    # to function from inital loop.     
    subnetwork_list = [None, None, None]
    
    for run_set_iterator, run in enumerate(run_sets):
        
        t0 = run.get("t0")
        dt = run.get("dt")
        nts = run.get("nts")

        if parity_sets:
            parity_sets[run_set_iterator]["dt"] = dt
            parity_sets[run_set_iterator]["nts"] = nts

#        if showtiming:
#            route_start_time = time.time()

        run_results = nwm_route(
            network.connections, 
            network.reverse_network, 
            network.waterbody_connections, 
            network._reaches_by_tw, ## check: def name is different from return self._ ..
            parallel_compute_method,
            compute_kernel,
            subnetwork_target_size,
            cpu_pool,
            network.t0,  ## check if t0 is being updated
            dt,
            nts,
            qts_subdivisions,
            network.independent_networks, 
            network.dataframe,
            network.q0,
            network._qlateral,
            data_assimilation.usgs_df,
            data_assimilation.lastobs_df,
            data_assimilation.reservoir_usgs_df,
            data_assimilation.reservoir_usgs_param_df,
            data_assimilation.reservoir_usace_df,
            data_assimilation.reservoir_usace_param_df,
            data_assimilation.assimilation_parameters,
            assume_short_ts,
            return_courant,
            network._waterbody_df, ## check:  network._waterbody_df ?? def name is different from return self._ ..
            waterbody_parameters,
            network._waterbody_types_df, ## check:  network._waterbody_types_df ?? def name is different from return self._ ..
            network.waterbody_type_specified,
            network.diffusive_network_data,
            network.topobathy_df,
            network.refactored_diffusive_domain,
            network.refactored_reaches,
            subnetwork_list,
            network.coastal_boundary_depth_df,
            network.unrefactored_topobathy_df,
        )
        
        # returns list, first item is run result, second item is subnetwork items
        subnetwork_list = run_results[1]
        run_results = run_results[0]
        
#        if showtiming:
#            route_end_time = time.time()
#            task_times['route_time'] += route_end_time - route_start_time

        # create initial conditions for next loop itteration
        network.new_nhd_q0(run_results)
        network.update_waterbody_water_elevation()               
        
        if run_set_iterator < len(run_sets) - 1:
            # update t0
            network.new_t0(dt,nts)
            
            # update forcing data
            network.assemble_forcings(run_sets[run_set_iterator + 1],
                                      forcing_parameters,
                                      hybrid_parameters,
                                      cpu_pool)
            
            # get reservoir DA initial parameters for next loop iteration
            data_assimilation.update(run_results,
                                     data_assimilation_parameters,
                                     run_parameters,
                                     network,
                                     da_sets[run_set_iterator + 1])
            
#            if showtiming:
#                forcing_end_time = time.time()
#                task_times['forcing_time'] += forcing_end_time - route_end_time
        
        # TODO move the conditional call to write_lite_restart to nwm_output_generator.
#        if showtiming:
#            output_start_time = time.time()
            
        if "lite_restart" in output_parameters:
            nhd_io.write_lite_restart(
                network.q0, 
                network._waterbody_df, 
                t0 + timedelta(seconds = dt * nts), 
                output_parameters['lite_restart']
            )
        '''
        nwm_output_generator(
            run,
            run_results,
            supernetwork_parameters,
            output_parameters,
            parity_parameters,
            restart_parameters,
            parity_sets[run_set_iterator] if parity_parameters else {},
            qts_subdivisions,
            compute_parameters.get("return_courant", False),
            cpu_pool,
            network._waterbody_df,
            network._waterbody_types_df,
            data_assimilation_parameters,
            data_assimilation.lastobs_df,
            network.link_gage_df,
            network.link_lake_crosswalk,
        )
        '''
#        if showtiming:
#            output_end_time = time.time()
#            task_times['output_time'] += output_end_time - output_start_time
    return run_results

        
def _handle_args_v03(argv):
    '''
    Handle command line input argument - filepath of configuration file
    '''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--custom-input-file",
        dest="custom_input_file",
        help="Path of a .yaml or .json file containing model configuration parameters. See doc/v3_doc.yaml",
    )
    
    return parser.parse_args(argv)

def _input_handler_v03(args):
    
    '''
    Read user inputs from configuration file and set logging level
    
    Arguments
    ---------
    Args (argparse.Namespace): Command line input arguments
    
    Returns
    -------
    log_parameters               (dict): Input parameters re logging
    preprocessing_parameters     (dict): Input parameters re preprocessing
    supernetwork_parameters      (dict): Input parameters re network extent
    waterbody_parameters         (dict): Input parameters re waterbodies
    compute_parameters           (dict): Input parameters re computation settings
    forcing_parameters           (dict): Input parameters re model forcings
    restart_parameters           (dict): Input parameters re model restart
    hybrid_parameters            (dict): Input parameters re MC/diffusive wave model
    output_parameters            (dict): Input parameters re output writing
    parity_parameters            (dict): Input parameters re parity assessment
    data_assimilation_parameters (dict): Input parameters re data assimilation
    '''
    # get name of user configuration file (e.g. test.yaml)
    custom_input_file = args.custom_input_file

    # read data from user configuration file
    (
        log_parameters,
        preprocessing_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        hybrid_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    ) = nhd_io.read_config_file(custom_input_file)

    '''
    # configure python logger
    log_level_set(log_parameters)
    LOG = logging.getLogger('')

    # if log level is at or below DEBUG, then check user inputs
    if LOG.level <= 10: # DEBUG
        check_inputs(
                log_parameters,
                preprocessing_parameters,
                supernetwork_parameters,
                waterbody_parameters,
                compute_parameters,
                forcing_parameters,
                restart_parameters,
                hybrid_parameters,
                output_parameters,
                parity_parameters,
                data_assimilation_parameters
                )
    ''' 
    return (
        log_parameters,
        preprocessing_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        hybrid_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    )

def nwm_route(
    downstream_connections,
    upstream_connections,
    waterbodies_in_connections,
    reaches_bytw,
    parallel_compute_method,
    compute_kernel,
    subnetwork_target_size,
    cpu_pool,
    t0,
    dt,
    nts,
    qts_subdivisions,
    independent_networks,
    param_df,
    q0,
    qlats,
    usgs_df,
    lastobs_df,
    reservoir_usgs_df,
    reservoir_usgs_param_df,
    reservoir_usace_df,
    reservoir_usace_param_df,
    da_parameter_dict,
    assume_short_ts,
    return_courant,
    waterbodies_df,
    waterbody_parameters,
    waterbody_types_df,
    waterbody_type_specified,
    diffusive_network_data,
    topobathy,
    refactored_diffusive_domain,
    refactored_reaches,
    subnetwork_list,
    coastal_boundary_depth_df,
    nonrefactored_topobathy,
):

    ################### Main Execution Loop across ordered networks
    
    
    start_time = time.time()

    if return_courant:
        LOG.info(
            f"executing routing computation, with Courant evaluation metrics returned"
        )
    else:
        LOG.info(f"executing routing computation ...")

    start_time_mc = time.time()
    results = compute_nhd_routing_v02(
        downstream_connections,
        upstream_connections,
        waterbodies_in_connections,
        reaches_bytw,
        compute_kernel,
        parallel_compute_method,
        subnetwork_target_size,  # The default here might be the whole network or some percentage...
        cpu_pool,
        t0,
        dt,
        nts,
        qts_subdivisions,
        independent_networks,
        param_df,
        q0,
        qlats,
        usgs_df,
        lastobs_df,
        reservoir_usgs_df,
        reservoir_usgs_param_df,
        reservoir_usace_df,
        reservoir_usace_param_df,
        da_parameter_dict,
        assume_short_ts,
        return_courant,
        waterbodies_df,
        waterbody_parameters,
        waterbody_types_df,
        waterbody_type_specified,
        subnetwork_list,
    )
    LOG.debug("MC computation complete in %s seconds." % (time.time() - start_time_mc))
    # returns list, first item is run result, second item is subnetwork items
    subnetwork_list = results[1]
    results = results[0]
    
    # run diffusive side of a hybrid simulation
    if diffusive_network_data:
        start_time_diff = time.time()
        '''
        # retrieve MC-computed streamflow value at upstream boundary of diffusive mainstem
        qvd_columns = pd.MultiIndex.from_product(
            [range(nts), ["q", "v", "d"]]
        ).to_flat_index()
        flowveldepth = pd.concat(
            [pd.DataFrame(r[1], index=r[0], columns=qvd_columns) for r in results],
            copy=False,
        )
        '''
        #upstream_boundary_flow={}
        #for tw,v in  diffusive_network_data.items():
        #    upstream_boundary_link     = diffusive_network_data[tw]['upstream_boundary_link']
        #    flow_              = flowveldepth.loc[upstream_boundary_link][0::3]
            # the very first value at time (0,q) is flow value at the first time step after initial time.
        #    upstream_boundary_flow[tw] = flow_         
          

        # call diffusive wave simulation and append results to MC results
        results.extend(
            compute_diffusive_routing(
                results,
                diffusive_network_data,
                cpu_pool,
                t0,
                dt,
                nts,
                q0,
                qlats,
                qts_subdivisions,
                usgs_df,
                lastobs_df,
                da_parameter_dict,
                waterbodies_df,
                topobathy,
                refactored_diffusive_domain,
                refactored_reaches,
                coastal_boundary_depth_df,
                nonrefactored_topobathy, 
            )
        )
        LOG.debug("Diffusive computation complete in %s seconds." % (time.time() - start_time_diff))

    LOG.debug("ordered reach computation complete in %s seconds." % (time.time() - start_time))

    return results, subnetwork_list

def create_output_dataframes(results, run_sets, waterbodies_df, link_lake_crosswalk):
    
    nts = run_sets[len(run_sets) - 1].get('nts')
    dt = run_sets[len(run_sets) - 1].get('dt')
    qts_subdivisions = run_sets[len(run_sets) - 1].get('qts_subdivisions')
    
    qvd_columns = pd.MultiIndex.from_product(
        [range(nts), ["q", "v", "d"]]
    ).to_flat_index()
    
    flowveldepth = pd.concat(
        [pd.DataFrame(r[1], index=r[0], columns=qvd_columns) for r in results], copy=False,
    )
    
    # create waterbody dataframe for output to netcdf file
    i_columns = pd.MultiIndex.from_product(
        [range(nts), ["i"]]
    ).to_flat_index()
    
    wbdy = pd.concat(
        [pd.DataFrame(r[6], index=r[0], columns=i_columns) for r in results],
        copy=False,
    )
    
    wbdy_id_list = waterbodies_df.index.values.tolist()
    '''
    flow_df = flowveldepth.loc[wbdy_id_list]
    wbdy = wbdy.loc[wbdy_id_list]
    
    timestep, variable = zip(*flow_df.columns.tolist())
    timestep_index = np.where(((np.array(list(set(list(timestep)))) + 1) * dt) % (dt * qts_subdivisions) == 0)
    ts_set = set(timestep_index[0].tolist())
    flow_df_col_index = [i for i, e in enumerate(timestep) if e in ts_set]
    flow_df = flow_df.iloc[:,flow_df_col_index]
    
    timestep, variable = zip(*wbdy.columns.tolist())
    timestep_index = np.where(((np.array(list(set(list(timestep)))) + 1) * dt) % (dt * qts_subdivisions) == 0)
    ts_set = set(timestep_index[0].tolist())
    wbdy_col_index = [i for i, e in enumerate(timestep) if e in ts_set]
    
    i_lakeout_df = wbdy.iloc[:,wbdy_col_index]
    q_lakeout_df = flow_df.iloc[:,0::3]
    d_lakeout_df = flow_df.iloc[:,2::3]
    '''
    i_lakeout_df = wbdy.loc[wbdy_id_list]
    q_lakeout_df = flowveldepth.loc[wbdy_id_list].iloc[:,0::3]
    d_lakeout_df = flowveldepth.loc[wbdy_id_list].iloc[:,2::3]
    # lakeout = pd.concat([i_df, q_df, d_df], axis=1)
    
    # replace waterbody lake_ids with outlet link ids
    flowveldepth = _reindex_lake_to_link_id(flowveldepth, link_lake_crosswalk)
    
    q_channel_df = flowveldepth.iloc[:,0::3]
    v_channel_df = flowveldepth.iloc[:,1::3]
    d_channel_df = flowveldepth.iloc[:,2::3]

    date_column_names = pd.date_range(start=run_sets[-1].get('t0') + timedelta(seconds=run_sets[-1].get('dt')), 
              end=run_sets[-1].get('final_timestamp'), 
              freq=str(run_sets[-1].get('dt'))+'S')
    
    i_lakeout_df.columns = date_column_names
    q_lakeout_df.columns = date_column_names
    d_lakeout_df.columns = date_column_names
    q_channel_df.columns = date_column_names
    v_channel_df.columns = date_column_names
    d_channel_df.columns = date_column_names
    
    return q_channel_df, v_channel_df, d_channel_df, i_lakeout_df, q_lakeout_df, d_lakeout_df
    
def _reindex_lake_to_link_id(target_df, crosswalk):
    '''
    Utility function for replacing lake ID index values
    with link ID values in a dataframe. This is used to 
    reinedex results dataframes
    
    Arguments:
    ----------
    - target_df (DataFrame): Data frame to be reinexed
    - crosswalk      (dict): Relates lake ids to outlet link ids
    
    Returns:
    --------
    - target_df (DataFrame): Re-indexed with link ids replacing 
                             lake ids
    '''

    # evaluate intersection of lake ids and target_df index values
    # i.e. what are the index positions of lake ids that need replacing?
    lakeids = np.fromiter(crosswalk.keys(), dtype = int)
    idxs = target_df.index.to_numpy()
    lake_index_intersect = np.intersect1d(
        idxs, 
        lakeids, 
        return_indices = True
    )

    # replace lake ids with link IDs in the target_df index array
    linkids = np.fromiter(crosswalk.values(), dtype = int)
    idxs[lake_index_intersect[1]] = linkids[lake_index_intersect[2]]

    # (re) set the target_df index
    target_df.set_index(idxs, inplace = True)
    
    return target_df