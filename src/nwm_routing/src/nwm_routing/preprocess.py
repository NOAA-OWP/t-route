import time
import pandas as pd
import xarray as xr
from datetime import datetime
from collections import defaultdict
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io
import numpy as np
import sys

def nwm_network_preprocess(
    supernetwork_parameters,
    waterbody_parameters,
    showtiming=False,
    verbose=False,
    debuglevel=0,
):

    if verbose:
        print("creating supernetwork connections set")
    if showtiming:
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

    if verbose:
        print("supernetwork connections set complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

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

        # Check if hybrid-usgs, hybrid-usace, or rfc type reservoirs are set to true
        wbtype = "hybrid_and_rfc"
        wb_params_hybrid_and_rfc = waterbody_parameters.get(
            wbtype, defaultdict(list)
        )

        wbtype = "level_pool"
        wb_params_level_pool = waterbody_parameters.get(
            wbtype, defaultdict(list)
        )

        # Determine if any data assimilation reservoirs are activated, and if so, read
        # the reservoir parameter file
        if (
            wb_params_hybrid_and_rfc["reservoir_persistence_usgs"]
            or wb_params_hybrid_and_rfc["reservoir_persistence_usace"]
            or wb_params_hybrid_and_rfc["reservoir_rfc_forecasts"]
        ):

            waterbody_type_specified = True

            waterbody_types_df = nhd_io.read_reservoir_parameter_file(
                wb_params_hybrid_and_rfc["reservoir_parameter_file"],
                wb_params_level_pool["level_pool_waterbody_id"],
                wbody_conn.values(),
            )

            # Remove duplicate lake_ids and rows
            waterbody_types_df = (
                waterbody_types_df.reset_index()
                .drop_duplicates(subset="lake_id")
                .set_index("lake_id")
            )

    else:
        # Declare empty dataframes
        waterbody_types_df = pd.DataFrame()
        waterbodies_df = pd.DataFrame()

    # STEP 2: Identify Independent Networks and Reaches by Network
    if showtiming:
        start_time = time.time()
    if verbose:
        print("organizing connections into reaches ...")

    network_break_segments = set()
    if break_network_at_waterbodies:
        network_break_segments = network_break_segments.union(wbody_conn.values())
    if break_network_at_gages:
        network_break_segments = network_break_segments.union(gages.keys())

    independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections,
        network_break_segments,
    )
    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

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
    )


def nwm_initial_warmstate_preprocess(
    break_network_at_waterbodies,
    restart_parameters,
    data_assimilation_parameters,
    segment_index,
    waterbodies_df,
    segment_list=None,
    wbodies_list=None,
    showtiming=False,
    verbose=False,
    debuglevel=0,
):

    if break_network_at_waterbodies:
        ## STEP 3c: Handle Waterbody Initial States
        # TODO: move step 3c into function in nnu, like other functions wrapped in main()
        if showtiming:
            start_time = time.time()
        if verbose:
            print("setting waterbody initial states ...")

        if restart_parameters.get("wrf_hydro_waterbody_restart_file", None):
            waterbodies_initial_states_df = nhd_io.get_reservoir_restart_from_wrf_hydro(
                restart_parameters["wrf_hydro_waterbody_restart_file"],
                restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file"],
                restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file_field_name"],
                restart_parameters["wrf_hydro_waterbody_crosswalk_filter_file"],
                restart_parameters[
                    "wrf_hydro_waterbody_crosswalk_filter_file_field_name"
                ],
            )
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

        if verbose:
            print("waterbody initial states complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))
            start_time = time.time()

    # STEP 4: Handle Channel Initial States, set T0, and initialize LastObs
    if showtiming:
        start_time = time.time()
    if verbose:
        print("setting channel initial states ...")

    q0 = nnu.build_channel_initial_state(restart_parameters, segment_index)

    # STEP 4a: Set Channel States and T0
    if restart_parameters.get("wrf_hydro_channel_restart_file", None):
        channel_initial_states_file = restart_parameters[
            "wrf_hydro_channel_restart_file"
        ]
        t0_str = nhd_io.get_param_str(channel_initial_states_file, "Restart_Time")
    else:
        t0_str = "2015-08-16_00:00:00"

    t0 = datetime.strptime(t0_str, "%Y-%m-%d_%H:%M:%S")

    # STEP 4b: Set LastObs
    lastobs_df, da_parameter_dict = nnu.build_data_assimilation_lastobs(
        data_assimilation_parameters
    )

    if verbose:
        print("channel initial states complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))
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
    lastobs_index,
    warmstate_t0 = None,
    showtiming=False,
    verbose=False,
    debuglevel=0,
):

    # TODO: Harmonize the t0 parameter -- this
    # configuration permits the Warm-State derived
    # t0 to carry as the default if no other value is
    # provided -- we need to confirm this is
    # desireable behavior.
    t0 = forcing_parameters.get("t0", warmstate_t0)
    nts = forcing_parameters.get("nts", None)
    dt = forcing_parameters.get("dt", None)
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", None)
    qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
    qlat_file_index_col = forcing_parameters.get("qlat_file_index_col", None)
    qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", None)

    # TODO: find a better way to deal with these defaults and overrides.
    run["t0"] = run.get("t0", t0)
    run["nts"] = run.get("nts", nts)
    run["dt"] = run.get("dt", dt)
    run["qts_subdivisions"] = run.get("qts_subdivisions", qts_subdivisions)
    run["qlat_input_folder"] = run.get("qlat_input_folder", qlat_input_folder)
    run["qlat_file_index_col"] = run.get("qlat_file_index_col", qlat_file_index_col)
    run["qlat_file_value_col"] = run.get("qlat_file_value_col", qlat_file_value_col)

    if data_assimilation_parameters:
        data_assimilation_folder = data_assimilation_parameters.get("data_assimilation_timeslices_folder", None)
        data_assimilation_csv = data_assimilation_parameters.get("data_assimilation_csv", None)
        lastobs_file = data_assimilation_parameters.get("wrf_hydro_lastobs_file", None)
        lastobs_start = data_assimilation_parameters.get("wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0)
        lastobs_type = data_assimilation_parameters.get("wrf_lastobs_type", "error-based")
        lastobs_crosswalk_file = data_assimilation_parameters.get("wrf_hydro_da_channel_ID_crosswalk_file", None)
        da_decay_coefficient = data_assimilation_parameters.get("da_decay_coefficient", 120)

        da_run["data_assimilation_timeslices_folder"] = da_run.get("data_assimilation_timeslices_folder", data_assimilation_folder)
        da_run["data_assimilation_csv"] = da_run.get("data_assimilation_csv", data_assimilation_csv)
        da_run["wrf_hydro_lastobs_file"] = da_run.get("wrf_hydro_lastobs_file", lastobs_file)
        da_run["wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0] = da_run.get("wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", lastobs_start)
        da_run["wrf_lastobs_type"] = da_run.get("wrf_lastobs_type", lastobs_type)
        da_run["wrf_hydro_da_channel_ID_crosswalk_file"] = da_run.get("wrf_hydro_da_channel_ID_crosswalk_file", lastobs_crosswalk_file)
        da_run["da_decay_coefficient"] = da_run.get("da_decay_coefficient", da_decay_coefficient)

    # STEP 5: Read (or set) QLateral Inputs
    if showtiming:
        start_time = time.time()
    if verbose:
        print("creating qlateral array ...")

    qlats_df = nnu.build_qlateral_array(
        run,
        segment_index,
    )

    if verbose:
        print("qlateral array complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 6
    data_assimilation_csv = da_run.get("data_assimilation_csv", None)
    data_assimilation_folder = da_run.get("data_assimilation_timeslices_folder", None)
    if data_assimilation_csv or data_assimilation_folder:

        if data_assimilation_folder and data_assimilation_csv:
            print(
                "Please select data_assimilation_parameters_folder + data_assimilation_filter or data_assimilation_csv not both."
            )

        if showtiming:
            start_time = time.time()
        if verbose:
            print("creating usgs time_slice data array ...")

        usgs_df = nnu.build_data_assimilation_usgs_df(da_run, run, lastobs_index)

        if verbose:
            print("usgs array complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))

    else:
        usgs_df = pd.DataFrame()

    # STEP 7
    coastal_boundary_elev = forcing_parameters.get("coastal_boundary_elev_data", None)
    coastal_ncdf = forcing_parameters.get("coastal_ncdf", None)

    if coastal_boundary_elev:
        print("creating coastal dataframe ...")
        coastal_df = nhd_io.build_coastal_dataframe(coastal_boundary_elev)

    if coastal_ncdf:
        print("creating coastal ncdf dataframe ...")
        coastal_ncdf_df = nhd_io.build_coastal_ncdf_dataframe(coastal_ncdf)

    # TODO: disentangle the implicit (run) and explicit (qlats_df, usgs_df) returns
    return qlats_df, usgs_df
    parallel_compute_method = compute_parameters.get("parallel_compute_method", None)
    subnetwork_target_size = compute_parameters.get("subnetwork_target_size", 1)
    cpu_pool = compute_parameters.get("cpu_pool", None)
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", 1)
    compute_kernel = compute_parameters.get("compute_kernel", "V02-caching")
    assume_short_ts = compute_parameters.get("assume_short_ts", False)
    return_courant = compute_parameters.get("return_courant", False)

def nwm_save_preprocessing(save_preprocessing,connections,rconn,wbody_conn,reaches_bytw,verbose,showtiming,debuglevel,independent_networks,param_df,waterbody_parameters,
    waterbody_types_df,waterbody_type_specified,diffusive_parameters,break_network_at_waterbodies,link_gage_df=None):
    preprocessed_outputs = {}
    preprocessed_outputs.update({'connections': connections,'rconn': rconn,'wbody_conn':wbody_conn,'reaches_bytw':reaches_bytw, 
    'independent_networks':independent_networks,'param_df':param_df,
    'waterbody_parameters':waterbody_parameters,'waterbody_types_df':waterbody_types_df,
    'waterbody_type_specified':waterbody_type_specified,'diffusive_parameters':diffusive_parameters,'showtiming':showtiming,'verbose':verbose,'debuglevel':debuglevel,
    'break_network_at_waterbodies':break_network_at_waterbodies,'link_gage_df':link_gage_df})
    np.save(save_preprocessing+'/preprocessed_outputs.npy', preprocessed_outputs)
    print("Saved preprocessed components to file.")
    sys.exit()
    #'q0':q0,'waterbodies_df':waterbodies_df,'lastobs_df':lastobs_df,'da_parameter_dict':da_parameter_dict,'run_sets':run_sets,'da_sets':da_sets,'parity_sets':parity_sets,
def nwm_load_preprocessing(load_preprocessing):
        preprocessed_outputs = np.load(load_preprocessing+'/preprocessed_outputs.npy',allow_pickle='TRUE').item()
        # run_sets = preprocessed_outputs.get('run_sets', None)
        connections = preprocessed_outputs.get('connections', None)
        rconn = preprocessed_outputs.get('rconn', None)
        wbody_conn = preprocessed_outputs.get('wbody_conn', None)
        reaches_bytw = preprocessed_outputs.get('reaches_bytw', None)
        # parallel_compute_method = preprocessed_outputs.get('parallel_compute_method', None)
        # compute_kernel = preprocessed_outputs.get('compute_kernel', None)
        # subnetwork_target_size = preprocessed_outputs.get('subnetwork_target_size', None)
        # cpu_pool = preprocessed_outputs.get('cpu_pool', None)
        # qts_subdivisions = preprocessed_outputs.get('qts_subdivisions', None)
        independent_networks = preprocessed_outputs.get('independent_networks', None)
        param_df = preprocessed_outputs.get('param_df', None)
        # q0 = preprocessed_outputs.get('q0', None)
        # qlats = preprocessed_outputs.get('qlats', None)
        # usgs_df = preprocessed_outputs.get('usgs_df', None)
        # lastobs_df = preprocessed_outputs.get('lastobs_df', None)
        # da_parameter_dict = preprocessed_outputs.get('da_parameter_dict', None)
        # assume_short_ts = preprocessed_outputs.get('assume_short_ts', None)
        # return_courant = preprocessed_outputs.get('return_courant', None)
        waterbodies_df = preprocessed_outputs.get('waterbodies_df', None)
        waterbody_parameters = preprocessed_outputs.get('waterbody_parameters', None)
        waterbody_type_specified = preprocessed_outputs.get('waterbody_type_specified', None)
        diffusive_parameters = preprocessed_outputs.get('diffusive_parameters', None)
        showtiming = preprocessed_outputs.get('showtiming', None)
        verbose = preprocessed_outputs.get('verbose', None)
        debuglevel = preprocessed_outputs.get('debuglevel', None)
        # parity_sets = preprocessed_outputs.get('parity_sets', None)
        waterbody_types_df = preprocessed_outputs.get('waterbody_types_df',None)
        # da_sets = preprocessed_outputs.get('da_sets',None)
        break_network_at_waterbodies = preprocessed_outputs.get('break_network_at_waterbodies',None)
        link_gage_df = preprocessed_outputs.get('link_gage_df',None)
        print("Loaded preprocessed components to file.")
        if showtiming:
            main_start_time = time.time()
        return (connections,rconn, wbody_conn,reaches_bytw,independent_networks,param_df,waterbodies_df,waterbody_parameters,waterbody_type_specified,diffusive_parameters
        ,showtiming,verbose,debuglevel,waterbody_types_df,break_network_at_waterbodies,link_gage_df)