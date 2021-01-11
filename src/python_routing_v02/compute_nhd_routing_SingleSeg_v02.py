#!/usr/bin/env python
# coding: utf-8
# example usage: python compute_nhd_routing_SingleSeg.py -v -t -w -n Mainstems_CONUS
# python compute_nhd_routing_SingleSeg_v02.py --test -t -v --debuglevel 1
# python compute_nhd_routing_SingleSeg_v02.py --test-full-pocono -t -v --debuglevel 1


# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
## Parallel execution
import os
import sys
import time
import numpy as np
import argparse
import pathlib
import glob
import pandas as pd
from functools import partial
from joblib import delayed, Parallel
from itertools import chain, islice
from operator import itemgetter


def _handle_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--debuglevel",
        help="Set the debuglevel",
        dest="debuglevel",
        choices=[0, 1, 2, 3],
        default=0,
        type=int,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose output (leave blank for quiet output)",
        dest="verbose",
        action="store_true",
    )
    parser.add_argument(
        "--qlat-dt",
        "--qlateral-time-step",
        help="Set the default qlateral timestep length",
        dest="qdt",
        default=3600,
    )
    parser.add_argument(
        "--qN",
        "--qts-subdivisions",
        help="number of simulation timesteps per qlateral timestep",
        dest="qts_subdivisions",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--dt",
        "--simulation-time-step",
        help="Set the default simulation timestep length",
        dest="dt",
        default=300,
    )
    parser.add_argument(
        "--nts",
        "--number-of-simulation-timesteps",
        help="Set the number of timesteps to execute. If used with ql_file or ql_folder, nts must be less than len(ql) x qN.",
        dest="nts",
        default=144,
        type=int,
    )
    parser.add_argument(
        "--test",
        "--run-pocono2-test-example",
        help="Use the data values stored in the repository for a test of the Pocono network",
        dest="run_pocono2_test",
        action="store_true",
    )
    parser.add_argument(
        "--test-full-pocono",
        "--run-pocono1-test-example",
        help="Use the data values stored in the repository for a test of the Mainstems_CONUS network",
        dest="run_pocono1_test",
        action="store_true",
    )
    parser.add_argument(
        "--sts",
        "--assume-short-ts",
        help="Use the previous timestep value for upstream flow",
        dest="assume_short_ts",
        action="store_true",
    )
    parser.add_argument(
        "-ocsv",
        "--write-output-csv",
        nargs="?",
        help="Write csv output files to this folder (omit flag for no csv writing)",
        dest="csv_output_folder",
        const="../../test/output/text",
    )
    parser.add_argument(
        "-t",
        "--showtiming",
        help="Set the showtiming (leave blank for no timing information)",
        dest="showtiming",
        action="store_true",
    )
    parser.add_argument(
        "-w",
        "--break-at-waterbodies",
        help="Use the waterbodies in the route-link dataset to divide the computation (leave blank for no splitting)",
        dest="break_network_at_waterbodies",
        action="store_true",
    )
    parser.add_argument(
        "--parallel",
        nargs="?",
        help="Use the parallel computation engine (omit flag for serial computation)",
        dest="parallel_compute_method",
        const="by-network",
    )
    parser.add_argument(
        "--cpu-pool",
        help="Assign the number of cores to multiprocess across.",
        dest="cpu_pool",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--compute-method",
        help="Use the cython version of the compute_network code (enter additional flag options for other compute_network possibilities).",
        dest="compute_method",
        default="standard cython compute network",
    )
    supernetwork_arg_group = parser.add_mutually_exclusive_group()
    supernetwork_arg_group.add_argument(
        "-n",
        "--supernetwork",
        help="Choose from among the pre-programmed supernetworks (Pocono_TEST1, Pocono_TEST2, LowerColorado_Conchos_FULL_RES, Brazos_LowerColorado_ge5, Brazos_LowerColorado_FULL_RES, Brazos_LowerColorado_Named_Streams, CONUS_ge5, Mainstems_CONUS, CONUS_Named_Streams, CONUS_FULL_RES_v20)",
        choices=[
            "Pocono_TEST1",
            "Pocono_TEST2",
            "LowerColorado_Conchos_FULL_RES",
            "Brazos_LowerColorado_ge5",
            "Brazos_LowerColorado_FULL_RES",
            "Brazos_LowerColorado_Named_Streams",
            "CONUS_ge5",
            "Mainstems_CONUS",
            "CONUS_Named_Streams",
            "CONUS_FULL_RES_v20",
            "CapeFear_FULL_RES",
            "Florence_FULL_RES",
        ],
        # TODO: accept multiple or a Path (argparse Action perhaps)
        # action='append',
        # nargs=1,
        dest="supernetwork",
        default="Pocono_TEST1",
    )
    supernetwork_arg_group.add_argument(
        "-f",
        "--custom-input-file",
        dest="custom_input_file",
        help="OR... please enter the path of a .yaml or .json file containing a custom supernetwork information. See for example test/input/yaml/CustomInput.yaml and test/input/json/CustomInput.json.",
    )
    parser.add_argument(
        "--wrf-hydro-channel-restart-file",
        dest="wrf_hydro_channel_restart_file",
        help="provide a WRF-Hydro channel warm state file (may be the same as waterbody restart file)",
    )
    parser.add_argument(
        "--wrf-hydro-channel-ID-crosswalk-file",
        dest="wrf_hydro_channel_ID_crosswalk_file",
        help="provide an xarray-readable file that defines the order of the outputs in the channel restart file. Specify the ID field with --wrf_hydro_channel_ID_crosswalk_file_field_name",
    )
    parser.add_argument(
        "--wrf-hydro-channel-ID-crosswalk-file-field-name",
        dest="wrf_hydro_channel_ID_crosswalk_file_field_name",
        help="Name of the column providing the channel segment IDs in the channel crosswalk file",
        default="ID",
    )
    parser.add_argument(
        "--wrf-hydro-channel-restart-upstream-flow-field-name",
        dest="wrf_hydro_channel_restart_upstream_flow_field_name",
        help="Name of the column providing the upstream flow at the beginning of the simulation.",
        default="qlink1",
    )
    parser.add_argument(
        "--wrf-hydro-channel-restart-downstream-flow-field-name",
        dest="wrf_hydro_channel_restart_downstream_flow_field_name",
        help="Name of the column providing the downstream flow at the beginning of the simulation.",
        default="qlink2",
    )
    parser.add_argument(
        "--wrf-hydro-channel-restart-depth-flow-field-name",
        dest="wrf_hydro_channel_restart_depth_flow_field_name",
        help="Name of the column providing the depth of flow at the beginning of the simulation.",
        default="hlink",
    )
    ql_arg_group = parser.add_mutually_exclusive_group()
    ql_arg_group.add_argument(
        "--qlc",
        "--constant_qlateral",
        help="Constant qlateral to apply to all time steps at all segments",
        dest="qlat_const",
        type=float,
        default=10,
    )
    ql_arg_group.add_argument(
        "--qlf",
        "--single_file_qlateral",
        help="QLaterals arranged with segment IDs as rows and timesteps as columns in a single .csv",
        dest="qlat_input_file",
    )
    ql_arg_group.add_argument(
        "--qlw",
        "--ql_wrf_hydro_folder",
        help="QLaterals in separate netcdf files as found in standard WRF-Hydro output",
        dest="qlat_input_folder",
    )
    ql_arg_group.add_argument(
        "--qlic",
        "--qlat_file_index_col",
        help="QLateral index column number",
        dest="qlat_file_index_col",
        default="feature_id",
    )
    ql_arg_group.add_argument(
        "--qlvc",
        "--qlat_file_value_col",
        help="QLateral value column number",
        dest="qlat_file_value_col",
        default="q_lateral",
    )
    parser.add_argument(
        "--qlat_file_pattern_filter",
        help="Provide a globbing pattern to identify files in the Wrf-Hydro qlateral output file folder",
        dest="qlat_file_pattern_filter",
        default="q_lateral",
    )
    parser.add_argument("--ql", help="QLat input data", dest="ql", default=None)
    # TODO: uncomment custominput file
    # supernetwork_arg_group = parser.add_mutually_exclusive_group()
    # supernetwork_arg_group.add_argument(
    #     "-f",
    #     "--custom-input-file",
    #     dest="custom_input_file",
    #     help="OR... please enter the path of a .yaml or .json file containing a custom supernetwork information. See for example test/input/yaml/CustomInput.yaml and test/input/json/CustomInput.json.",
    # )
    return parser.parse_args()


ENV_IS_CL = False
if ENV_IS_CL:
    root = pathlib.Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = pathlib.Path("../..").resolve()
    # sys.path.append(r"../python_framework_v02")

    # TODO: automate compile for the package scripts
    sys.path.append("fast_reach")

## network and reach utilities
import troute.nhd_network_utilities_v02 as nnu
import mc_reach
import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io


def writetoFile(file, writeString):
    file.write(writeString)
    file.write("\n")


def constant_qlats(index_dataset, nsteps, qlat):
    q = np.full((len(index_dataset.index), nsteps), qlat, dtype="float32")
    ql = pd.DataFrame(q, index=index_dataset.index, columns=range(nsteps))
    return ql


def compute_nhd_routing_v02(
    reaches_bytw,
    compute_func,
    parallel_compute_method,
    cpu_pool,
    nts,
    qts_subdivisions,
    independent_networks,
    param_df,
    qlats,
    q0,
    assume_short_ts,
):

    if parallel_compute_method == "by-network":
        with Parallel(n_jobs=cpu_pool, backend="threading") as parallel:
            jobs = []
            for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
                r = list(chain.from_iterable(reach_list))
                param_df_sub = param_df.loc[
                    r, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
                ].sort_index()
                qlat_sub = qlats.loc[r].sort_index()
                q0_sub = q0.loc[r].sort_index()
                # da_sub = da.loc[r].sort_index()
                jobs.append(
                    delayed(compute_func)(
                        nts,
                        qts_subdivisions,
                        reach_list,
                        independent_networks[tw],
                        param_df_sub.index.values,
                        param_df_sub.columns.values,
                        param_df_sub.values,
                        qlat_sub.values,
                        q0_sub.values,
                        # da_sub.values,
                        assume_short_ts,
                    )
                )
            results = parallel(jobs)

    else:  # Execute in serial
        results = []
        for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
            r = list(chain.from_iterable(reach_list))
            param_df_sub = param_df.loc[
                r, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
            ].sort_index()
            qlat_sub = qlats.loc[r].sort_index()
            q0_sub = q0.loc[r].sort_index()
            # print(q0_sub)
            # print(qlat_sub)
            # da_sub = da.loc[r].sort_index()
            results.append(
                compute_func(
                    nts,
                    qts_subdivisions,
                    reach_list,
                    independent_networks[tw],
                    param_df_sub.index.values,
                    param_df_sub.columns.values,
                    param_df_sub.values,
                    qlat_sub.values,
                    q0_sub.values,
                    # da_sub.values,
                    assume_short_ts,
                )
            )
        # print(results)
    return results
   


def _input_handler():

    args = _handle_args()

    custom_input_file = args.custom_input_file
    supernetwork_parameters = None
    waterbody_parameters = {}
    forcing_parameters = {}
    restart_parameters = {}
    output_parameters = {}
    run_parameters = {}
    data_assimilation_parameters = {}

    if custom_input_file:
        (
            supernetwork_parameters,
            waterbody_parameters,
            forcing_parameters,
            restart_parameters,
            output_parameters,
            run_parameters,
            data_assimilation_parameters,
        ) = nhd_io.read_custom_input(custom_input_file)
        # TODO: uncomment custominput file
        #     qlat_const = forcing_parameters.get("qlat_const", None)
        #     qlat_input_file = forcing_parameters.get("qlat_input_file", None)
        #     qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
        #     qlat_file_pattern_filter = forcing_parameters.get(
        #         "qlat_file_pattern_filter", None
        #     )
        #     qlat_file_index_col = forcing_parameters.get("qlat_file_index_col", None)
        #     qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", None)
        # else:
        # TODO: uncomment custominput file

    else:
        run_parameters["assume_short_ts"] = args.assume_short_ts
        run_parameters["parallel_compute_method"] = args.parallel_compute_method
        run_parameters["cpu_pool"] = args.cpu_pool
        run_parameters["showtiming"] = args.showtiming

        run_parameters["debuglevel"] = debuglevel = -1 * args.debuglevel
        run_parameters["verbose"] = verbose = args.verbose
        
        test_folder = pathlib.Path(root, "test")
        geo_input_folder = test_folder.joinpath("input", "geo")

        run_pocono2_test = args.run_pocono2_test
        run_pocono1_test = args.run_pocono1_test

        if run_pocono2_test:
            if verbose:
                print("running test case for Pocono_TEST2 domain")
            # Overwrite the following test defaults
            supernetwork = "Pocono_TEST2"
            waterbody_parameters["break_network_at_waterbodies"] = False
            run_parameters["qts_subdivisions"] = qts_subdivisions = 1
            run_parameters["dt"] = 300 / qts_subdivisions
            run_parameters["nts"] = 144 * qts_subdivisions
            output_parameters["csv_output"] = {
                "csv_output_folder": os.path.join(root, "test", "output", "text")
            }
            output_parameters["nc_output_folder"] = os.path.join(
                root, "test", "output", "text"
            )
            # test 1. Take lateral flow from re-formatted wrf-hydro output from Pocono Basin simulation
            forcing_parameters["qlat_input_file"] = os.path.join(
                root, r"test/input/geo/PoconoSampleData2/Pocono_ql_testsamp1_nwm_mc.csv"
            )

        elif run_pocono1_test:
            # NOTE: The test case for the Pocono basin was derived from this
            # resource on HydroShare, developed by aaraney and sourced from the
            # wrf_hydro_nwm_public repository on GitHub
            # see: https://www.hydroshare.org/resource/03ca354200e540018d44183598890448/
            # By downloading aaraney's docker job scheduler repo from GitHub, one can
            # execute the WRF-Hydro model that generated the test results
            # see: https://github.com/aaraney/NWM-Dockerized-Job-Scheduler
            if verbose:
                print("running test case for Pocono_TEST1 domain")
            # Overwrite the following test defaults

            NWM_test_path = os.path.join(
                root, "test/input/geo/NWM_2.1_Sample_Datasets/Pocono_TEST1/"
            )
            # lakeparm_file = os.path.join(
            #     NWM_test_path, "primary_domain", "DOMAIN", "LAKEPARM.nc",
            # )
            routelink_file = os.path.join(
                NWM_test_path, "primary_domain", "DOMAIN", "Route_Link.nc",
            )
            time_string = "2017-12-31_06-00_DOMAIN1"
            wrf_hydro_restart_file = os.path.join(
                NWM_test_path, "example_RESTART", "HYDRO_RST." + time_string
            )
            supernetwork_parameters = {
                "title_string": "Custom Input Example (using Pocono Test Example datafile)",
                "geo_file_path": routelink_file,
                "columns": {
                    "key": "link",
                    "downstream": "to",
                    "dx": "Length",
                    "n": "n",  # TODO: rename to `manningn`
                    "ncc": "nCC",  # TODO: rename to `mannningncc`
                    "s0": "So",  # TODO: rename to `bedslope`
                    "bw": "BtmWdth",  # TODO: rename to `bottomwidth`
                    "waterbody": "NHDWaterbodyComID",
                    "tw": "TopWdth",  # TODO: rename to `topwidth`
                    "twcc": "TopWdthCC",  # TODO: rename to `topwidthcc`
                    "musk": "MusK",
                    "musx": "MusX",
                    "cs": "ChSlp",  # TODO: rename to `sideslope`
                },
                "waterbody_null_code": -9999,
                "terminal_code": 0,
                "driver_string": "NetCDF",
                "layer_string": 0,
            }
            # waterbody_parameters = {
            #     "level_pool": {
            #         "level_pool_waterbody_parameter_file_path": lakeparm_file,
            #         "level_pool_waterbody_id": "lake_id",
            #         "level_pool_waterbody_area": "LkArea",
            #         "level_pool_weir_elevation": "WeirE",
            #         "level_pool_waterbody_max_elevation": "LkMxE",
            #         "level_pool_outfall_weir_coefficient": "WeirC",
            #         "level_pool_outfall_weir_length": "WeirL",
            #         "level_pool_overall_dam_length": "DamL",
            #         "level_pool_orifice_elevation": "OrificeE",
            #         "level_pool_orifice_coefficient": "OrificeC",
            #         "level_pool_orifice_area": "OrificeA",
            #     }
            # }
            # break_network_at_waterbodies = True
            run_parameters["qts_subdivisions"] = qts_subdivisions = 12
            run_parameters["dt"] = 3600 / qts_subdivisions
            run_parameters["nts"] = 24 * qts_subdivisions
            data_assimilation_parameters["wrf_hydro_channel_ID_routelink_file"] = routelink_subset_file
            output_parameters["csv_output"] = None
            output_parameters["nc_output_folder"] = None
            # build a time string to specify input date
            restart_parameters["wrf_hydro_channel_restart_file"] = wrf_hydro_restart_file
            restart_parameters["wrf_hydro_channel_ID_crosswalk_file"] = routelink_file
            restart_parameters[
                "wrf_hydro_channel_ID_crosswalk_file_field_name"
            ] = "link"
            restart_parameters[
                "wrf_hydro_channel_restart_upstream_flow_field_name"
            ] = "qlink1"
            restart_parameters[
                "wrf_hydro_channel_restart_downstream_flow_field_name"
            ] = "qlink2"
            restart_parameters[
                "wrf_hydro_channel_restart_depth_flow_field_name"
            ] = "hlink"
            # restart_parameters["wrf_hydro_waterbody_restart_file"] = wrf_hydro_restart_file
            # restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file"] = lakeparm_file
            # restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file_field_name"] = "lake_id"
            # restart_parameters["wrf_hydro_waterbody_crosswalk_filter_file"] = routelink_file
            # restart_parameters["wrf_hydro_waterbody_crosswalk_filter_file_field_name"] = "NHDWaterbodyComID"
            # restart_parameters["wrf_hydro_waterbody_crosswalk_file_output_order_field= "AscendingIndex"
            forcing_parameters["qlat_input_folder"] = os.path.join(
                root,
                "test/input/geo/NWM_2.1_Sample_Datasets/Pocono_TEST1/example_CHRTOUT/",
            )
            forcing_parameters["qlat_file_pattern_filter"] = "/*.CHRTOUT_DOMAIN1"
            forcing_parameters["qlat_file_index_col"] = "feature_id"
            forcing_parameters["qlat_file_value_col"] = "q_lateral"

        else:
            run_parameters["dt"] = args.dt
            run_parameters["nts"] = args.nts
            run_parameters["qts_subdivisions"] = args.qts_subdivisions
            run_parameters["compute_method"] = args.compute_method
            waterbody_parameters[
                "break_network_at_waterbodies"
            ] = args.break_network_at_waterbodies
            output_parameters["csv_output_folder"] = args.csv_output_folder
            data_assimilation_parameters["wrf_hydro_channel_ID_routelink_file"] = args.wrf_hydro_channel_ID_crosswalk_file
            restart_parameters[
                "wrf_hydro_channel_restart_file"
            ] = args.wrf_hydro_channel_restart_file
            restart_parameters[
                "wrf_hydro_channel_ID_crosswalk_file"
            ] = args.wrf_hydro_channel_ID_crosswalk_file
            restart_parameters[
                "wrf_hydro_channel_ID_crosswalk_file_field_name"
            ] = args.wrf_hydro_channel_ID_crosswalk_file_field_name
            restart_parameters[
                "wrf_hydro_channel_restart_upstream_flow_field_name"
            ] = args.wrf_hydro_channel_restart_upstream_flow_field_name
            restart_parameters[
                "wrf_hydro_channel_restart_downstream_flow_field_name"
            ] = args.wrf_hydro_channel_restart_downstream_flow_field_name
            restart_parameters[
                "wrf_hydro_channel_restart_depth_flow_field_name"
            ] = args.wrf_hydro_channel_restart_depth_flow_field_name

            forcing_parameters["qlat_const"] = float(args.qlat_const)
            forcing_parameters["qlat_input_folder"] = args.qlat_input_folder
            forcing_parameters["qlat_input_file"] = args.qlat_input_file
            forcing_parameters[
                "qlat_file_pattern_filter"
            ] = args.qlat_file_pattern_filter
            forcing_parameters["qlat_file_index_col"] = args.qlat_file_index_col
            forcing_parameters["qlat_file_value_col"] = args.qlat_file_value_col

            supernetwork = args.supernetwork

        # STEP 0.5: Obtain Supernetwork Parameters for test cases
        if not supernetwork_parameters:
            supernetwork_parameters = nnu.set_supernetwork_parameters(
                supernetwork=supernetwork,
                geo_input_folder=geo_input_folder,
                verbose=False,
                debuglevel=debuglevel,
            )

    return (
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        data_assimilation_parameters,
    )


def main():

    (
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        data_assimilation_parameters,
    ) = _input_handler()

    dt = run_parameters.get("dt", None)
    nts = run_parameters.get("nts", None)
    verbose = run_parameters.get("verbose", None)
    showtiming = run_parameters.get("showtiming", None)
    debuglevel = run_parameters.get("debuglevel", 0)

    if verbose:
        print("creating supernetwork connections set")
    if showtiming:
        start_time = time.time()

    # STEP 1: Build basic network connections graph
    connections, wbodies, param_df = nnu.build_connections(supernetwork_parameters, dt)

    if verbose:
        print("supernetwork connections set complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 2: Identify Independent Networks and Reaches by Network
    if showtiming:
        start_time = time.time()
    if verbose:
        print("organizing connections into reaches ...")

    independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections
    )

    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 4: Handle Channel Initial States
    if showtiming:
        start_time = time.time()
    if verbose:
        print("setting channel initial states ...")

    q0 = nnu.build_channel_initial_state(restart_parameters, param_df.index)
    
    if verbose:
        print("channel initial states complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))
        start_time = time.time()

    # STEP 5: Read (or set) QLateral Inputs
    if showtiming:
        start_time = time.time()
    if verbose:
        print("creating qlateral array ...")

    qlats = nnu.build_qlateral_array(forcing_parameters, connections.keys(), nts)

    if verbose:
        print("qlateral array complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    #STEP 6

    
    usgs_timeslices_folder = os.path.join(
        root,
        "test/input/geo/nudgingTimeSliceObs/",
    )
    routelink_subset_folder = os.path.join(
        root,
        "test/input/geo/routelink/",
    )
    # "../../test/input/geo/routelink/routeLink_subset.nc"
    routelink_subset_file = routelink_subset_folder+"routeLink_subset.nc"
    data_assimilation_parameters["wrf_hydro_channel_ID_routelink_file"] = routelink_subset_file
    # usgs_file_pattern_filter = "*.usgsTimeSlice.ncdf"

    # usgs_files = glob.glob(usgs_timeslices_folder + usgs_file_pattern_filter)
    # file_name = "2020-03-19_18:00:00.15min.usgsTimeSlice.ncdf"

    usgs_df = nhd_io.get_usgs_from_wrf_hydro(data_assimilation_parameters["wrf_hydro_channel_ID_routelink_file"],
    usgs_timeslices_folder,
    )

    print(usgs_df)
    # da = nnu.build_channel_initial_state(data_assimilation_parameters["wrf_hydro_channel_ID_routelink_file"], usgs_df.index)
    ################### Main Execution Loop across ordered networks
    if showtiming:
        main_start_time = time.time()
    if verbose:
        print(f"executing routing computation ...")

    if run_parameters.get("compute_method", None) == "standard cython compute network":
        compute_func = mc_reach.compute_network
    else:
        compute_func = mc_reach.compute_network

    results = compute_nhd_routing_v02(
        reaches_bytw,
        compute_func,
        run_parameters.get("parallel_compute_method", None),
        run_parameters.get("cpu_pool", None),
        run_parameters.get("nts", 1),
        run_parameters.get("qts_subdivisions", 1),
        independent_networks,
        param_df,
        qlats,
        q0,
        run_parameters.get("assume_short_ts", False),
    )

    csv_output_folder = output_parameters.get("csv_output_folder", None)
    if (debuglevel <= -1) or csv_output_folder:
        qvd_columns = pd.MultiIndex.from_product(
            [range(nts), ["q", "v", "d"]]
        ).to_flat_index()
        flowveldepth = pd.concat(
            [pd.DataFrame(d, index=i, columns=qvd_columns) for i, d in results],
            copy=False,
        )

        if csv_output_folder:
            flowveldepth = flowveldepth.sort_index()
            output_path = pathlib.Path(csv_output_folder).resolve()
            flowveldepth.to_csv(output_path.joinpath(f"{args.supernetwork}.csv"))

        if debuglevel <= -1:
            print(flowveldepth)

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # print(q0)
if __name__ == "__main__":
    main()
