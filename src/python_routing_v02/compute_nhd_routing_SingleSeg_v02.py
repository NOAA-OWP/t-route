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
import xarray as xr


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
            "Florence_FULL_RES"
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


def main():
    args = _handle_args()

    custom_input_file = args.custom_input_file

    if custom_input_file:
        (
            supernetwork_parameters,
            waterbody_parameters,
            forcing_parameters,
            restart_parameters,
            output_parameters,
            run_parameters,
        ) = nio.read_custom_input(custom_input_file)
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
        dt = args.dt
        nts = args.nts
        qts_subdivisions = args.qts_subdivisions
        debuglevel = -1 * args.debuglevel
        verbose = args.verbose
        showtiming = args.showtiming
        supernetwork = args.supernetwork
        break_network_at_waterbodies = args.break_network_at_waterbodies
        csv_output_folder = args.csv_output_folder
        assume_short_ts = args.assume_short_ts
        test_folder = pathlib.Path(root, "test")
        geo_input_folder = test_folder.joinpath("input", "geo")

        parallel_compute_method = args.parallel_compute_method

        cpu_pool = args.cpu_pool
        compute_method = args.compute_method
        wrf_hydro_channel_restart_file = args.wrf_hydro_channel_restart_file
        wrf_hydro_channel_ID_crosswalk_file = args.wrf_hydro_channel_ID_crosswalk_file
        wrf_hydro_channel_ID_crosswalk_file_field_name = (
            args.wrf_hydro_channel_ID_crosswalk_file_field_name
        )
        wrf_hydro_channel_restart_upstream_flow_field_name = (
            args.wrf_hydro_channel_restart_upstream_flow_field_name
        )
        wrf_hydro_channel_restart_downstream_flow_field_name = (
            args.wrf_hydro_channel_restart_downstream_flow_field_name
        )
        wrf_hydro_channel_restart_depth_flow_field_name = (
            args.wrf_hydro_channel_restart_depth_flow_field_name
        )
        qlat_const = float(args.qlat_const)
        qlat_input_folder = args.qlat_input_folder
        qlat_input_file = args.qlat_input_file
        qlat_file_pattern_filter = args.qlat_file_pattern_filter
        qlat_file_index_col = args.qlat_file_index_col
        qlat_file_value_col = args.qlat_file_value_col
        # print(forcing_parameters,qlat_const)
        # TODO: Make these commandline args
        """##NHD Subset (Brazos/Lower Colorado)"""
        # supernetwork = 'Brazos_LowerColorado_Named_Streams'
        # supernetwork = 'Brazos_LowerColorado_ge5'
        # supernetwork = 'Pocono_TEST1'
        """##NHD CONUS order 5 and greater"""
        # supernetwork = 'CONUS_ge5'
        """These are large -- be careful"""
        # supernetwork = 'Mainstems_CONUS'
        # supernetwork = 'CONUS_FULL_RES_v20'
        # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
        # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches

    run_pocono2_test = args.run_pocono2_test
    run_pocono1_test = args.run_pocono1_test

    if run_pocono2_test:
        if verbose:
            print("running test case for Pocono_TEST2 domain")
        # Overwrite the following test defaults
        supernetwork = "Pocono_TEST2"
        break_network_at_waterbodies = False
        qts_subdivisions = 1  # change qts_subdivisions = 1 as  default
        dt = 300 / qts_subdivisions
        nts = 144 * qts_subdivisions
        csv_output = {"csv_output_folder": os.path.join(root, "test", "output", "text")}
        nc_output_folder = os.path.join(root, "test", "output", "text")
        # test 1. Take lateral flow from re-formatted wrf-hydro output from Pocono Basin simulation
        qlat_input_file = os.path.join(
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

        supernetwork = "Pocono_TEST2"
        assume_short_ts = True
        
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
            "cols_as_text": False,
            "key_col": 16,  # "link",
            "downstream_col": 22,  # "to",
            "length_col": 3,  # "Length",
            "manningn_col": 18,  # "n",
            "manningncc_col": 19,  # "nCC",
            "slope_col": 8,  # "So",
            "bottomwidth_col": 0,  # "BtmWdth",
            "topwidth_col": 9,  # "TopWdth",
            "topwidthcc_col": 10,  # "TopWdthCC",
            "waterbody_col": 6,  # "NHDWaterbodyComID",
            "waterbody_null_code": -9999,
            "MusK_col": 4,  # "MusK",
            "MusX_col": 5,  # "MusX",
            "ChSlp_col": 1,  # "ChSlp",
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
        csv_output = None
        nc_output_folder = None
        # build a time string to specify input date
        wrf_hydro_channel_restart_file = wrf_hydro_restart_file
        wrf_hydro_channel_ID_crosswalk_file = routelink_file
        wrf_hydro_channel_ID_crosswalk_file_field_name = "link"
        wrf_hydro_channel_restart_upstream_flow_field_name = "qlink1"
        wrf_hydro_channel_restart_downstream_flow_field_name = "qlink2"
        wrf_hydro_channel_restart_depth_flow_field_name = "hlink"
        # wrf_hydro_waterbody_restart_file = wrf_hydro_restart_file
        # wrf_hydro_waterbody_ID_crosswalk_file = lakeparm_file
        # wrf_hydro_waterbody_ID_crosswalk_file_field_name = "lake_id"
        # wrf_hydro_waterbody_crosswalk_filter_file = routelink_file
        # wrf_hydro_waterbody_crosswalk_filter_file_field_name = "NHDWaterbodyComID"
        # wrf_hydro_waterbody_crosswalk_file_output_order_field= "AscendingIndex"
        qlat_input_folder = os.path.join(
            root, "test/input/geo/NWM_2.1_Sample_Datasets/Pocono_TEST1/example_CHRTOUT/"
        )
        qlat_file_pattern_filter = "/*.CHRTOUT_DOMAIN1"
        qlat_file_index_col = "feature_id"
        qlat_file_value_col = "q_lateral"
        
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # time domain - develop from temporal extent of qlat data
        qlat_files = glob.glob(qlat_input_folder + qlat_file_pattern_filter, recursive=True)
        ql = nhd_io.get_ql_from_wrf_hydro(qlat_files,qlat_file_index_col)
        
        wrf_time = ql.columns.astype("datetime64[ns]")
        dt_wrf = (wrf_time[1] - wrf_time[0])
        sim_duration = (wrf_time[-1] + dt_wrf) - wrf_time[0]
        
        dt = 300.0  # routing simulation timestep (seconds) 
        dt_routing = pd.Timedelta(str(dt) + 'seconds')
        nts = round(sim_duration / dt_routing)
        qts_subdivisions = dt_wrf/dt_routing
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # get wrf-simulated flows
        flow_wrf = nhd_io.get_ql_from_wrf_hydro(qlat_files,qlat_file_index_col, "streamflow")
#         li = []
#         for filename in qlat_files:
#             with xr.open_dataset(filename) as ds:
#                 df1 = ds.to_dataframe()
#             li.append(df1)

#         frame = pd.concat(li, axis=0, ignore_index=False)
#         mod = frame.reset_index()
#         flow_wrf = mod.pivot(index=qlat_file_index_col, columns="time", values="streamflow")
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    if verbose:
        print("creating supernetwork connections set")
    if showtiming:
        start_time = time.time()

    # STEP 1
    network_data = nnu.set_supernetwork_data(
        supernetwork=supernetwork,
        geo_input_folder=geo_input_folder,
        verbose=False,
        debuglevel=debuglevel,
    )

    cols = network_data["columns"]
    param_df = nhd_io.read(network_data["geo_file_path"])
    param_df = param_df[list(cols.values())]
    param_df = param_df.set_index(cols["key"])

    if "mask_file_path" in network_data:
        data_mask = nhd_io.read_mask(
            network_data["mask_file_path"],
            layer_string=network_data["mask_layer_string"],
        )
        param_df = param_df.filter(data_mask.iloc[:, network_data["mask_key"]], axis=0)

    param_df = param_df.sort_index()
    param_df = nhd_io.replace_downstreams(param_df, cols["downstream"], 0)

    connections = nhd_network.extract_connections(param_df, cols["downstream"])
    wbodies = nhd_network.extract_waterbodies(
        param_df, cols["waterbody"], network_data["waterbody_null_code"]
    )

    #     # initial conditions, assume to be zero
    #     # TODO: Allow optional reading of initial conditions from WRF
    #     q0 = pd.DataFrame(
    #         0, index=param_df.index, columns=["qu0", "qd0", "h0"], dtype="float32"
    #     )

    if verbose:
        print("supernetwork connections set complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming:
        start_time = time.time()
    if verbose:
        print("organizing connections into reaches ...")

    rconn = nhd_network.reverse_network(connections)
    independent_networks = nhd_network.reachable_network(rconn)
    reaches_bytw = {}
    for tw, net in independent_networks.items():
        path_func = partial(nhd_network.split_at_junction, net)
        reaches_bytw[tw] = nhd_network.dfs_decomposition(net, path_func)

    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    if showtiming:
        start_time = time.time()

    param_df["dt"] = dt
    param_df = param_df.rename(columns=nnu.reverse_dict(cols))
    param_df = param_df.astype("float32")

    # datasub = data[['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']]

    # STEP 4: Handle Channel Initial States
    if showtiming:
        start_time = time.time()
    if verbose:
        print("setting channel initial states ...")

    if wrf_hydro_channel_restart_file:

        q0 = nhd_io.get_stream_restart_from_wrf_hydro(
            wrf_hydro_channel_restart_file,
            wrf_hydro_channel_ID_crosswalk_file,
            wrf_hydro_channel_ID_crosswalk_file_field_name,
            wrf_hydro_channel_restart_upstream_flow_field_name,
            wrf_hydro_channel_restart_downstream_flow_field_name,
            wrf_hydro_channel_restart_depth_flow_field_name,
        )
    else:
        # Set cold initial state
        q0 = pd.DataFrame(
            0, index=connections.keys(), columns=["qu0", "qd0", "h0",], dtype="float32"
        )

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

    if qlat_input_folder:
        qlat_files = glob.glob(qlat_input_folder + qlat_file_pattern_filter)
        qlat_df = nhd_io.get_ql_from_wrf_hydro(
            qlat_files=qlat_files,
            index_col=qlat_file_index_col,
            value_col=qlat_file_value_col,
        )

        qlat_df = qlat_df[qlat_df.index.isin(connections.keys())]
        df_length = len(qlat_df.columns)

    # TODO: These three lines seem extraneous
    #        for x in range(df_length, 144):
    #            qlat_df[str(x)] = 0
    #            qlat_df = qlat_df.astype("float32")

    elif qlat_input_file:
        qlat_df = nhd_io.get_ql_from_csv(qlat_input_file)

    else:
        qlat_df = pd.DataFrame(
            qlat_const, index=connections.keys(), columns=range(nts), dtype="float32",
        )

    qlats = qlat_df

    if verbose:
        print("qlateral array complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))
        start_time = time.time()

    if compute_method == "standard cython compute network":
        compute_func = mc_reach.compute_network
    else:
        compute_func = mc_reach.compute_network

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
                    assume_short_ts,
                )
            )

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
        
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if verbose and run_pocono1_test:
        print("parity check against WRF Hydro data")
    
        # create a multi-index DataFrame with flow, depth, and velocity simulations
        fdv_columns = pd.MultiIndex.from_product([range(nts), ["q", "v", "d"]])
        flowveldepth = pd.concat(
            [pd.DataFrame(d, index=i, columns=fdv_columns) for i, d in results], copy=False
        )
        flowveldepth = flowveldepth.sort_index()

        flows = flowveldepth.loc[:, (slice(None), "q")]
        flows = flows.T.reset_index(level = [0,1])
        flows.rename(columns = {"level_0": "Timestep", "level_1": "Parameter"}, inplace = True)
        flows['Time (d)'] = ((flows.Timestep + 1) * dt)/(24*60*60)
        flows = flows.set_index('Time (d)')
        
        
        # specify calendar times for each simulation
        time_routing = pd.period_range((wrf_time[0]-dt_wrf+dt_routing), periods=nts, freq=dt_routing).astype("datetime64[ns]")
        time_wrf = flow_wrf.columns.values

        compare_node = 4186169
        
        # print flow values at common timesteps
        # create table of flow results at common timesteps
        trt = pd.DataFrame(flows.loc[:,compare_node ].values, index = time_routing, columns = ["flow, t-route (cms)"])
        wrf = pd.DataFrame(flow_wrf.loc[compare_node ,:].values, index = time_wrf, columns = ["flow, wrf (cms)"])

        compare = pd.concat([wrf, trt], axis=1, sort=False, join = 'inner')
        
        print(compare)


if __name__ == "__main__":
    main()
