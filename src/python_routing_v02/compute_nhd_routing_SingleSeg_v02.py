#!/usr/bin/env python
# coding: utf-8
# example usage: python compute_nhd_routing_SingleSeg.py -v -t -w -n Mainstems_CONUS


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
        "--nts",
        "--number-of-qlateral-timesteps",
        help="Set the number of timesteps to execute. If used with ql_file or ql_folder, nts must be less than len(ql) x qN.",
        dest="nts",
        default=144,
        type=int,
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
        "--compute_method",
        help="Use the cython version of the compute_network code (enter additional flag options for other compute_network possibilities).",
        dest="compute_method",
        default="standard cython compute network",
    )
    parser.add_argument(
        "-n",
        "--supernetwork",
        help="Choose from among the pre-programmed supernetworks (Pocono_TEST1, Pocono_TEST2, LowerColorado_Conchos_FULL_RES, Brazos_LowerColorado_ge5, Brazos_LowerColorado_FULL_RES, Brazos_LowerColorado_Named_Streams, CONUS_ge5, Mainstems_CONUS, CONUS_Named_Streams, CONUS_FULL_RES_v20",
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
    parser.add_argument(
        "--wrf_hydro_channel_restart_file",
        dest="wrf_hydro_channel_restart_file",
        help="provide a WRF-Hydro channel warm state file (may be the same as waterbody restart file)",
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

    nts = args.nts
    debuglevel = -1 * args.debuglevel
    verbose = args.verbose
    showtiming = args.showtiming
    supernetwork = args.supernetwork
    break_network_at_waterbodies = args.break_network_at_waterbodies
    csv_output_folder = args.csv_output_folder
    assume_short_ts = args.assume_short_ts
    # TODO: uncomment custominput file
    # custom_input_file = args.custom_input_file
    test_folder = pathlib.Path(root, "test")
    geo_input_folder = test_folder.joinpath("input", "geo")

    # TODO: uncomment custominput file
    # if custom_input_file:
    #     (
    #         supernetwork_parameters,
    #         waterbody_parameters,
    #         forcing_parameters,
    #         restart_parameters,
    #         output_parameters,
    #         run_parameters,
    #     ) = nhd_io.read_custom_input(custom_input_file)

    #     qlat_const = forcing_parameters.get("qlat_const", None)
    #     qlat_input_file = forcing_parameters.get("qlat_input_file", None)
    #     qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
    #     qlat_file_pattern_filter = forcing_parameters.get(
    #         "qlat_file_pattern_filter", None
    #     )
    #     qlat_file_index_col = forcing_parameters.get("qlat_file_index_col", None)
    #     qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", None)
    # else:
    wrf_hydro_channel_restart_file = args.wrf_hydro_channel_restart_file
    # TODO: uncomment custominput file
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

    if verbose:
        print("creating supernetwork connections set")
    if showtiming:
        start_time = time.time()

    # STEP 1
    network_data = nnu.set_supernetwork_data(
        supernetwork=args.supernetwork,
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

    # initial conditions, assume to be zero
    # TODO: Allow optional reading of initial conditions from WRF
    q0 = pd.DataFrame(
        0, index=param_df.index, columns=["qu0", "qd0", "h0"], dtype="float32"
    )

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

    param_df["dt"] = 300.0
    param_df = param_df.rename(columns=nnu.reverse_dict(cols))
    param_df = param_df.astype("float32")

    # datasub = data[['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']]
    
    # STEP 4: Handle Channel Initial States
    if showtiming:
        start_time = time.time()
    if verbose:
        print("setting channel initial states ...")

    if wrf_hydro_channel_restart_file:

        channel_initial_states_df = nhd_io.get_stream_restart_from_wrf_hydro(
            wrf_hydro_channel_restart_file,
            wrf_hydro_channel_ID_crosswalk_file,
            wrf_hydro_channel_ID_crosswalk_file_field_name,
            wrf_hydro_channel_restart_upstream_flow_field_name,
            wrf_hydro_channel_restart_downstream_flow_field_name,
            wrf_hydro_channel_restart_depth_flow_field_name,
        )
    else:
        # TODO: Consider adding option to read cold state from route-link file
        channel_initial_us_flow_const = 0.0
        channel_initial_ds_flow_const = 0.0
        channel_initial_depth_const = 0.0
        # Set initial states from cold-state
        channel_initial_states_df = pd.DataFrame(
            0, index=connections.keys(), columns=["qu0", "qd0", "h0",], dtype="float32"
        )
        channel_initial_states_df["qu0"] = channel_initial_us_flow_const
        channel_initial_states_df["qd0"] = channel_initial_ds_flow_const
        channel_initial_states_df["h0"] = channel_initial_depth_const
        channel_initial_states_df["index"] = range(len(channel_initial_states_df))

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

        for x in range(df_length, 144):
            qlat_df[str(x)] = 0
            qlat_df = qlat_df.astype("float32")

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

    parallel_compute_method = args.parallel_compute_method

    cpu_pool = args.cpu_pool
    compute_method = args.compute_method

    if compute_method == "standard cython compute network":
        compute_func = mc_reach.compute_network
    else:
        compute_func = mc_reach.compute_network

    if parallel_compute_method == "by-network":
        with Parallel(n_jobs=cpu_pool, backend="threading") as parallel:
            jobs = []
            for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
                r = list(chain.from_iterable(reach_list))
                channel_initial_states_df = param_df.loc[
                    r, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
                ].sort_index()
                qlat_sub = qlats.loc[r].sort_index()
                q0_sub = q0.loc[r].sort_index()
                jobs.append(
                    delayed(compute_func)(
                        nts,
                        reach_list,
                        independent_networks[tw],
                        channel_initial_states_df.index.values,
                        channel_initial_states_df.columns.values,
                        channel_initial_states_df.values,
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
            channel_initial_states_df = param_df.loc[
                r, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
            ].sort_index()
            qlat_sub = qlats.loc[r].sort_index()
            q0_sub = q0.loc[r].sort_index()
            results.append(
                compute_func(
                    nts,
                    reach_list,
                    independent_networks[tw],
                    channel_initial_states_df.index.values,
                    channel_initial_states_df.columns.values,
                    channel_initial_states_df.values,
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


if __name__ == "__main__":
    main()
