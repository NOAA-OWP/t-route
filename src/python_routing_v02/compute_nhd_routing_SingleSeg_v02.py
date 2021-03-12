#!/usr/bin/env python
# coding: utf-8
# example usage: python compute_nhd_routing_SingleSeg.py -v -t -w -n Mainstems_CONUS
# python compute_nhd_routing_SingleSeg_v02.py --test pocono1 -t -v --debuglevel 1
# python compute_nhd_routing_SingleSeg_v02.py -f ../../test/input/yaml/CustomInput.yaml
# python compute_nhd_routing_SingleSeg_v02.py -f ../../test/input/yaml/FlorenceBenchmark.yaml


# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
## Parallel execution
import os
import sys
import time
from datetime import datetime
import numpy as np
import argparse
import pathlib
import glob
import pandas as pd
from collections import defaultdict
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

    # change this so after --test, the user enters a test choice
    parser.add_argument(
        "--test",
        help="Select a test case, routing results will be compared against WRF hydro for parity",
        choices=["pocono1"],
        dest="test_case",
    )

    parser.add_argument(
        "--sts",
        "--assume-short-ts",
        help="Use the previous timestep value for upstream flow",
        dest="assume_short_ts",
        action="store_true",
    )
    parser.add_argument(
        "--courant",
        "--return-courant-metrics",
        help="Return Courant evaluation metrics for each segment/timestep",
        dest="return_courant",
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
        "--subnet-size",
        help="Set the target size (number of segments) for grouped subnetworks.",
        dest="subnetwork_target_size",
        default=-1,
        type=int,
    )
    parser.add_argument(
        "--cpu-pool",
        help="Assign the number of cores to multiprocess across.",
        dest="cpu_pool",
        type=int,
        default=-1,
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
    # TODO: Refine exclusivity of ql args (currently not going to accept more than one arg; more than one is needed for qlw, for instance.)
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
import build_tests  # TODO: Determine whether and how to incorporate this into setup.py


def writetoFile(file, writeString):
    file.write(writeString)
    file.write("\n")


def constant_qlats(index_dataset, nsteps, qlat):
    q = np.full((len(index_dataset.index), nsteps), qlat, dtype="float32")
    ql = pd.DataFrame(q, index=index_dataset.index, columns=range(nsteps))
    return ql


def compute_nhd_routing_v02(
    connections,
    rconn,
    reaches_bytw,
    compute_func,
    parallel_compute_method,
    subnetwork_target_size,
    cpu_pool,
    nts,
    qts_subdivisions,
    independent_networks,
    param_df,
    qlats,
    q0,
    assume_short_ts,
    return_courant,
):

    start_time = time.time()
    if parallel_compute_method == "by-subnetwork-jit-clustered":
        networks_with_subnetworks_ordered_jit = nhd_network.build_subnetworks(
            connections, rconn, subnetwork_target_size
        )
        subnetworks_only_ordered_jit = defaultdict(dict)
        subnetworks = defaultdict(dict)
        for tw, ordered_network in networks_with_subnetworks_ordered_jit.items():
            intw = independent_networks[tw]
            for order, subnet_sets in ordered_network.items():
                subnetworks_only_ordered_jit[order].update(subnet_sets)
                for subn_tw, subnetwork in subnet_sets.items():
                    subnetworks[subn_tw] = {k: intw[k] for k in subnetwork}

        reaches_ordered_bysubntw = defaultdict(dict)
        for order, ordered_subn_dict in subnetworks_only_ordered_jit.items():
            for subn_tw, subnet in ordered_subn_dict.items():
                conn_subn = {k: connections[k] for k in subnet if k in connections}
                rconn_subn = {k: rconn[k] for k in subnet if k in rconn}
                path_func = partial(nhd_network.split_at_junction, rconn_subn)
                reaches_ordered_bysubntw[order][
                    subn_tw
                ] = nhd_network.dfs_decomposition(rconn_subn, path_func)

        cluster_threshold = 0.65  # When a job has a total segment count 65% of the target size, compute it
        # Otherwise, keep adding reaches.

        reaches_ordered_bysubntw_clustered = defaultdict(dict)

        for order in subnetworks_only_ordered_jit:
            cluster = 0
            reaches_ordered_bysubntw_clustered[order][cluster] = {
                "segs": [],
                "upstreams": {},
                "tw": [],
                "subn_reach_list": [],
            }
            for twi, (subn_tw, subn_reach_list) in enumerate(
                reaches_ordered_bysubntw[order].items(), 1
            ):
                segs = list(chain.from_iterable(subn_reach_list))
                reaches_ordered_bysubntw_clustered[order][cluster]["segs"].extend(segs)
                reaches_ordered_bysubntw_clustered[order][cluster]["upstreams"].update(
                    subnetworks[subn_tw]
                )

                reaches_ordered_bysubntw_clustered[order][cluster]["tw"].append(subn_tw)
                reaches_ordered_bysubntw_clustered[order][cluster][
                    "subn_reach_list"
                ].extend(subn_reach_list)

                if (
                    len(reaches_ordered_bysubntw_clustered[order][cluster]["segs"])
                    >= cluster_threshold * subnetwork_target_size
                ) and (
                    twi
                    < len(reaches_ordered_bysubntw[order])
                    # i.e., we haven't reached the end
                    # TODO: perhaps this should be a while condition...
                ):
                    cluster += 1
                    reaches_ordered_bysubntw_clustered[order][cluster] = {
                        "segs": [],
                        "upstreams": {},
                        "tw": [],
                        "subn_reach_list": [],
                    }

        if 1 == 1:
            print("JIT Preprocessing time %s seconds." % (time.time() - start_time))
            print("starting Parallel JIT calculation")

        start_para_time = time.time()
        # if 1 == 1:
        with Parallel(n_jobs=cpu_pool, backend="threading") as parallel:
            results_subn = defaultdict(list)
            flowveldepth_interorder = {}

            for order in range(max(subnetworks_only_ordered_jit.keys()), -1, -1):
                jobs = []
                for cluster, clustered_subns in reaches_ordered_bysubntw_clustered[
                    order
                ].items():
                    segs = clustered_subns["segs"]
                    offnetwork_upstreams = set()
                    segs_set = set(segs)
                    for seg in segs:
                        for us in rconn[seg]:
                            if us not in segs_set:
                                offnetwork_upstreams.add(us)

                    segs.extend(offnetwork_upstreams)
                    param_df_sub = param_df.loc[
                        segs, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
                    ].sort_index()
                    if order < max(subnetworks_only_ordered_jit.keys()):
                        for us_subn_tw in offnetwork_upstreams:
                            subn_tw_sortposition = param_df_sub.index.get_loc(
                                us_subn_tw
                            )
                            flowveldepth_interorder[us_subn_tw][
                                "position_index"
                            ] = subn_tw_sortposition
                    qlat_sub = qlats.loc[segs].sort_index()
                    q0_sub = q0.loc[segs].sort_index()
                    subn_reach_list = clustered_subns["subn_reach_list"]
                    upstreams = clustered_subns["upstreams"]

                    # results_subn[order].append(
                    #     compute_func(
                    jobs.append(
                        delayed(compute_func)(
                            nts,
                            qts_subdivisions,
                            subn_reach_list,
                            upstreams,
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            qlat_sub.values,
                            q0_sub.values,
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
                            {
                                us: fvd
                                for us, fvd in flowveldepth_interorder.items()
                                if us in offnetwork_upstreams
                            },
                            assume_short_ts,
                            return_courant,
                        )
                    )
                results_subn[order] = parallel(jobs)

                if order > 0:  # This is not needed for the last rank of subnetworks
                    flowveldepth_interorder = {}
                    for ci, (cluster, clustered_subns) in enumerate(
                        reaches_ordered_bysubntw_clustered[order].items()
                    ):
                        for subn_tw in clustered_subns["tw"]:
                            # TODO: This index step is necessary because we sort the segment index
                            # TODO: I think there are a number of ways we could remove the sorting step
                            #       -- the binary search could be replaced with an index based on the known topology
                            flowveldepth_interorder[subn_tw] = {}
                            subn_tw_sortposition = (
                                results_subn[order][ci][0].tolist().index(subn_tw)
                            )
                            flowveldepth_interorder[subn_tw]["results"] = results_subn[
                                order
                            ][ci][1][subn_tw_sortposition]
                            # what will it take to get just the tw FVD values into an array to pass to the next loop?
                            # There will be an empty array initialized at the top of the loop, then re-populated here.
                            # we don't have to bother with populating it after the last group

        results = []
        for order in subnetworks_only_ordered_jit:
            results.extend(results_subn[order])

        if 1 == 1:
            print("PARALLEL TIME %s seconds." % (time.time() - start_para_time))

    elif parallel_compute_method == "by-subnetwork-jit":
        networks_with_subnetworks_ordered_jit = nhd_network.build_subnetworks(
            connections, rconn, subnetwork_target_size
        )
        subnetworks_only_ordered_jit = defaultdict(dict)
        subnetworks = defaultdict(dict)
        for tw, ordered_network in networks_with_subnetworks_ordered_jit.items():
            intw = independent_networks[tw]
            for order, subnet_sets in ordered_network.items():
                subnetworks_only_ordered_jit[order].update(subnet_sets)
                for subn_tw, subnetwork in subnet_sets.items():
                    subnetworks[subn_tw] = {k: intw[k] for k in subnetwork}

        reaches_ordered_bysubntw = defaultdict(dict)
        for order, ordered_subn_dict in subnetworks_only_ordered_jit.items():
            for subn_tw, subnet in ordered_subn_dict.items():
                conn_subn = {k: connections[k] for k in subnet if k in connections}
                rconn_subn = {k: rconn[k] for k in subnet if k in rconn}
                path_func = partial(nhd_network.split_at_junction, rconn_subn)
                reaches_ordered_bysubntw[order][
                    subn_tw
                ] = nhd_network.dfs_decomposition(rconn_subn, path_func)

        if 1 == 1:
            print("JIT Preprocessing time %s seconds." % (time.time() - start_time))
            print("starting Parallel JIT calculation")

        start_para_time = time.time()
        with Parallel(n_jobs=cpu_pool, backend="threading") as parallel:
            results_subn = defaultdict(list)
            flowveldepth_interorder = {}

            for order in range(max(subnetworks_only_ordered_jit.keys()), -1, -1):
                jobs = []
                for twi, (subn_tw, subn_reach_list) in enumerate(
                    reaches_ordered_bysubntw[order].items(), 1
                ):
                    # TODO: Confirm that a list here is best -- we are sorting,
                    # so a set might be sufficient/better
                    segs = list(chain.from_iterable(subn_reach_list))
                    offnetwork_upstreams = set()
                    segs_set = set(segs)
                    for seg in segs:
                        for us in rconn[seg]:
                            if us not in segs_set:
                                offnetwork_upstreams.add(us)

                    segs.extend(offnetwork_upstreams)
                    param_df_sub = param_df.loc[
                        segs, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
                    ].sort_index()
                    if order < max(subnetworks_only_ordered_jit.keys()):
                        for us_subn_tw in offnetwork_upstreams:
                            subn_tw_sortposition = param_df_sub.index.get_loc(
                                us_subn_tw
                            )
                            flowveldepth_interorder[us_subn_tw][
                                "position_index"
                            ] = subn_tw_sortposition
                    qlat_sub = qlats.loc[segs].sort_index()
                    q0_sub = q0.loc[segs].sort_index()
                    jobs.append(
                        delayed(compute_func)(
                            nts,
                            qts_subdivisions,
                            subn_reach_list,
                            subnetworks[subn_tw],
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            qlat_sub.values,
                            q0_sub.values,
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
                            {
                                us: fvd
                                for us, fvd in flowveldepth_interorder.items()
                                if us in offnetwork_upstreams
                            },
                            assume_short_ts,
                            return_courant,
                        )
                    )

                results_subn[order] = parallel(jobs)

                if order > 0:  # This is not needed for the last rank of subnetworks
                    flowveldepth_interorder = {}
                    for twi, subn_tw in enumerate(reaches_ordered_bysubntw[order]):
                        # TODO: This index step is necessary because we sort the segment index
                        # TODO: I think there are a number of ways we could remove the sorting step
                        #       -- the binary search could be replaced with an index based on the known topology
                        flowveldepth_interorder[subn_tw] = {}
                        subn_tw_sortposition = (
                            results_subn[order][twi][0].tolist().index(subn_tw)
                        )
                        flowveldepth_interorder[subn_tw]["results"] = results_subn[
                            order
                        ][twi][1][subn_tw_sortposition]
                        # what will it take to get just the tw FVD values into an array to pass to the next loop?
                        # There will be an empty array initialized at the top of the loop, then re-populated here.
                        # we don't have to bother with populating it after the last group

        results = []
        for order in subnetworks_only_ordered_jit:
            results.extend(results_subn[order])

        if 1 == 1:
            print("PARALLEL TIME %s seconds." % (time.time() - start_para_time))

    elif parallel_compute_method == "by-network":
        with Parallel(n_jobs=cpu_pool, backend="threading") as parallel:
            jobs = []
            for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
                segs = list(chain.from_iterable(reach_list))
                param_df_sub = param_df.loc[
                    segs, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
                ].sort_index()
                qlat_sub = qlats.loc[segs].sort_index()
                q0_sub = q0.loc[segs].sort_index()
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
                        {},
                        assume_short_ts,
                        return_courant,
                    )
                )
            results = parallel(jobs)

    else:  # Execute in serial
        results = []
        # import pdb; pdb.set_trace()
        for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
            
            segs = list(chain.from_iterable(reach_list))
            param_df_sub = param_df.loc[
                segs, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
            ].sort_index()
            qlat_sub = qlats.loc[segs].sort_index()
            q0_sub = q0.loc[segs].sort_index()
            # import pdb; pdb.set_trace()
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
                    {},
                    assume_short_ts,
                    return_courant,
                )
            )

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
    parity_parameters = {}

    if custom_input_file:
        (
            supernetwork_parameters,
            waterbody_parameters,
            forcing_parameters,
            restart_parameters,
            output_parameters,
            run_parameters,
            parity_parameters,
        ) = nhd_io.read_custom_input(custom_input_file)
        run_parameters["debuglevel"] *= -1

    else:
        run_parameters["assume_short_ts"] = args.assume_short_ts
        run_parameters["return_courant"] = args.return_courant
        run_parameters["parallel_compute_method"] = args.parallel_compute_method
        run_parameters["subnetwork_target_size"] = args.subnetwork_target_size
        run_parameters["cpu_pool"] = args.cpu_pool
        run_parameters["showtiming"] = args.showtiming

        run_parameters["debuglevel"] = debuglevel = -1 * args.debuglevel
        run_parameters["verbose"] = verbose = args.verbose

        output_parameters["csv_output"] = {}
        output_parameters["csv_output"]["csv_output_folder"] = args.csv_output_folder

        test_folder = pathlib.Path(root, "test")
        geo_input_folder = test_folder.joinpath("input", "geo")

        test_case = args.test_case

        if test_case:

            # call test case assemble function
            (
                supernetwork_parameters,
                run_parameters,
                output_parameters,
                restart_parameters,
                forcing_parameters,
                parity_parameters,
            ) = build_tests.build_test_parameters(
                test_case,
                supernetwork_parameters,
                run_parameters,
                output_parameters,
                restart_parameters,
                forcing_parameters,
                parity_parameters,
            )

        else:
            run_parameters["dt"] = args.dt
            run_parameters["nts"] = args.nts
            run_parameters["qts_subdivisions"] = args.qts_subdivisions
            run_parameters["compute_method"] = args.compute_method

            waterbody_parameters[
                "break_network_at_waterbodies"
            ] = args.break_network_at_waterbodies

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
        parity_parameters,
    )
total_timeslice_files = 24
runs_to_be_completed = 3
ts_iterator = 0
restart_file_number = 2
file_run_size = int(total_timeslice_files/runs_to_be_completed)
def main():
    (
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        parity_parameters,
    ) = _input_handler()

    dt = run_parameters.get("dt", None)
    nts = run_parameters.get("nts", None)
    verbose = run_parameters.get("verbose", None)
    showtiming = run_parameters.get("showtiming", None)
    debuglevel = run_parameters.get("debuglevel", 0)
    global ts_iterator
    global restart_file_number
    global runs_to_be_completed
    global file_run_size
    # import pdb; pdb.set_trace()
    print(restart_file_number)
    if showtiming:
        main_start_time = time.time()

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

    if ts_iterator == 0:
        q0 = nnu.build_channel_initial_state(restart_parameters, supernetwork_parameters, param_df.index)
    
    else:
        q0_file_name = restart_parameters["wrf_hydro_channel_restart_file"][:-1] + str(ts_iterator+1) + ".csv" 
        q0 = nnu.restart_file_csv(q0_file_name)

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
    qlats = nnu.build_qlateral_array(
        forcing_parameters,
        connections.keys(),
        nts,
        ts_iterator,
        file_run_size,
        supernetwork_parameters,
        run_parameters.get("qts_subdivisions", 1),
    )

    if verbose:
        print("qlateral array complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    ################### Main Execution Loop across ordered networks
    if showtiming:
        start_time = time.time()
    if verbose:
        if run_parameters.get("return_courant", False):
            print(
                f"executing routing computation, with Courant evaluation metrics returned"
            )
        else:
            print(f"executing routing computation ...")

    if run_parameters.get("compute_method", None) == "standard cython compute network":
        compute_func = mc_reach.compute_network
    else:
        compute_func = mc_reach.compute_network
    # import pdb; pdb.set_trace()
    results = compute_nhd_routing_v02(
        connections,
        rconn,
        reaches_bytw,
        compute_func,
        run_parameters.get("parallel_compute_method", None),
        run_parameters.get("subnetwork_target_size", 1),
        # The default here might be the whole network or some percentage...
        run_parameters.get("cpu_pool", None),
        run_parameters.get("nts", 1),
        run_parameters.get("qts_subdivisions", 1),
        independent_networks,
        param_df,
        qlats,
        q0,
        run_parameters.get("assume_short_ts", False),
        run_parameters.get("return_courant", False),
    )

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    ################### Output Handling

    if showtiming:
        start_time = time.time()
    if verbose:
        print(f"Handling output ...")

    if output_parameters["csv_output"]:
        csv_output_folder = output_parameters["csv_output"].get(
            "csv_output_folder", None
        )

    if (debuglevel <= -1) or output_parameters["csv_output"]:

        qvd_columns = pd.MultiIndex.from_product(
            [range(nts), ["q", "v", "d"]]
        ).to_flat_index()
        if run_parameters.get("return_courant", False):

            flowveldepth = pd.concat(
                [pd.DataFrame(d, index=i, columns=qvd_columns) for i, d, c in results],
                copy=False,
            )
        else:
            flowveldepth = pd.concat(
                [pd.DataFrame(d, index=i, columns=qvd_columns) for i, d in results],
                copy=False,
            )
        # Output new restart CSV
        restart_flows = flowveldepth.iloc[:,-3:]
        restart_flows.index.name = 'link'
        restart_flows.columns = ['qu0', 'qd0', 'h0']
        output_iteration = restart_parameters["wrf_hydro_channel_restart_file"][:-1] + str(ts_iterator+2) + ".csv"
        restart_flows.to_csv(output_iteration)


        if run_parameters.get("return_courant", False):
            courant_columns = pd.MultiIndex.from_product(
                [range(nts), ["cn", "ck", "X"]]
            ).to_flat_index()
            courant = pd.concat(
                [
                    pd.DataFrame(c, index=i, columns=courant_columns)
                    for i, d, c in results
                ],
                copy=False,
            )

        if csv_output_folder:
            # create filenames
            # TO DO: create more descriptive filenames
            if supernetwork_parameters.get("title_string", None):
                filename_fvd = (
                    "flowveldepth_" + supernetwork_parameters["title_string"] + ".csv"
                )
                filename_courant = (
                    "courant_" + supernetwork_parameters["title_string"] + ".csv"
                )
            else:
                run_time_stamp = datetime.now().isoformat()
                filename_fvd = "flowveldepth_" + run_time_stamp + ".csv"
                filename_courant = "courant_" + run_time_stamp + ".csv"

            output_path = pathlib.Path(csv_output_folder).resolve()
            flowveldepth = flowveldepth.sort_index()
            flowveldepth.to_csv(output_path.joinpath(filename_fvd))

            if run_parameters.get("return_courant", False):
                courant = courant.sort_index()
                courant.to_csv(output_path.joinpath(filename_courant))

        if debuglevel <= -1:
            print(flowveldepth)

    if verbose:
        print("output complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))
    

    ################### Parity Check

    if (
        "parity_check_input_folder" in parity_parameters
        or "parity_check_file" in parity_parameters
    ):

        if verbose:
            print(
                "conducting parity check, comparing WRF Hydro results against t-route results"
            )
        if showtiming:
            start_time = time.time()
        build_tests.parity_check(
            parity_parameters,
            run_parameters,
            ts_iterator,
            file_run_size,
            supernetwork_parameters,
            run_parameters["nts"],
            run_parameters["dt"],
            results,
        )

        if verbose:
            print("parity check complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))
    
    if ts_iterator == (runs_to_be_completed-1):
        if verbose:
            print("process complete")
        if showtiming:
            print("%s seconds." % (time.time() - main_start_time))
    else:
        ts_iterator = (ts_iterator + 1)
        restart_file_number = (restart_file_number + 1) 
        main()




if __name__ == "__main__":
    main()
