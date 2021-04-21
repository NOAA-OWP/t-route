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
import sys
import time
from datetime import datetime
import numpy as np
import argparse
from pathlib import Path
import glob
import pandas as pd
from collections import defaultdict
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
        "--compute-kernel",
        nargs="?",
        help="Use the cython version of the compute_network code [options: 'V02-caching'; 'V02-structured'; 'V02-structured-obj' ... ).",
        dest="compute_kernel",
        default="VO2-caching",
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

    parser.add_argument(
        "--data_assimilation_folder_path",
        help="Provide a path to a folder containing the usgs time slice files",
        dest="data_assimilation_parameters_folder",
        default=None,
    )
    parser.add_argument(
        "--data_assimilation_filter",
        help="Provide a glob pattern filter for ncdf files (e.g., 2020-03-21*.usgsTimeSlice.ncdf)",
        dest="data_assimilation_filter",
        default=None,
    )
    parser.add_argument(
        "--data_assimilation_csv",
        help="Provide a csv with the timeslices prepared for use",
        dest="data_assimilation_csv",
        default=None,
    )

    return parser.parse_args()


ENV_IS_CL = False
if ENV_IS_CL:
    root = Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = Path("../..").resolve()
    # sys.path.append(r"../python_framework_v02")

    # TODO: automate compile for the package scripts
    sys.path.append("fast_reach")

## network and reach utilities
import troute.nhd_network_utilities_v02 as nnu
import fast_reach
import diffusive
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
    wbodies,
    reaches_bytw,
    compute_func,
    parallel_compute_method,
    subnetwork_target_size,
    cpu_pool,
    dt,
    nts,
    qts_subdivisions,
    independent_networks,
    param_df,
    q0,
    qlats,
    usgs_df,
    last_obs_df,
    assume_short_ts,
    return_courant,
    waterbodies_df,
    diffusive_parameters=None,
    waterbody_parameters,
):

    param_df["dt"] = dt
    param_df = param_df.astype("float32")

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
                        segs,
                        ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
                    if order < max(subnetworks_only_ordered_jit.keys()):
                        for us_subn_tw in offnetwork_upstreams:
                            subn_tw_sortposition = param_df_sub.index.get_loc(
                                us_subn_tw
                            )
                            flowveldepth_interorder[us_subn_tw][
                                "position_index"
                            ] = subn_tw_sortposition

                    subn_reach_list = clustered_subns["subn_reach_list"]
                    upstreams = clustered_subns["upstreams"]

                    if not usgs_df.empty:
                        usgs_segs = list(usgs_df.index.intersection(param_df_sub.index))
                        nudging_positions_list = param_df_sub.index.get_indexer(
                            usgs_segs
                        )
                        usgs_df_sub = usgs_df.loc[usgs_segs]
                        usgs_df_sub.drop(
                            usgs_df_sub.columns[range(0, 1)], axis=1, inplace=True
                        )
                    else:
                        usgs_df_sub = pd.DataFrame()
                        nudging_positions_list = []

                    last_obs_sub = pd.DataFrame()

                    qlat_sub = qlats.loc[param_df_sub.index]
                    q0_sub = q0.loc[param_df_sub.index]

                    # TODO: Wire in the proper reservoir distinction
                    # At present, in by-subnetwork-jit/jit-clustered, these next two lines
                    # only produce a dummy list, but...
                    # Eventually, the wiring for reservoir simulation needs to be added.
                    subn_reach_type_list = [0 for reaches in subn_reach_list]
                    subn_reach_list_with_type = list(
                        zip(subn_reach_list, subn_reach_type_list)
                    )

                    # results_subn[order].append(
                    #     compute_func(
                    jobs.append(
                        delayed(compute_func)(
                            nts,
                            qts_subdivisions,
                            subn_reach_list_with_type,
                            upstreams,
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            q0_sub.values.astype("float32"),
                            qlat_sub.values.astype("float32"),
                            [],  # lake_segs
                            np.empty(
                                shape=(0, 0), dtype="float64"
                            ),  # waterbodies_df_sub.values
                            usgs_df_sub.values.astype("float32"),
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
                            np.array(nudging_positions_list, dtype="int32"),
                            last_obs_sub.values.astype("float32"),
                            {
                                us: fvd
                                for us, fvd in flowveldepth_interorder.items()
                                if us in offnetwork_upstreams
                            },
                            assume_short_ts,
                            return_courant,
                            diffusive_parameters,
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
                        segs,
                        ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
                    if order < max(subnetworks_only_ordered_jit.keys()):
                        for us_subn_tw in offnetwork_upstreams:
                            subn_tw_sortposition = param_df_sub.index.get_loc(
                                us_subn_tw
                            )
                            flowveldepth_interorder[us_subn_tw][
                                "position_index"
                            ] = subn_tw_sortposition

                    if not usgs_df.empty:
                        usgs_segs = list(usgs_df.index.intersection(param_df_sub.index))
                        nudging_positions_list = param_df_sub.index.get_indexer(
                            usgs_segs
                        )
                        usgs_df_sub = usgs_df.loc[usgs_segs]
                        usgs_df_sub.drop(
                            usgs_df_sub.columns[range(0, 1)], axis=1, inplace=True
                        )
                    else:
                        usgs_df_sub = pd.DataFrame()
                        nudging_positions_list = []

                    last_obs_sub = pd.DataFrame()

                    qlat_sub = qlats.loc[param_df_sub.index]
                    q0_sub = q0.loc[param_df_sub.index]

                    # At present, in by-subnetwork-jit/jit-clustered, these next two lines
                    # only produce a dummy list, but...
                    # Eventually, the wiring for reservoir simulation needs to be added.
                    subn_reach_type_list = [0 for reaches in subn_reach_list]
                    subn_reach_list_with_type = list(
                        zip(subn_reach_list, subn_reach_type_list)
                    )

                    jobs.append(
                        delayed(compute_func)(
                            nts,
                            qts_subdivisions,
                            subn_reach_list_with_type,
                            subnetworks[subn_tw],
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            q0_sub.values.astype("float32"),
                            qlat_sub.values.astype("float32"),
                            [],  # lake_segs
                            np.empty(
                                shape=(0, 0), dtype="float64"
                            ),  # waterbodies_df_sub.values
                            usgs_df_sub.values.astype("float32"),
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
                            np.array(nudging_positions_list, dtype="int32"),
                            {
                                us: fvd
                                for us, fvd in flowveldepth_interorder.items()
                                if us in offnetwork_upstreams
                            },
                            assume_short_ts,
                            return_courant,
                            diffusive_parameters,
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
                # The X_sub lines use SEGS...
                # which is now invalid with the wbodies included.
                # So we define "common_segs" to identify regular routing segments
                # and wbodies_segs for the waterbody reaches/segments
                segs = list(chain.from_iterable(reach_list))
                common_segs = param_df.index.intersection(segs)
                # Assumes everything else is a waterbody...
                wbodies_segs = set(segs).symmetric_difference(common_segs)

                # If waterbody parameters exist
                if not waterbodies_df.empty:

                    lake_segs = list(waterbodies_df.index.intersection(segs))

                    waterbodies_df_sub = waterbodies_df.loc[
                        lake_segs,
                        [
                            "LkArea",
                            "LkMxE",
                            "OrificeA",
                            "OrificeC",
                            "OrificeE",
                            "WeirC",
                            "WeirE",
                            "WeirL",
                            "ifd",
                            "qd0",
                            "h0",
                        ],
                    ]

                else:
                    lake_segs = []
                    waterbodies_df_sub = pd.DataFrame()

                param_df_sub = param_df.loc[
                    common_segs,
                    ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                ].sort_index()

                if not usgs_df.empty:
                    usgs_segs = list(usgs_df.index.intersection(param_df_sub.index))
                    nudging_positions_list = param_df_sub.index.get_indexer(usgs_segs)
                    usgs_df_sub = usgs_df.loc[usgs_segs]
                    usgs_df_sub.drop(
                        usgs_df_sub.columns[range(0, 1)], axis=1, inplace=True
                    )
                else:
                    usgs_df_sub = pd.DataFrame()
                    nudging_positions_list = []

                last_obs_sub = pd.DataFrame()

                reaches_list_with_type = []

                for reaches in reach_list:
                    if set(reaches) & wbodies_segs:
                        reach_type = 1  # type 1 for waterbody/lake
                    else:
                        reach_type = 0  # type 0 for reach

                    reach_and_type_tuple = (reaches, reach_type)

                    reaches_list_with_type.append(reach_and_type_tuple)

                # qlat_sub = qlats.loc[common_segs].sort_index()
                # q0_sub = q0.loc[common_segs].sort_index()
                qlat_sub = qlats.loc[param_df_sub.index]
                q0_sub = q0.loc[param_df_sub.index]

                param_df_sub = param_df_sub.reindex(
                    param_df_sub.index.tolist() + lake_segs
                ).sort_index()
                qlat_sub = qlat_sub.reindex(param_df_sub.index)
                q0_sub = q0_sub.reindex(param_df_sub.index)

                jobs.append(
                    delayed(compute_func)(
                        nts,
                        qts_subdivisions,
                        reaches_list_with_type,
                        independent_networks[tw],
                        param_df_sub.index.values.astype("int64"),
                        param_df_sub.columns.values,
                        param_df_sub.values,
                        q0_sub.values.astype("float32"),
                        qlat_sub.values.astype("float32"),
                        lake_segs,
                        waterbodies_df_sub.values,
                        waterbody_parameters,
                        usgs_df_sub.values.astype("float32"),
                        np.array(nudging_positions_list, dtype="int32"),
                        last_obs_sub.values.astype("float32"),
                        {},
                        assume_short_ts,
                        return_courant,
                        diffusive_parameters,
                    )
                )
            results = parallel(jobs)

    else:  # Execute in serial
        results = []
        for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
            # The X_sub lines use SEGS...
            # which becomes invalid with the wbodies included.
            # So we define "common_segs" to identify regular routing segments
            # and wbodies_segs for the waterbody reaches/segments
            segs = list(chain.from_iterable(reach_list))
            common_segs = param_df.index.intersection(segs)
            # Assumes everything else is a waterbody...
            wbodies_segs = set(segs).symmetric_difference(common_segs)

            # If waterbody parameters exist
            if not waterbodies_df.empty:

                lake_segs = list(waterbodies_df.index.intersection(segs))

                waterbodies_df_sub = waterbodies_df.loc[
                    lake_segs,
                    [
                        "LkArea",
                        "LkMxE",
                        "OrificeA",
                        "OrificeC",
                        "OrificeE",
                        "WeirC",
                        "WeirE",
                        "WeirL",
                        "ifd",
                        "qd0",
                        "h0",
                    ],
                ]

            else:
                lake_segs = []
                waterbodies_df_sub = pd.DataFrame()

            param_df_sub = param_df.loc[
                common_segs,
                ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
            ].sort_index()

            if not usgs_df.empty:
                usgs_segs = list(usgs_df.index.intersection(param_df_sub.index))
                nudging_positions_list = param_df_sub.index.get_indexer(usgs_segs)
                usgs_df_sub = usgs_df.loc[usgs_segs]
                usgs_df_sub.drop(usgs_df_sub.columns[range(0, 1)], axis=1, inplace=True)
            else:
                usgs_df_sub = pd.DataFrame()
                nudging_positions_list = []

            if not last_obs_df.empty:
                pass
            #     lastobs_segs = list(last_obs_df.index.intersection(param_df_sub.index))
            #     nudging_positions_list = param_df_sub.index.get_indexer(lastobs_segs)
            #     last_obs_sub = last_obs_df.loc[lastobs_segs]
            else:
                last_obs_sub = pd.DataFrame()
            #     nudging_positions_list = []

            # qlat_sub = qlats.loc[common_segs].sort_index()
            # q0_sub = q0.loc[common_segs].sort_index()
            qlat_sub = qlats.loc[param_df_sub.index]
            q0_sub = q0.loc[param_df_sub.index]

            param_df_sub = param_df_sub.reindex(
                param_df_sub.index.tolist() + lake_segs
            ).sort_index()
            qlat_sub = qlat_sub.reindex(param_df_sub.index)
            q0_sub = q0_sub.reindex(param_df_sub.index)

            reach_type_list = [
                1 if (set(reaches) & wbodies_segs) else 0 for reaches in reach_list
            ]
            reaches_list_with_type = list(zip(reach_list, reach_type_list))
            """
            reaches_list_with_type = []

            for reaches in reach_list:
                if (set(reaches) & wbodies_segs):
                    reach_type = 1 # type 1 for waterbody/lake
                else:
                    reach_type = 0 # type 0 for reach

                reach_and_type_tuple = (reaches, reach_type)

                reaches_list_with_type.append(reach_and_type_tuple)
            """

            results.append(
                compute_func(
                    nts,
                    qts_subdivisions,
                    reaches_list_with_type,
                    independent_networks[tw],
                    param_df_sub.index.values.astype("int64"),
                    param_df_sub.columns.values,
                    param_df_sub.values,
                    q0_sub.values.astype("float32"),
                    qlat_sub.values.astype("float32"),
                    lake_segs,
                    waterbodies_df_sub.values,
                    waterbody_parameters,
                    usgs_df_sub.values.astype("float32"),
                    np.array(nudging_positions_list, dtype="int32"),
                    last_obs_sub.values.astype("float32"),
                    {},
                    assume_short_ts,
                    return_courant,
                    diffusive_parameters,
                )
            )

    return results


def _input_handler():

    args = _handle_args()

    custom_input_file = args.custom_input_file
    log_parameters = {}
    supernetwork_parameters = None
    waterbody_parameters = {}
    compute_parameters = {}
    forcing_parameters = {}
    restart_parameters = {}
    output_parameters = {}
    parity_parameters = {}
    data_assimilation_parameters = {}
    diffusive_parameters = {}
    coastal_parameters = {}

    if custom_input_file:
        (
            log_parameters,
            supernetwork_parameters,
            waterbody_parameters,
            compute_parameters,
            forcing_parameters,
            restart_parameters,
            diffusive_parameters,
            coastal_parameters,
            output_parameters,
            parity_parameters,
            data_assimilation_parameters,
        ) = nhd_io.read_custom_input_new(custom_input_file)

    else:
        # TODO: Fix CLI input path
        compute_parameters["assume_short_ts"] = args.assume_short_ts
        compute_parameters["return_courant"] = args.return_courant
        compute_parameters["parallel_compute_method"] = args.parallel_compute_method
        compute_parameters["subnetwork_target_size"] = args.subnetwork_target_size
        compute_parameters["cpu_pool"] = args.cpu_pool
        compute_parameters["compute_kernel"] = args.compute_kernel

        log_parameters["showtiming"] = args.showtiming
        log_parameters["debuglevel"] = debuglevel = -1 * args.debuglevel
        log_parameters["verbose"] = verbose = args.verbose

        output_parameters["csv_output"] = {}
        output_parameters["csv_output"]["csv_output_folder"] = args.csv_output_folder

        test_folder = Path(root, "test")
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
            compute_parameters["dt"] = args.dt
            compute_parameters["nts"] = args.nts
            compute_parameters["qts_subdivisions"] = args.qts_subdivisions
            compute_parameters["compute_kernel"] = args.compute_kernel

            waterbody_parameters[
                "break_network_at_waterbodies"
            ] = args.break_network_at_waterbodies

            data_assimilation_parameters[
                "data_assimilation_parameters_folder"
            ] = args.data_assimilation_parameters_folder
            data_assimilation_parameters[
                "data_assimilation_filter"
            ] = args.data_assimilation_filter
            data_assimilation_parameters[
                "data_assimilation_csv"
            ] = args.data_assimilation_csv

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
        log_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        diffusive_parameters,
        coastal_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    )


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
    connections, param_df, wbodies, gages = nnu.build_connections(
        supernetwork_parameters,
    )

    break_network_at_waterbodies = waterbody_parameters.get(
        "break_network_at_waterbodies", False
    )
    break_network_at_gages = supernetwork_parameters.get(
        "break_network_at_gages", False
    )

    if (
        not wbodies
    ):  # Turn off any further reservoir processing if the network contains no waterbodies
        break_network_at_waterbodies = False

    if break_network_at_waterbodies:
        connections = nhd_network.replace_waterbodies_connections(connections, wbodies)

    if verbose:
        print("supernetwork connections set complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    ################################
    ## STEP 3a: Read waterbody parameter file
    # waterbodies_values = supernetwork_values[12]
    # waterbodies_segments = supernetwork_values[13]
    # connections_tailwaters = supernetwork_values[4]

    if break_network_at_waterbodies:
        # Read waterbody parameters
        waterbodies_df = nhd_io.read_waterbody_df(
            waterbody_parameters, {"level_pool": wbodies.values()}
        )

        # Remove duplicate lake_ids and rows
        waterbodies_df = (
            waterbodies_df.reset_index()
            .drop_duplicates(subset="lake_id")
            .set_index("lake_id")
        )

        #Check if hybrid-usgs, hybrid-usace, or rfc type reservoirs are set to true
        wbtype="level_pool"
        wb_params = waterbody_parameters[wbtype]

        if wb_params["reservoir_persistence_usgs"] or wb_params["reservoir_persistence_usace"] \
        or wb_params["reservoir_rfc_forecasts"]:

            #reservoir_types_df = nhd_io.read_level_pool_waterbody_df(wb_params["reservoir_parameter_file"], \
            #    wb_params["level_pool_waterbody_id"], {"level_pool": wbodies.values()},) 
            
            #reservoir_types_df = nhd_io.read_level_pool_waterbody_df(wb_params["reservoir_parameter_file"], \
            #    wb_params["level_pool_waterbody_id"],) 

            print ("wb_params[level_pool_waterbody_i]")
            print (wb_params["level_pool_waterbody_id"])
            print ({"level_pool": wbodies.values()})
            print ("============")
            print (wb_params["reservoir_parameter_file"])

            #reservoir_types_df = nhd_io.read_reservoir_parameter_file(wb_params["reservoir_parameter_file"], \
            #    wb_params["level_pool_waterbody_id"], {"level_pool": wbodies.values()},) 
            reservoir_types_df = nhd_io.read_reservoir_parameter_file(wb_params["reservoir_parameter_file"], \
                wb_params["level_pool_waterbody_id"], wbodies.values(),) 

            #reservoir_types_df = nhd_io.read_reservoir_parameter_file(wb_params["reservoir_parameter_file"])

            print ("=======Reading hybrid----------------")
            print (reservoir_types_df)
            print ("============")


        else:
            print ("========hybrid off -------------------------------")



    else:
        waterbodies_df = pd.DataFrame()

    # STEP 2: Identify Independent Networks and Reaches by Network
    if showtiming:
        start_time = time.time()
    if verbose:
        print("organizing connections into reaches ...")

    network_break_segments = set()
    if break_network_at_waterbodies:
        network_break_segments = network_break_segments.union(wbodies.values())
    if break_network_at_gages:
        network_break_segments = network_break_segments.union(gages.keys())

    independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections, network_break_segments,
    )
    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    return (
        connections,
        param_df,
        wbodies,
        waterbodies_df,
        break_network_at_waterbodies,
        independent_networks,
        reaches_bytw,
        rconn,
    )


def nwm_initial_warmstate_preprocess(
    break_network_at_waterbodies,
    restart_parameters,
    param_df,
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
            waterbodies_initial_depth_const = -1.0
            # Set initial states from cold-state
            waterbodies_initial_states_df = pd.DataFrame(
                0, index=waterbodies_df.index, columns=["qd0", "h0",], dtype="float32"
            )
            # TODO: This assignment could probably by done in the above call
            waterbodies_initial_states_df["qd0"] = waterbodies_initial_ds_flow_const
            waterbodies_initial_states_df["h0"] = waterbodies_initial_depth_const
            waterbodies_initial_states_df["index"] = range(
                len(waterbodies_initial_states_df)
            )

        if verbose:
            print("waterbody initial states complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))
            start_time = time.time()

    # STEP 4: Handle Channel Initial States
    if showtiming:
        start_time = time.time()
    if verbose:
        print("setting channel initial states ...")

    q0 = nnu.build_channel_initial_state(
        restart_parameters, supernetwork_parameters, param_df.index
    )

    if verbose:
        print("channel initial states complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))
        start_time = time.time()

    waterbodies_df = pd.merge(
        waterbodies_df, waterbodies_initial_states_df, on="lake_id"
    )

    last_obs_file = data_assimilation_parameters.get("wrf_hydro_last_obs_file", None)
    last_obs_df = pd.DataFrame()

    return waterbodies_df, q0, last_obs_df
    # TODO: This returns a full dataframe (waterbodies_df) with the
    # merged initial states for waterbodies, but only the
    # initial state values (q0; not merged with the channel properties)
    # for the channels --
    # That is because that is how they are used downstream. Need to
    # trace that back and decide if there is one of those two ways
    # that is optimal and make both returns that way.


def nwm_forcing_preprocess(
    connections,
    break_network_at_waterbodies,
    run,
    forcing_parameters,
    data_assimilation_parameters,
    showtiming=False,
    verbose=False,
    debuglevel=0,
):

    nts = forcing_parameters.get("nts", None)
    dt = forcing_parameters.get("dt", None)
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", None)
    qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
    qlat_file_index_col = forcing_parameters.get("qlat_file_index_col", None)
    qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", None)

    # TODO: find a better way to deal with these defaults and overrides.
    run["nts"] = run.get("nts", nts)
    run["dt"] = run.get("dt", dt)
    run["qts_subdivisions"] = run.get("qts_subdivisions", qts_subdivisions)
    run["qlat_input_folder"] = run.get("qlat_input_folder", qlat_input_folder)
    run["qlat_file_index_col"] = run.get("qlat_file_index_col", qlat_file_index_col)
    run["qlat_file_value_col"] = run.get("qlat_file_value_col", qlat_file_value_col)

    # STEP 5: Read (or set) QLateral Inputs
    if showtiming:
        start_time = time.time()
    if verbose:
        print("creating qlateral array ...")

    qlats_df = nnu.build_qlateral_array(
        run, connections.keys(), supernetwork_parameters,
    )

    if verbose:
        print("qlateral array complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 6
    data_assimilation_csv = data_assimilation_parameters.get(
        "data_assimilation_csv", None
    )
    data_assimilation_folder = data_assimilation_parameters.get(
        "data_assimilation_timeslices_folder", None
    )
    if data_assimilation_csv or data_assimilation_folder:

        if data_assimilation_folder and data_assimilation_csv:
            print(
                "Please select data_assimilation_parameters_folder + data_assimilation_filter or data_assimilation_csv not both."
            )

        if showtiming:
            start_time = time.time()
        if verbose:
            print("creating usgs time_slice data array ...")

        usgs_df, _ = nnu.build_data_assimilation(data_assimilation_parameters)

        if verbose:
            print("usgs array complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))

    else:
        usgs_df = pd.DataFrame()

    # STEP 7
    coastal_boundary_elev = coastal_parameters.get("coastal_boundary_elev_data", None)
    coastal_ncdf = coastal_parameters.get("coastal_ncdf", None)

    if coastal_boundary_elev:
        print("creating coastal dataframe ...")
        coastal_df = nhd_io.build_coastal_dataframe(coastal_boundary_elev)

    if coastal_ncdf:
        print("creating coastal ncdf dataframe ...")
        coastal_ncdf_df = nhd_io.build_coastal_ncdf_dataframe(coastal_ncdf)

    # TODO: disentangle the implicit (run) and explicit (qlats_df, usgs_df) returns
    return qlats_df, usgs_df


def nwm_route(
    downstream_connections,
    upstream_connections,
    waterbodies_in_connections,
    reaches_bytw,
    parallel_compute_method,
    subnetwork_target_size,
    cpu_pool,
    dt,
    nts,
    qts_subdivisions,
    independent_networks,
    param_df,
    q0,
    qlats,
    usgs_df,
    last_obs_df,
    assume_short_ts,
    return_courant,
    waterbodies_df,
    diffusive_parameters,
    showtiming=False,
    verbose=False,
    debuglevel=0,
):

    ################### Main Execution Loop across ordered networks
    if showtiming:
        start_time = time.time()
    if verbose:
        if compute_parameters.get("return_courant", False):
            print(
                f"executing routing computation, with Courant evaluation metrics returned"
            )
        else:
            print(f"executing routing computation ...")

    if compute_parameters.get("compute_kernel", None) == "diffusive":
        compute_func = diffusive.compute_diffusive_tst
    elif compute_parameters.get("compute_kernel", None) == "V02-caching":
        compute_func = fast_reach.compute_network
    elif compute_parameters.get("compute_kernel", None) == "V02-structured":
        compute_func = fast_reach.compute_network_structured
    elif compute_parameters.get("compute_kernel", None) == "V02-structured-obj":
        compute_func = fast_reach.compute_network_structured_obj
    else:
        compute_func = fast_reach.compute_network

    # TODO: Remove below. --compute-kernel=V02-structured-obj did not work on command line
    # compute_func = fast_reach.compute_network_structured_obj

    results = compute_nhd_routing_v02(
        connections,
        rconn,
        wbodies,
        reaches_bytw,
        compute_func,
        parallel_compute_method,
        subnetwork_target_size,  # The default here might be the whole network or some percentage...
        cpu_pool,
        dt,
        nts,
        qts_subdivisions,
        independent_networks,
        param_df,
        q0,
        qlats,
        usgs_df,
        last_obs_df,
        compute_parameters.get("assume_short_ts", False),
        compute_parameters.get("return_courant", False),
        waterbodies_df,
        diffusive_parameters,
        waterbody_parameters,
    )

    with open("mainstems_conus.txt", "w") as filehandle:
        for listitem in reaches_bytw.keys():
            filehandle.write("%s\n" % listitem)

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    return results


def nwm_output_generator(
    results,
    output_parameters,
    parity_parameters,
    parity_set,
    return_courant,
    showtiming=False,
    verbose=False,
    debuglevel=0,
):

    parity_check_input_folder = parity_parameters.get("parity_check_input_folder", None)
    parity_check_file_index_col = parity_parameters.get(
        "parity_check_file_index_col", None
    )
    parity_check_file_value_col = parity_parameters.get(
        "parity_check_file_value_col", None
    )
    parity_check_compare_node = parity_parameters.get("parity_check_compare_node", None)

    # TODO: find a better way to deal with these defaults and overrides.
    parity_set["parity_check_input_folder"] = parity_set.get(
        "parity_check_input_folder", parity_check_input_folder
    )
    parity_set["parity_check_file_index_col"] = parity_set.get(
        "parity_check_file_index_col", parity_check_file_index_col
    )
    parity_set["parity_check_file_value_col"] = parity_set.get(
        "parity_check_file_value_col", parity_check_file_value_col
    )
    parity_set["parity_check_compare_node"] = parity_set.get(
        "parity_check_compare_node", parity_check_compare_node
    )

    ################### Output Handling
    if showtiming:
        start_time = time.time()
    if verbose:
        print(f"Handling output ...")

    csv_output = output_parameters.get("csv_output", None)
    if csv_output:
        csv_output_folder = output_parameters["csv_output"].get(
            "csv_output_folder", None
        )
        csv_output_segments = csv_output.get("csv_output_segments", None)

    if (debuglevel <= -1) or csv_output:

        qvd_columns = pd.MultiIndex.from_product(
            [range(nts), ["q", "v", "d"]]
        ).to_flat_index()

        flowveldepth = pd.concat(
            [pd.DataFrame(r[1], index=r[0], columns=qvd_columns) for r in results],
            copy=False,
        )

        if return_courant:
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

        # directory containing WRF Hydro restart files
        wrf_hydro_restart_dir = output_parameters.get(
            "wrf_hydro_channel_restart_directory", None
        )
        if wrf_hydro_restart_dir:

            wrf_hydro_channel_restart_new_extension = output_parameters.get(
                "wrf_hydro_channel_restart_new_extension", "TRTE"
            )

            # list of WRF Hydro restart files
            wrf_hydro_restart_files = list(
                Path(wrf_hydro_restart_dir).glob(
                    output_parameters["wrf_hydro_channel_restart_pattern_filter"]
                    + "[!"
                    + wrf_hydro_channel_restart_new_extension
                    + "]"
                )
            )

            if len(wrf_hydro_restart_files) > 0:
                nhd_io.write_channel_restart_to_wrf_hydro(
                    flowveldepth,
                    wrf_hydro_restart_files,
                    restart_parameters.get("wrf_hydro_channel_restart_file"),
                    run_parameters.get("dt"),
                    run_parameters.get("nts"),
                    restart_parameters.get("wrf_hydro_channel_ID_crosswalk_file"),
                    restart_parameters.get(
                        "wrf_hydro_channel_ID_crosswalk_file_field_name"
                    ),
                    wrf_hydro_channel_restart_new_extension,
                )
            else:
                # print error and/or raise exception
                print(
                    "WRF Hydro restart files not found - Aborting restart write sequence"
                )
                raise AssertionError

        chrtout_folder = output_parameters.get("wrf_hydro_channel_output_folder", None)
        if chrtout_folder:
            wrf_hydro_channel_output_new_extension = output_parameters.get(
                "wrf_hydro_channel_output_new_extension", "TRTE"
            )
            chrtout_files = list(
                Path(chrtout_folder).glob(
                    output_parameters["wrf_hydro_channel_output_file_pattern_filter"]
                )
            )
            nhd_io.write_q_to_wrf_hydro(
                flowveldepth, chrtout_files, run_parameters["qts_subdivisions"]
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

            output_path = Path(csv_output_folder).resolve()

            flowveldepth = flowveldepth.sort_index()
            flowveldepth.to_csv(output_path.joinpath(filename_fvd))

            if return_courant:
                courant = courant.sort_index()
                courant.to_csv(output_path.joinpath(filename_courant))

            usgs_df_filtered = usgs_df[usgs_df.index.isin(csv_output_segments)]
            usgs_df_filtered.to_csv(output_path.joinpath("usgs_df.csv"))

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
        or "parity_check_waterbody_file" in parity_parameters
    ):

        if verbose:
            print(
                "conducting parity check, comparing WRF Hydro results against t-route results"
            )
        if showtiming:
            start_time = time.time()

        build_tests.parity_check(
            parity_set, run_results,
        )

        if verbose:
            print("parity check complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))


def new_nwm_q0(run_results):
    return pd.concat(
        # TODO: we only need two fields, technically, and the restart file produced by WRF-Hydro
        # actually contains a field qu0, which is never used for restart (the qu0 can be obtained
        # as the qd0 from the topologically upstream segments, just like during the calculation).
        # In any case, the qu0 currently in the WRF-Hydro output is populated with the same value
        # as the qd0.
        # [pd.DataFrame(d[:,-3::2], index=i, columns=["qd0", "h0"]) for i, d in run_results],
        # [pd.DataFrame(r[1][:,-3:], index=r[0], columns=["qu0", "v0", "h0"]) for r in run_results],
        [
            pd.DataFrame(
                r[1][:, [-3, -3, -1]], index=r[0], columns=["qu0", "qd0", "h0"]
            )
            for r in run_results
        ],
        copy=False,
    )


if __name__ == "__main__":

    (
        log_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        diffusive_parameters,
        coastal_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    ) = _input_handler()

    verbose = log_parameters.get("verbose", None)
    showtiming = log_parameters.get("showtiming", None)
    debuglevel = log_parameters.get("debuglevel", 0)

    if showtiming:
        main_start_time = time.time()

    (
        connections,
        param_df,
        wbodies,
        waterbodies_df,
        break_network_at_waterbodies,
        independent_networks,
        reaches_bytw,
        rconn,
    ) = nwm_network_preprocess(
        supernetwork_parameters,
        waterbody_parameters,
        showtiming=showtiming,
        verbose=verbose,
        debuglevel=debuglevel,
    )

    waterbodies_df, q0, last_obs_df = nwm_initial_warmstate_preprocess(
        break_network_at_waterbodies,
        restart_parameters,
        param_df,
        waterbodies_df,
        segment_list=None,
        wbodies_list=None,
        showtiming=showtiming,
        verbose=verbose,
        debuglevel=debuglevel,
    )

    # The inputs below assume a very pedantic setup
    # with each run set explicitly defined, so...
    # TODO: Make this more flexible.
    run_sets = forcing_parameters.get("qlat_forcing_sets", False)

    # TODO: Data Assimilation will be something like the parity block
    # if DA:
    #     da_sets = [BIG LIST OF DA BLOCKS]

    if "wrf_hydro_parity_check" in output_parameters:
        parity_sets = parity_parameters.get("parity_check_compare_file_sets", False)
    else:
        parity_sets = []

    parallel_compute_method = (compute_parameters.get("parallel_compute_method", None),)
    subnetwork_target_size = (compute_parameters.get("subnetwork_target_size", 1),)
    cpu_pool = (compute_parameters.get("cpu_pool", None),)
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", 1)

    qlats, usgs_df = nwm_forcing_preprocess(
        connections,
        break_network_at_waterbodies,
        run_sets[0],
        forcing_parameters,
        data_assimilation_parameters,
        showtiming,
        verbose,
        debuglevel,
    )

    for run_set_iterator, run in enumerate(run_sets):

        dt = run.get("dt")
        nts = run.get("nts")
        if parity_sets:
            parity_sets[run_set_iterator]["dt"] = dt
            parity_sets[run_set_iterator]["nts"] = nts

        run_results = nwm_route(
            connections,
            rconn,
            wbodies,
            reaches_bytw,
            parallel_compute_method,
            subnetwork_target_size,
            cpu_pool,
            dt,
            nts,
            qts_subdivisions,
            independent_networks,
            param_df,
            q0,
            qlats,
            usgs_df,
            last_obs_df,
            compute_parameters.get("assume_short_ts", False),
            compute_parameters.get("return_courant", False),
            waterbodies_df,
            diffusive_parameters,
            showtiming,
            verbose,
            debuglevel,
        )

        if (
            run_set_iterator < len(run_sets) - 1
        ):  # No forcing to prepare for the last loop
            qlats, usgs_df = nwm_forcing_preprocess(
                connections,
                break_network_at_waterbodies,
                run_sets[run_set_iterator + 1],
                forcing_parameters,
                data_assimilation_parameters,
                showtiming,
                verbose,
                debuglevel,
            )

            # q0 = run_results
            q0 = new_nwm_q0(run_results)

        nwm_output_generator(
            run_results,
            output_parameters,
            parity_parameters,
            parity_sets[run_set_iterator],
            compute_parameters.get("return_courant", False),
            showtiming,
            verbose,
            debuglevel,
        )

    # nwm_final_output_generator()

    if verbose:
        print("process complete")
    if showtiming:
        print("%s seconds." % (time.time() - main_start_time))

        """
        Asynchronous execution Psuedocode
        Sync1: Prepare first warmstate from files
        Sync1: Prepare first forcing from files

        For first forcing set
            Sync2a: run model
            Sync2b: begin preparing next forcing
            Sync3a - AFTER Sync2a, prepare next warmstate (last state of model run)
            Sync3b: write any output from Sync2a
            Loop has to wait for Sync2a+b+Sync3a, does not have to wait for Sync3b
                  if next forcing prepared
        """
