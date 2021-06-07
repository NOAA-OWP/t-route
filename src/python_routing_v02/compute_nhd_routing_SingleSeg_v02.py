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
        "--compute-method",
        nargs="?",
        help="Use the cython version of the compute_network code [options: 'V02-caching'; 'V02-structured'; 'V02-structured-obj' ... ).",
        dest="compute_method",
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
    nts,
    qts_subdivisions,
    independent_networks,
    param_df,
    q0,
    qlats,
    usgs_df,
    assume_short_ts,
    return_courant,
    waterbodies_df,
    diffusive_parameters=None,
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

                    qlat_sub = qlats.loc[param_df_sub.index]
                    q0_sub = q0.loc[param_df_sub.index]

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
        
        # Construct a dictionary of subnetwork segments sorted by subnetwork order
        # subnetworks_only_ordered_jit (dict) 
        #      (subnetwork_order (int): (subnetwork_tailwater (int): subnetwork_segments (set))) 
        subnetworks_only_ordered_jit = defaultdict(dict)
        subnetworks = defaultdict(dict)
        for tw, ordered_network in networks_with_subnetworks_ordered_jit.items():
            intw = independent_networks[tw]
            for order, subnet_sets in ordered_network.items():
                subnetworks_only_ordered_jit[order].update(subnet_sets)
                for subn_tw, subnetwork in subnet_sets.items():
                    subnetworks[subn_tw] = {k: intw[k] for k in subnetwork}

        # Construct a dictionary of subnetwork reaches sorted by subnetwork order
        # reaches_ordered_bysubntw (dict) 
        #      (subnetwork_order (int): (subnetwork_tailwater (int): subnetwork_reaches (list of lists))) 
        reaches_ordered_bysubntw = defaultdict(dict)
        for order, ordered_subn_dict in subnetworks_only_ordered_jit.items():
            for subn_tw, subnet in ordered_subn_dict.items():
                conn_subn = {k: connections[k] for k in subnet if k in connections}
                rconn_subn = {k: rconn[k] for k in subnet if k in rconn}
                
                # break subnetwork reaches at waterbieds if reservoirs are included
                if waterbodies_df.empty:
                    path_func = partial(nhd_network.split_at_junction, rconn_subn)
                else:
                    path_func = partial(
                        nhd_network.split_at_waterbodies_and_junctions,
                        list(waterbodies_df.index.values),
                        rconn_subn
                        )
                    
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
                    
                    # a list of all segments in the subnetwork
                    segs = list(chain.from_iterable(subn_reach_list))
                    
                    # identify upstream neighbors of this subnetwork
                    offnetwork_upstreams = set()
                    segs_set = set(segs)
                    for seg in segs:
                        for us in rconn[seg]:
                            if us not in segs_set:
                                offnetwork_upstreams.add(us)

                    # add upstream neighbors to common_segs list
                    segs.extend(offnetwork_upstreams)
                    
                    # Define "common_segs" to identify regular routing segments
                    common_segs = list(param_df.index.intersection(segs))
                    
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
                        
                    # extract routing/geometry parameters for stream segments in the subnetwork
                    param_df_sub = param_df.loc[
                        common_segs,
                        ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
                    
                    # need this to include water bodies in param_df_sub, otherwise if 
                    # us_subn_tw is water body is will not be found in param_df_sub
                    param_df_sub_super = param_df_sub.reindex(
                        param_df_sub.index.tolist() + lake_segs
                    ).sort_index()
                        
                    if order < max(subnetworks_only_ordered_jit.keys()):
                        for us_subn_tw in offnetwork_upstreams:
                            subn_tw_sortposition = param_df_sub_super.index.get_loc(
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

                    # classify reaches as either waterbodies (1) or rivers (2)
                    subn_reach_list_with_type = []
                    for reaches in subn_reach_list:
                        if set(reaches) & wbodies_segs:
                            reach_type = 1  # type 1 for waterbody/lake
                        else:
                            reach_type = 0  # type 0 for reach

                        reach_and_type_tuple = (reaches, reach_type)

                        subn_reach_list_with_type.append(reach_and_type_tuple)


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
                            subn_reach_list_with_type,
                            subnetworks[subn_tw],
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            q0_sub.values.astype("float32"),
                            qlat_sub.values.astype("float32"),
                            lake_segs,  # lake_segs
                            waterbodies_df_sub.values,
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
                        usgs_df_sub.values.astype("float32"),
                        np.array(nudging_positions_list, dtype="int32"),
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
                    usgs_df_sub.values.astype("float32"),
                    np.array(nudging_positions_list, dtype="int32"),
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
    supernetwork_parameters = None
    waterbody_parameters = {}
    forcing_parameters = {}
    restart_parameters = {}
    output_parameters = {}
    run_parameters = {}
    parity_parameters = {}
    data_assimilation_parameters = {}
    diffusive_parameters = {}
    coastal_parameters = {}

    if custom_input_file:
        (
            supernetwork_parameters,
            waterbody_parameters,
            forcing_parameters,
            restart_parameters,
            output_parameters,
            run_parameters,
            parity_parameters,
            data_assimilation_parameters,
            diffusive_parameters,
            coastal_parameters,
        ) = nhd_io.read_custom_input(custom_input_file)

    else:
        run_parameters["assume_short_ts"] = args.assume_short_ts
        run_parameters["return_courant"] = args.return_courant
        run_parameters["parallel_compute_method"] = args.parallel_compute_method
        run_parameters["subnetwork_target_size"] = args.subnetwork_target_size
        run_parameters["cpu_pool"] = args.cpu_pool
        run_parameters["showtiming"] = args.showtiming
        run_parameters["compute_method"] = args.compute_method

        run_parameters["debuglevel"] = debuglevel = -1 * args.debuglevel
        run_parameters["verbose"] = verbose = args.verbose

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
            run_parameters["dt"] = args.dt
            run_parameters["nts"] = args.nts
            run_parameters["qts_subdivisions"] = args.qts_subdivisions
            run_parameters["compute_method"] = args.compute_method

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
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        parity_parameters,
        data_assimilation_parameters,
        diffusive_parameters,
        coastal_parameters,
    )


def main():

    (
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        parity_parameters,
        data_assimilation_parameters,
        diffusive_parameters,
        coastal_parameters,
    ) = _input_handler()

    dt = run_parameters.get("dt", None)
    nts = run_parameters.get("nts", None)
    verbose = run_parameters.get("verbose", None)
    showtiming = run_parameters.get("showtiming", None)
    debuglevel = run_parameters.get("debuglevel", 0)
    break_network_at_waterbodies = run_parameters.get(
        "break_network_at_waterbodies", False
    )

    if showtiming:
        main_start_time = time.time()

    if verbose:
        print("creating supernetwork connections set")
    if showtiming:
        start_time = time.time()

    # STEP 1: Build basic network connections graph
    connections, param_df = nnu.build_connections(supernetwork_parameters, dt,)
    wbodies = nnu.build_waterbodies(param_df, supernetwork_parameters, "waterbody")
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
        waterbodies_df_reduced = (
            waterbodies_df.reset_index()
            .drop_duplicates(subset="lake_id")
            .set_index("lake_id")
        )
    else:
        waterbodies_df_reduced = pd.DataFrame()

    # STEP 2: Identify Independent Networks and Reaches by Network
    if showtiming:
        start_time = time.time()
    if verbose:
        print("organizing connections into reaches ...")

    independent_networks, reaches_bytw, rconn = nnu.organize_independent_networks(
        connections,
        list(waterbodies_df_reduced.index.values)
        if break_network_at_waterbodies
        else None,
    )
    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

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

        waterbodies_df_reduced = pd.merge(
            waterbodies_df_reduced, waterbodies_initial_states_df, on="lake_id"
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

    # STEP 5: Read (or set) QLateral Inputs
    if showtiming:
        start_time = time.time()
    if verbose:
        print("creating qlateral array ...")

    qlats = nnu.build_qlateral_array(
        forcing_parameters,
        connections.keys(),
        supernetwork_parameters,
        nts,
        run_parameters.get("qts_subdivisions", 1),
    )

    if verbose:
        print("qlateral array complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 6
    data_assimilation_csv = data_assimilation_parameters.get(
        "data_assimilation_csv", None
    )
    data_assimilation_filter = data_assimilation_parameters.get(
        "data_assimilation_filter", None
    )
    if data_assimilation_csv or data_assimilation_filter:
        if showtiming:
            start_time = time.time()
        if verbose:
            print("creating usgs time_slice data array ...")

        usgs_df = nnu.build_data_assimilation(data_assimilation_parameters)

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

    if run_parameters.get("compute_kernel", None) == "diffusive":
        compute_func = diffusive.compute_diffusive_tst
    elif run_parameters.get("compute_method", None) == "V02-caching":
        compute_func = fast_reach.compute_network
    elif run_parameters.get("compute_method", None) == "V02-structured":
        compute_func = fast_reach.compute_network_structured
    elif run_parameters.get("compute_method", None) == "V02-structured-obj":
        compute_func = fast_reach.compute_network_structured_obj
    else:
        compute_func = fast_reach.compute_network

    # TODO: Remove below. --compute-method=V02-structured-obj did not work on command line
    # compute_func = fast_reach.compute_network_structured_obj

    results = compute_nhd_routing_v02(
        connections,
        rconn,
        wbodies,
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
        q0,
        qlats,
        usgs_df,
        run_parameters.get("assume_short_ts", False),
        run_parameters.get("return_courant", False),
        waterbodies_df_reduced,
        diffusive_parameters,
    )

    with open("mainstems_conus.txt", "w") as filehandle:
        for listitem in reaches_bytw.keys():
            filehandle.write("%s\n" % listitem)

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

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

            if run_parameters.get("return_courant", False):
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
            parity_parameters, run_parameters, results,
        )

        if verbose:
            print("parity check complete")
        if showtiming:
            print("... in %s seconds." % (time.time() - start_time))

    if verbose:
        print("process complete")
    if showtiming:
        print("%s seconds." % (time.time() - main_start_time))


if __name__ == "__main__":
    main()
