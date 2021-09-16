from collections import defaultdict
from itertools import chain
from functools import partial
from tlz import concat
from joblib import delayed, Parallel
from datetime import datetime, timedelta
import time
import pandas as pd
import numpy as np

import troute.nhd_network as nhd_network
from troute.routing.fast_reach.mc_reach import (
    compute_network,
    compute_network_structured,
    compute_network_structured_obj,
)
from troute.routing.fast_reach import diffusive
from troute.routing.fast_reach import diffusive_cnt

_compute_func_map = defaultdict(
    compute_network,
    {
        "diffusive_cnt": diffusive_cnt.compute_diffusive_tst,
        "diffusive": diffusive.compute_diffusive_tst,
        "V02-caching": compute_network,
        "V02-diffusive-dummy": compute_network,
        "V02-structured": compute_network_structured,
        "V02-structured-obj": compute_network_structured_obj,
    },
)


def _format_qlat_start_time(qlat_start_time):
    if not isinstance(qlat_start_time,datetime):
        try:
            return datetime.strptime(qlat_start_time, '%Y-%m-%d %H:%M:%S')
        except:  # TODO: make sure this doesn't introduce a silent error
            return datetime.now()

    else:
        return qlat_start_time


def _build_reach_type_list(reach_list, wbodies_segs):

    reach_type_list = [
                1 if (set(reaches) & wbodies_segs) else 0 for reaches in reach_list
            ]

    return list(zip(reach_list, reach_type_list))


def _prep_da_dataframes(
    usgs_df,
    lastobs_df,
    param_df_sub_idx,
    exclude_segments=None,
    ):
    """
    Produce, based on the segments in the param_df_sub_idx (which is a subset
    representing a subnetwork of the larger collection of all segments),
    a subset of the relevant usgs gage observation time series
    and the relevant last-valid gage observation from any
    prior model execution.
    
    exclude_segments (list): segments to exclude from param_df_sub when searching for gages
                             This catches and excludes offnetwork upstreams segments from being
                             realized as locations for DA substitution. Else, by-subnetwork
                             parallel executions fail.

    Cases to consider:
    USGS_DF, LAST_OBS
    Yes, Yes: Analysis and Assimilation; Last_Obs used to fill gaps in the front of the time series
    No, Yes: Forecasting mode;
    Yes, No; Cold-start case;
    No, No: Open-Loop;

    For both cases where USGS_DF is present, there is a sub-case where the length of the observed
    time series is as long as the simulation.

    """
    
    subnet_segs = param_df_sub_idx
    # segments in the subnetwork ONLY, no offnetwork upstreams included
    if exclude_segments:
        subnet_segs = param_df_sub_idx.difference(set(exclude_segments))
    
    # NOTE: Uncomment to easily test no observations...
    # usgs_df = pd.DataFrame()
    if not usgs_df.empty and not lastobs_df.empty:
        # index values for last obs are not correct, but line up correctly with usgs values. Switched
        lastobs_segs = list(lastobs_df.index.intersection(subnet_segs))
        lastobs_df_sub = lastobs_df.loc[lastobs_segs]
        usgs_segs = list(usgs_df.index.intersection(subnet_segs))
        da_positions_list_byseg = param_df_sub_idx.get_indexer(usgs_segs)
        usgs_df_sub = usgs_df.loc[usgs_segs]
    elif usgs_df.empty and not lastobs_df.empty:
        lastobs_segs = list(lastobs_df.index.intersection(subnet_segs))
        lastobs_df_sub = lastobs_df.loc[lastobs_segs]
        # Create a completely empty list of gages -- the .shape[1] attribute
        # will be == 0, and that will trigger a reference to the lastobs.
        # in the compute kernel below.
        usgs_df_sub = pd.DataFrame(index=lastobs_df_sub.index,columns=[])
        usgs_segs = lastobs_segs
        da_positions_list_byseg = param_df_sub_idx.get_indexer(lastobs_segs)
    elif not usgs_df.empty and lastobs_df.empty:
        usgs_segs = list(usgs_df.index.intersection(subnet_segs))
        da_positions_list_byseg = param_df_sub_idx.get_indexer(usgs_segs)
        usgs_df_sub = usgs_df.loc[usgs_segs]
        lastobs_df_sub = pd.DataFrame(index=usgs_df_sub.index,columns=["discharge","time","model_discharge"])
    else:
        usgs_df_sub = pd.DataFrame()
        lastobs_df_sub = pd.DataFrame()
        da_positions_list_byseg = []

    return usgs_df_sub, lastobs_df_sub, da_positions_list_byseg


def _prep_da_positions_byreach(reach_list, gage_index):
    """
    produce a list of indexes of the reach_list identifying reaches with gages
    and a corresponding list of indexes of the gage_list of the gages in
    the order they are found in the reach_list.
    """
    reach_key = []
    reach_gage = []
    for i, r in enumerate(reach_list):
        for s in r:
            if s in gage_index:
                reach_key.append(i)
                reach_gage.append(s)
    gage_reach_i = gage_index.get_indexer(reach_gage)

    return reach_key, gage_reach_i


def compute_nhd_routing_v02(
    connections,
    rconn,
    wbody_conn,
    reaches_bytw,
    compute_func_name,
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
    lastobs_df,
    da_parameter_dict,
    assume_short_ts,
    return_courant,
    waterbodies_df,
    waterbody_parameters,
    waterbody_types_df,
    waterbody_type_specified,
    diffusive_parameters=None,
):

    da_decay_coefficient = da_parameter_dict.get("da_decay_coefficient", 0)
    param_df["dt"] = dt
    param_df = param_df.astype("float32")

    start_time = time.time()
    compute_func = _compute_func_map[compute_func_name]
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
                if waterbodies_df.empty:
                    path_func = partial(nhd_network.split_at_junction, rconn_subn)
                else:
                    path_func = partial(
                        nhd_network.split_at_waterbodies_and_junctions,
                        set(waterbodies_df.index.values),
                        rconn_subn
                        )
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
                    
                    common_segs = list(param_df.index.intersection(segs))
                    wbodies_segs = set(segs).symmetric_difference(common_segs)
                    
                    #Declare empty dataframe
                    waterbody_types_df_sub = pd.DataFrame()
                
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
                        
                        #If reservoir types other than Level Pool are active
                        if not waterbody_types_df.empty:
                            waterbody_types_df_sub = waterbody_types_df.loc[
                                lake_segs,
                                [
                                    "reservoir_type",
                                ],
                            ]

                    else:
                        lake_segs = []
                        waterbodies_df_sub = pd.DataFrame()
                    
                    param_df_sub = param_df.loc[
                        common_segs,
                        ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
                    
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

                    subn_reach_list = clustered_subns["subn_reach_list"]
                    upstreams = clustered_subns["upstreams"]

                    subn_reach_list_with_type = _build_reach_type_list(subn_reach_list, wbodies_segs)

                    qlat_sub = qlats.loc[param_df_sub.index]
                    q0_sub = q0.loc[param_df_sub.index]
                    
                    #Determine model_start_time from qlat_start_time
                    qlat_start_time = list(qlat_sub)[0]

                    qlat_time_step_seconds = qts_subdivisions * dt

                    qlat_start_time_datetime_object =  _format_qlat_start_time(qlat_start_time)

                    model_start_time_datetime_object = qlat_start_time_datetime_object \
                    - timedelta(seconds=qlat_time_step_seconds)

                    model_start_time = model_start_time_datetime_object.strftime('%Y-%m-%d_%H:%M:%S')

                    param_df_sub = param_df_sub.reindex(
                        param_df_sub.index.tolist() + lake_segs
                    ).sort_index()

                    usgs_df_sub, lastobs_df_sub, da_positions_list_byseg = _prep_da_dataframes(usgs_df, lastobs_df, param_df_sub.index, offnetwork_upstreams)
                    da_positions_list_byreach, da_positions_list_bygage = _prep_da_positions_byreach(subn_reach_list, lastobs_df_sub.index)

                    qlat_sub = qlat_sub.reindex(param_df_sub.index)
                    q0_sub = q0_sub.reindex(param_df_sub.index)

                    # results_subn[order].append(
                    #     compute_func(
                    jobs.append(
                        delayed(compute_func)(
                            nts,
                            dt,
                            qts_subdivisions,
                            subn_reach_list_with_type,
                            upstreams,
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            q0_sub.values.astype("float32"),
                            qlat_sub.values.astype("float32"),
                            lake_segs, 
                            waterbodies_df_sub.values,
                            waterbody_parameters,
                            waterbody_types_df_sub.values.astype("int32"),
                            waterbody_type_specified,
                            model_start_time,
                            usgs_df_sub.values.astype("float32"),
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
                            np.array(da_positions_list_byseg, dtype="int32"),
                            np.array(da_positions_list_byreach, dtype="int32"),
                            np.array(da_positions_list_bygage, dtype="int32"),
                            lastobs_df_sub.get(
                                "lastobs_discharge",
                                pd.Series(index=lastobs_df_sub.index, name="Null"),
                            ).values.astype("float32"),
                            lastobs_df_sub.get(
                                "time_since_lastobs",
                                pd.Series(index=lastobs_df_sub.index, name="Null"),
                            ).values.astype("float32"),
                            da_decay_coefficient,
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
                if waterbodies_df.empty:
                    path_func = partial(nhd_network.split_at_junction, rconn_subn)
                else:
                    path_func = partial(
                        nhd_network.split_at_waterbodies_and_junctions,
                        set(waterbodies_df.index.values),
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
                    
                    common_segs = list(param_df.index.intersection(segs))
                    wbodies_segs = set(segs).symmetric_difference(common_segs)
                    
                    #Declare empty dataframe
                    waterbody_types_df_sub = pd.DataFrame()
                
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
                        
                        #If reservoir types other than Level Pool are active
                        if not waterbody_types_df.empty:
                            waterbody_types_df_sub = waterbody_types_df.loc[
                                lake_segs,
                                [
                                    "reservoir_type",
                                ],
                            ]

                    else:
                        lake_segs = []
                        waterbodies_df_sub = pd.DataFrame()
                    
                    param_df_sub = param_df.loc[
                        common_segs,
                        ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
                    
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

                    subn_reach_list_with_type = _build_reach_type_list(subn_reach_list, wbodies_segs)

                    qlat_sub = qlats.loc[param_df_sub.index]
                    q0_sub = q0.loc[param_df_sub.index]
                    
                    #Determine model_start_time from qlat_start_time
                    qlat_start_time = list(qlat_sub)[0]

                    qlat_time_step_seconds = qts_subdivisions * dt

                    qlat_start_time_datetime_object =  _format_qlat_start_time(qlat_start_time)

                    model_start_time_datetime_object = qlat_start_time_datetime_object \
                    - timedelta(seconds=qlat_time_step_seconds)

                    model_start_time = model_start_time_datetime_object.strftime('%Y-%m-%d_%H:%M:%S')

                    param_df_sub = param_df_sub.reindex(
                        param_df_sub.index.tolist() + lake_segs
                    ).sort_index()

                    usgs_df_sub, lastobs_df_sub, da_positions_list_byseg = _prep_da_dataframes(usgs_df, lastobs_df, param_df_sub.index, offnetwork_upstreams)
                    da_positions_list_byreach, da_positions_list_bygage = _prep_da_positions_byreach(subn_reach_list, lastobs_df_sub.index)

                    qlat_sub = qlat_sub.reindex(param_df_sub.index)
                    q0_sub = q0_sub.reindex(param_df_sub.index)

                    jobs.append(
                        delayed(compute_func)(
                            nts,
                            dt,
                            qts_subdivisions,
                            subn_reach_list_with_type,
                            subnetworks[subn_tw],
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            q0_sub.values.astype("float32"),
                            qlat_sub.values.astype("float32"),
                            lake_segs,
                            waterbodies_df_sub.values,
                            waterbody_parameters,
                            waterbody_types_df_sub.values.astype("int32"),
                            waterbody_type_specified,
                            model_start_time,
                            usgs_df_sub.values.astype("float32"),
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
                            np.array(da_positions_list_byseg, dtype="int32"),
                            np.array(da_positions_list_byreach, dtype="int32"),
                            np.array(da_positions_list_bygage, dtype="int32"),
                            lastobs_df_sub.get(
                                "lastobs_discharge",
                                pd.Series(index=lastobs_df_sub.index, name="Null"),
                            ).values.astype("float32"),
                            lastobs_df_sub.get(
                                "time_since_lastobs",
                                pd.Series(index=lastobs_df_sub.index, name="Null"),
                            ).values.astype("float32"),
                            da_decay_coefficient,
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

    elif parallel_compute_method == "by-subnetwork-diffusive":
        reaches_ordered_bysubntw, subnetworks, subnetworks_only_ordered_jit = nhd_network.build_subnetworks_btw_reservoirs(
            connections, rconn, wbody_conn, independent_networks, sources=None
        )

        if 1 == 1:
            print("JIT Preprocessing time %s seconds." % (time.time() - start_time))
            print("starting Parallel JIT calculation")

        start_para_time = time.time()
        with Parallel(n_jobs=cpu_pool, backend="threading") as parallel:
            results_subn = defaultdict(list)
            flowveldepth_interorder = {}

            for order in range(max(reaches_ordered_bysubntw.keys()), -1, -1):
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

                    common_segs = list(param_df.index.intersection(segs))
                    wbodies_segs = set(segs).symmetric_difference(common_segs)
                    
                    # Declare empty dataframe
                    waterbody_types_df_sub = pd.DataFrame()

                    # Set compute_func_switch to compute_func.
                    # compute_func for this function should be set to "diffusive" 
                    compute_func_switch = compute_func

                    # Can comment out above statement and uncomment below
                    # if need to run compute_network_structured in mc_reach
                    # for every subnetwork for debugging purposes.
                    #compute_func_switch = compute_network_structured

                    if not waterbodies_df.empty:
                        lake_segs = list(waterbodies_df.index.intersection(segs))

                        if subn_tw in waterbodies_df.index:
                            # Since this subn_tw is a resevoir, set compute_func_switch
                            # to compute_network_structured.
                            compute_func_switch = compute_network_structured

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
                        
                        #If reservoir types other than Level Pool are active
                        if not waterbody_types_df.empty:
                            waterbody_types_df_sub = waterbody_types_df.loc[
                                lake_segs,
                                [
                                    "reservoir_type",
                                ],
                            ]

                    else:
                        lake_segs = []
                        waterbodies_df_sub = pd.DataFrame()
                    
                    param_df_sub = param_df.loc[
                        common_segs,
                        ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                    ].sort_index()
                    
                    param_df_sub_super = param_df_sub.reindex(
                        param_df_sub.index.tolist() + lake_segs
                    ).sort_index()

                    if order < max(reaches_ordered_bysubntw.keys()):
                        for us_subn_tw in offnetwork_upstreams:
                            subn_tw_sortposition = param_df_sub_super.index.get_loc(
                                us_subn_tw
                            )
                            flowveldepth_interorder[us_subn_tw][
                                "position_index"
                            ] = subn_tw_sortposition

                    subn_reach_list_with_type = _build_reach_type_list(subn_reach_list, wbodies_segs)

                    qlat_sub = qlats.loc[param_df_sub.index]
                    q0_sub = q0.loc[param_df_sub.index]
                    
                    #Determine model_start_time from qlat_start_time
                    qlat_start_time = list(qlat_sub)[0]

                    qlat_time_step_seconds = qts_subdivisions * dt

                    qlat_start_time_datetime_object =  _format_qlat_start_time(qlat_start_time)

                    model_start_time_datetime_object = qlat_start_time_datetime_object \
                    - timedelta(seconds=qlat_time_step_seconds)

                    model_start_time = model_start_time_datetime_object.strftime('%Y-%m-%d_%H:%M:%S')

                    param_df_sub = param_df_sub.reindex(
                        param_df_sub.index.tolist() + lake_segs
                    ).sort_index()

                    usgs_df_sub, lastobs_df_sub, da_positions_list_byseg = _prep_da_dataframes(usgs_df, lastobs_df, param_df_sub.index, offnetwork_upstreams)
                    da_positions_list_byreach, da_positions_list_bygage = _prep_da_positions_byreach(subn_reach_list, lastobs_df_sub.index)

                    qlat_sub = qlat_sub.reindex(param_df_sub.index)
                    q0_sub = q0_sub.reindex(param_df_sub.index)

                    jobs.append(
                        delayed(compute_func_switch)(
                            nts,
                            dt,
                            qts_subdivisions,
                            subn_reach_list_with_type,
                            subnetworks[subn_tw],
                            param_df_sub.index.values,
                            param_df_sub.columns.values,
                            param_df_sub.values,
                            q0_sub.values.astype("float32"),
                            qlat_sub.values.astype("float32"),
                            lake_segs,
                            waterbodies_df_sub.values,
                            waterbody_parameters,
                            waterbody_types_df_sub.values.astype("int32"),
                            waterbody_type_specified,
                            model_start_time,
                            usgs_df_sub.values.astype("float32"),
                            np.array(da_positions_list_byseg, dtype="int32"),
                            np.array(da_positions_list_byreach, dtype="int32"),
                            np.array(da_positions_list_bygage, dtype="int32"),
                            lastobs_df_sub.get("lastobs_discharge", pd.Series(index=lastobs_df_sub.index, name="Null")).values.astype("float32"),
                            lastobs_df_sub.get("time_since_lastobs", pd.Series(index=lastobs_df_sub.index, name="Null")).values.astype("float32"),
                            da_decay_coefficient,
                            # flowveldepth_interorder,  # obtain keys and values from this dataset
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

                #Declare empty dataframe
                waterbody_types_df_sub = pd.DataFrame()

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

                    #If reservoir types other than Level Pool are active
                    if not waterbody_types_df.empty:
                        waterbody_types_df_sub = waterbody_types_df.loc[
                            lake_segs,
                            [
                                "reservoir_type",
                            ],
                        ]

                else:
                    lake_segs = []
                    waterbodies_df_sub = pd.DataFrame()

                param_df_sub = param_df.loc[
                    common_segs,
                    ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
                ].sort_index()

                reaches_list_with_type = _build_reach_type_list(reach_list, wbodies_segs)

                # qlat_sub = qlats.loc[common_segs].sort_index()
                # q0_sub = q0.loc[common_segs].sort_index()
                qlat_sub = qlats.loc[param_df_sub.index]
                q0_sub = q0.loc[param_df_sub.index]

                #Determine model_start_time from qlat_start_time
                qlat_start_time = list(qlat_sub)[0]

                qlat_time_step_seconds = qts_subdivisions * dt

                qlat_start_time_datetime_object =  _format_qlat_start_time(qlat_start_time)

                model_start_time_datetime_object = qlat_start_time_datetime_object \
                - timedelta(seconds=qlat_time_step_seconds)

                model_start_time = model_start_time_datetime_object.strftime('%Y-%m-%d_%H:%M:%S')

                param_df_sub = param_df_sub.reindex(
                    param_df_sub.index.tolist() + lake_segs
                ).sort_index()

                usgs_df_sub, lastobs_df_sub, da_positions_list_byseg = _prep_da_dataframes(usgs_df, lastobs_df, param_df_sub.index)
                da_positions_list_byreach, da_positions_list_bygage = _prep_da_positions_byreach(reach_list, lastobs_df_sub.index)

                qlat_sub = qlat_sub.reindex(param_df_sub.index)
                q0_sub = q0_sub.reindex(param_df_sub.index)

                jobs.append(
                    delayed(compute_func)(
                        nts,
                        dt,
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
                        waterbody_types_df_sub.values.astype("int32"),
                        waterbody_type_specified,
                        model_start_time,
                        usgs_df_sub.values.astype("float32"),
                        np.array(da_positions_list_byseg, dtype="int32"),
                        np.array(da_positions_list_byreach, dtype="int32"),
                        np.array(da_positions_list_bygage, dtype="int32"),
                        lastobs_df_sub.get("lastobs_discharge", pd.Series(index=lastobs_df_sub.index, name="Null")).values.astype("float32"),
                        lastobs_df_sub.get("time_since_lastobs", pd.Series(index=lastobs_df_sub.index, name="Null")).values.astype("float32"),
                        da_decay_coefficient,
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

            #Declare empty dataframe
            waterbody_types_df_sub = pd.DataFrame()

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

                #If reservoir types other than Level Pool are active
                if not waterbody_types_df.empty:
                    waterbody_types_df_sub = waterbody_types_df.loc[
                        lake_segs,
                        [
                            "reservoir_type",
                        ],
                    ]

            else:
                lake_segs = []
                waterbodies_df_sub = pd.DataFrame()

            param_df_sub = param_df.loc[
                common_segs,
                ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"],
            ].sort_index()

            reaches_list_with_type = _build_reach_type_list(reach_list, wbodies_segs)

            # qlat_sub = qlats.loc[common_segs].sort_index()
            # q0_sub = q0.loc[common_segs].sort_index()
            qlat_sub = qlats.loc[param_df_sub.index]
            q0_sub = q0.loc[param_df_sub.index]

            #Determine model_start_time from qlat_start_time
            qlat_start_time = list(qlat_sub)[0]

            qlat_time_step_seconds = qts_subdivisions * dt

            qlat_start_time_datetime_object =  _format_qlat_start_time(qlat_start_time)

            model_start_time_datetime_object = qlat_start_time_datetime_object \
            - timedelta(seconds=qlat_time_step_seconds)

            model_start_time = model_start_time_datetime_object.strftime('%Y-%m-%d_%H:%M:%S')

            param_df_sub = param_df_sub.reindex(
                param_df_sub.index.tolist() + lake_segs
            ).sort_index()

            usgs_df_sub, lastobs_df_sub, da_positions_list_byseg = _prep_da_dataframes(usgs_df, lastobs_df, param_df_sub.index)
            da_positions_list_byreach, da_positions_list_bygage = _prep_da_positions_byreach(reach_list, lastobs_df_sub.index)

            qlat_sub = qlat_sub.reindex(param_df_sub.index)
            q0_sub = q0_sub.reindex(param_df_sub.index)

            results.append(
                compute_func(
                    nts,
                    dt,
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
                    waterbody_types_df_sub.values.astype("int32"),
                    waterbody_type_specified,
                    model_start_time,
                    usgs_df_sub.values.astype("float32"),
                    np.array(da_positions_list_byseg, dtype="int32"),
                    np.array(da_positions_list_byreach, dtype="int32"),
                    np.array(da_positions_list_bygage, dtype="int32"),
                    lastobs_df_sub.get("lastobs_discharge", pd.Series(index=lastobs_df_sub.index, name="Null")).values.astype("float32"),
                    lastobs_df_sub.get("time_since_lastobs", pd.Series(index=lastobs_df_sub.index, name="Null")).values.astype("float32"),
                    da_decay_coefficient,
                    {},
                    assume_short_ts,
                    return_courant,
                    diffusive_parameters,
                )
            )

    return results
