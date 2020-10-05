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
import pandas as pd
from functools import partial
from joblib import delayed, Parallel
from itertools import chain, islice
from operator import itemgetter
import logging


def _handle_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--debuglevel",
        help="Set the debuglevel",
        dest="debuglevel",
        choices=[0, -1, -2, -3],
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
        "--assume-short-ts",
        help="Use the previous timestep value for upstream flow",
        dest="assume_short_ts",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Write output files (leave blank for no writing)",
        dest="write_output",
        action="store_true",
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
        ],
        # TODO: accept multiple or a Path (argparse Action perhaps)
        # action='append',
        # nargs=1,
        dest="supernetwork",
        default="Pocono_TEST1",
    )
    parser.add_argument("--ql", help="QLat input data", dest="ql", default=None)

    return parser.parse_args()


ENV_IS_CL = False
if ENV_IS_CL:
    root = pathlib.Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = pathlib.Path("../..").resolve()
    sys.path.append(r"../python_framework_v02")

    # TODO: automate compile for the package scripts
    # sys.path.append(r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS")

## network and reach utilities
import nhd_network_utilities_v02 as nnu
import mc_reach
import nhd_network
import nhd_io


def writetoFile(file, writeString):
    file.write(writeString)
    file.write("\n")


def constant_qlats(data, nsteps, qlat):
    q = np.full((len(data.index), nsteps), qlat, dtype="float32")
    ql = pd.DataFrame(q, index=data.index, columns=range(nsteps))
    return ql


def main():

    args = _handle_args()

    nts = 144
    # debuglevel = -1 * args.debuglevel
    # verbose = args.verbose
    showtiming = args.showtiming
    supernetwork = args.supernetwork
    break_network_at_waterbodies = args.break_network_at_waterbodies
    write_output = args.write_output
    assume_short_ts = args.assume_short_ts

    test_folder = pathlib.Path(root, "test")
    geo_input_folder = test_folder.joinpath("input", "geo")

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


    LOG.info("creating supernetwork connections set")
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
    data = nhd_io.read(network_data["geo_file_path"])
    data = data[list(cols.values())]
    data = data.set_index(cols["key"])

    if "mask_file_path" in network_data:
        data_mask = nhd_io.read_mask(
            network_data["mask_file_path"],
            layer_string=network_data["mask_layer_string"],
        )
        data = data.filter(data_mask.iloc[:, network_data["mask_key"]], axis=0)

    data = data.sort_index()
    data = nhd_io.replace_downstreams(data, cols["downstream"], 0)

    if args.ql:
        qlats = nhd_io.read_qlat(args.ql)
    else:
        qlats = constant_qlats(data, nts, 10.0)

    connections = nhd_network.extract_connections(data, cols["downstream"])
    wbodies = nhd_network.extract_waterbodies(
        data, cols["waterbody"], network_data["waterbody_null_code"]
    )

    
    LOG.info("supernetwork connections set complete")
    if showtiming:
        LOG.info("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming:
        start_time = time.time()
    
    LOG.info("organizing connections into reaches ...")

    rconn = nhd_network.reverse_network(connections)
    subnets = nhd_network.reachable_network(rconn)
    subreaches = {}
    for tw, net in subnets.items():
        path_func = partial(nhd_network.split_at_junction, net)
        subreaches[tw] = nhd_network.dfs_decomposition(net, path_func)

    LOG.info("reach organization complete")
    if showtiming:
        LOG.info("... in %s seconds." % (time.time() - start_time))

    if showtiming:
        start_time = time.time()

    data["dt"] = 300.0
    data = data.rename(columns=nnu.reverse_dict(cols))
    data = data.astype("float32")

    # datasub = data[['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']]

    parallelcompute = False
    if parallelcompute:
        with Parallel(n_jobs=-1, backend="threading") as parallel:
            jobs = []
            for twi, (tw, reach) in enumerate(subreaches.items(), 1):
                r = list(chain.from_iterable(reach))
                data_sub = data.loc[
                    r, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
                ].sort_index()
                qlat_sub = qlats.loc[r].sort_index()
                jobs.append(
                    delayed(mc_reach.compute_network)(
                        nts,
                        reach,
                        subnets[tw],
                        data_sub.index.values,
                        data_sub.columns.values,
                        data_sub.values,
                        qlat_sub.values,
                    )
                )
            results = parallel(jobs)
    else:
        results = []
        for twi, (tw, reach) in enumerate(subreaches.items(), 1):
            r = list(chain.from_iterable(reach))
            data_sub = data.loc[
                r, ["dt", "bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0"]
            ].sort_index()
            qlat_sub = qlats.loc[r].sort_index()
            results.append(
                mc_reach.compute_network(
                    nts,
                    reach,
                    subnets[tw],
                    data_sub.index.values,
                    data_sub.columns.values,
                    data_sub.values,
                    qlat_sub.values,
                )
            )

    fdv_columns = pd.MultiIndex.from_product(
        [range(nts), ["q", "v", "d"]]
    ).to_flat_index()
    flowveldepth = pd.concat(
        [pd.DataFrame(d, index=i, columns=fdv_columns) for i, d in results], copy=False
    )
    flowveldepth = flowveldepth.sort_index()
    flowveldepth.to_csv(f"{args.supernetwork}.csv")
    print(flowveldepth)


    LOG.info("ordered reach computation complete")
    if showtiming:
        LOG.info("... in %s seconds." % (time.time() - start_time))


if __name__ == "__main__":
    args = _handle_args()
    debuglevel = -1 * int(args.debuglevel)
    verbose = args.verbose
    # create LOG
    # logging.basicConfig(filename='INFO.log',level=logging.DEBUG)
    LOG = logging.getLogger('log')
    # # switch to debug for all, warning gives minor printouts
    if  verbose:
        LOG.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
    elif debuglevel == 1:
        LOG.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
    elif debuglevel == 2:
        LOG.setLevel(logging.WARNING)
        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)
    else:
        LOG.setLevel(logging.CRITICAL)
        ch = logging.StreamHandler()
        ch.setLevel(logging.CRITICAL)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to LOG
    LOG.addHandler(ch)
    main()
