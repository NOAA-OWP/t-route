#!/usr/bin/env python
# coding: utf-8


# -*- coding: utf-8 -*-
"""NHD Network traversal

A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
## Parallel execution
import multiprocessing
import os
import sys
import time
import numpy as np
from operator import itemgetter
from itertools import chain, islice
from functools import partial
from joblib import delayed, Parallel
import random
import pandas as pd

ENV_IS_CL = False
if ENV_IS_CL:
    root = "/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/"
elif not ENV_IS_CL:
    root = os.path.dirname(os.path.dirname(os.path.abspath("")))
    sys.path.append(r"../python_framework")
    sys.path.append(r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS")

## Muskingum Cunge
COMPILE = False
if COMPILE:
    import subprocess

    print("Recompiling Fortran module")
    fortran_compile_call = [
        "f2py3",
        "-c",
        "MCsingleSegStime_f2py_NOLOOP.f90",
        "-m",
        "mc_sseg_stime_NOLOOP",
    ]
    subprocess.run(
        fortran_compile_call,
        cwd=r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS",
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
# from mc_sseg_stime_NOLOOP import muskingcungenwm
# import mc_sseg_stime_NOLOOP as mc
# import single_seg as mc
# from mc_reach import muskingcunge as muskingcungenwm
import mc_reach

# muskingcungenwm = mc.muskingcunge_module.muskingcungenwm

data_values = None
WRITE_OUTPUT = False

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru
import nhd_io
import nhd_network


def writetoFile(file, writeString):
    file.write(writeString)
    file.write("\n")


first = itemgetter(0)
second = itemgetter(1)


def constant_qlats(data, nsteps, qlat):
    q = np.full((len(data.index), nsteps), qlat, dtype='float32')
    ql = pd.DataFrame(q, index=data.index, columns=range(nsteps))
    return ql

def replace_downstreams(data, downstream_col, terminal_code):
    ds0_mask = data[downstream_col] == terminal_code
    new_data = data.copy()
    new_data.loc[ds0_mask, downstream_col] = ds0_mask.index[ds0_mask]

    # Also set negative any nodes in downstream col not in data.index
    new_data.loc[~data[downstream_col].isin(data.index), downstream_col] *= -1
    return new_data

def main():

    global data_values

    verbose = True
    debuglevel = 0
    showtiming = True
    nts = 1440

    test_folder = os.path.join(root, r"test")
    geo_input_folder = os.path.join(test_folder, r"input", r"geo")

    # TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    #supernetwork = "Pocono_TEST2"
    """##NHD CONUS order 5 and greater"""
    # supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    supernetwork = 'Mainstems_CONUS'
    # supernetwork = 'CONUS_FULL_RES_v20'
    # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
    # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches

    if verbose:
        print("creating supernetwork connections set")
    if showtiming:
        start_time = time.time()
    # STEP 1
    network_data = nnu.set_supernetwork_data(
        supernetwork=supernetwork, geo_input_folder=geo_input_folder
    )

    cols = [v for c, v in network_data.items() if c.endswith("_col")]
    data = nhd_io.read(network_data["geo_file_path"])
    data = data[cols]
    data = data.set_index(network_data["key_col"])

    if "mask_file_path" in network_data:
        data_mask = nhd_io.read_mask(
            network_data["mask_file_path"],
            layer_string=network_data["mask_layer_string"],
        )
        data = data.filter(data_mask.iloc[:, network_data["mask_key"]], axis=0)

    data = data.sort_index()
    data = replace_downstreams(data, network_data['downstream_col'], 0)

    if supernetwork == "Pocono_TEST2":
        qlats = pd.read_csv('../../test/input/geo/PoconoSampleData2/Pocono_ql_testsamp1_nwm_mc.txt', index_col='ntt')
        qlats = qlats.drop(columns=['nt'])
        qlats.columns = qlats.columns.astype(int)
        qlats = qlats.sort_index(axis='columns').sort_index(axis='index')
        qlats = qlats.drop(columns=qlats.columns.difference(data.index)).T
        qlats = qlats.astype('float32')
    else:
        qlats = constant_qlats(data, nts, 10.0)


    connections = nhd_network.extract_connections(data, network_data["downstream_col"])
    rconn = nhd_network.reverse_network(connections)
    # rconn_annoated = translate_network_to_index(rconn, data.index)
    if verbose:
        print("supernetwork connections set complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming:
        start_time = time.time()
    if verbose:
        print("organizing connections into reaches ...")

    subnets = nhd_network.reachable_network(rconn)
    subreaches = {}
    for tw, net in subnets.items():
        path_func = partial(nhd_network.split_at_junction, net)
        reach = nhd_network.dfs_decomposition(
            nhd_network.reverse_network(net), path_func
        )
        subreaches[tw] = reach

    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    if showtiming:
        start_time = time.time()

    # flowdepthvel = {node: {'flow': {'prev': 0, 'curr': 0}
    #    , 'depth': {'prev': 0, 'curr': 0}
    #    , 'vel': {'prev': 0, 'curr': 0}
    #    , 'qlat': {'prev': 0, 'curr': 0}} for node in nhd_network.nodes(connections)}

    parallelcompute = False

    # Data column ordering is very important as we directly lookup values.
    # The column order *must* be:
    # 0: bw, 1: tw, 2: twcc, 3: dx, 4: n_manning 5: n_manning_cc, 6: cs, 7: s0, 8: qlat
    data['dt'] = 300.0

    data = data.rename(columns={'Length': 'dx', 'TopWdth': 'tw', 'TopWdthCC': 'twcc',
        'BtmWdth': 'bw', 'nCC': 'ncc', 'So': 's0', 'ChSlp': 'cs'})
    data = data.astype('float32')    
    #datasub = data[['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']]
    
    #qlats = qlats.loc[:, :nts]
    compute_start = time.time()
    if parallelcompute:
        if verbose:
            print("executing computation on ordered reaches ...")
        with Parallel(
            n_jobs=-1, pre_dispatch="all", backend="threading", verbose=5
        ) as parallel:
            jobs = []
            for twi, (tw, reach) in enumerate(subreaches.items(), 1):
                r = list(chain.from_iterable(reach))
                #assert r[-1] == tw  # always be True
                #assert len(data.index.intersection(r)) == len(r)
                data_sub = data.loc[r, ['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']].sort_index()
                qlat_sub = qlats.loc[r].sort_index()
                jobs.append(
                    delayed(mc_reach.compute_network)(
                        nts, reach, subnets[tw], data_sub.index.values, data_sub.columns.values, data_sub.values, qlat_sub.values
                    )
                )
            random.shuffle(jobs)
            rets = parallel(jobs)
            #for findex, fdv in rets:
            #    flowdepthvel[findex] = fdv

    else:
        rets = []
        for twi, (tw, reach) in enumerate(subreaches.items(), 1):
            r = list(chain.from_iterable(reach))
            #assert r[-1] == tw  # always be True
            #assert len(data.index.intersection(r)) == len(r)
            data_sub = data.loc[r, ['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']].sort_index()
            #TODO: Do the data_sub = data.loc as a preprocessing step, then dump the pointer to data
            qlat_sub = qlats.loc[r].sort_index()
            rets.append(
                mc_reach.compute_network(
                    nts, reach, subnets[tw], data_sub.index.values, data_sub.columns.values, data_sub.values, qlat_sub.values))
            # TODO: rets could be dumped to files
            #findex, fdv = mc_reach.compute_network(ts, reach, subnets[tw], data_idx, data_values, qlat_values)
            #flowdepthvel[findex] = fdv
            if verbose:
                print(
                    f"tailwater: {tw} completed",
                    end = "", 
                )
            # NOTE: Mississippi River tailwater is {22811611,}:

            if showtiming:
                print(
                    f"... in {time.time()-compute_start} seconds ({twi}/{len(subreaches)})"
                )

    print("Computation time: ", time.time() - compute_start)
    fdv_columns = pd.MultiIndex.from_product([range(nts), ['q', 'v', 'd']], names=['timestep', 'qvd'])
    flowveldepth = pd.concat([pd.DataFrame(d, index=i, columns=fdv_columns) for i, d in rets])
    print(flowveldepth)
    #with np.printoptions(precision=6, suppress=True, linewidth=180, edgeitems=5):
    #    print(flowdepthvel)


# if __name__ == '__main__':
#    main()
main()
