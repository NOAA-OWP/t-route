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
from itertools import chain
from functools import partial
from joblib import delayed, Parallel

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL:
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.path.append(r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS')

## Muskingum Cunge
COMPILE = False
if COMPILE:
    import subprocess
    print("Recompiling Fortran module")
    fortran_compile_call = ['f2py3', '-c', 'MCsingleSegStime_f2py_NOLOOP.f90', '-m', 'mc_sseg_stime_NOLOOP']
    subprocess.run(fortran_compile_call, cwd=r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
#from mc_sseg_stime_NOLOOP import muskingcungenwm
#import mc_sseg_stime_NOLOOP as mc
#import single_seg as mc
#from mc_reach import muskingcunge as muskingcungenwm
import mc_reach

#muskingcungenwm = mc.muskingcunge_module.muskingcungenwm

data_values = None
WRITE_OUTPUT = False

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru
import nhd_io
import nhd_network


def writetoFile(file, writeString):
    file.write(writeString)
    file.write('\n')

first = itemgetter(0)
second = itemgetter(1)

# def compute_reach(
#         reach,
#         qup, quc,
#         flowdepthvel,
#         assume_short_ts=False,
# ):
#     for i, current_segment in reach:
#         dt = 60.0
#         bw = data_values[current_segment, 0]
#         tw = data_values[current_segment, 1]
#         twcc = data_values[current_segment, 2]
#         dx = data_values[current_segment, 3]
#         n_manning = data_values[current_segment, 4]
#         n_manning_cc = data_values[current_segment, 5]
#         cs = data_values[current_segment, 6]
#         s0 = data_values[current_segment, 7]
#
#         flowdepthvel[i, 7] = qlat = 10.0
#         #qdp = flowdepthvel[current_segment, 0]
#         #depthp = flowdepthvel[current_segment, 1]
#         #velp = flowdepthvel[current_segment, 2]
#         qdp, depthp, velp = flowdepthvel[i, 0:3]
#
#         flowdepthvel[i, :4] = flowdepthvel[i, 4:]
#
#         if assume_short_ts:
#             quc = qup
#
#         # qdc, velc, depthc
#         #args = tuple(map(c_float, (dt, qup, quc, qdp, qlat, dx, bw, tw, twcc, n_manning)))
#         qdc, depthc, velc = muskingcungenwm(
#             dt,
#             qup,
#             quc,
#             qdp,
#             qlat,
#             dx,
#             bw,
#             tw,
#             twcc,
#             n_manning,
#             n_manning_cc,
#             cs,
#             s0,
#             velp,
#             depthp
#         )
#
#         flowdepthvel[i, 4:7] = qdc, depthc, velc
#         quc = qdc
#         qup = qdp
#
#
# # ### Psuedocode
# #
#
# def compute_network(nsteps, reaches, connections, assume_short_ts=False):
#     findex = np.array(sorted(chain.from_iterable(reaches), key=first))
#     flowdepthvel = np.zeros((len(findex), 8), dtype='float32')
#     qup_quc = np.zeros(2, dtype='float32')
#
#     idx_view = findex[:,0]
#
#     global_idxs = []
#     local_idxs = []
#     upstream_segs = []
#     for reach in reaches:
#         x, y = list(zip(*reach))
#         global_idxs.append(x)
#         local_idxs.append(np.searchsorted(idx_view, x).tolist())
#         us_segs = list(map(first, connections.get(y[0], {}).get('children', ())))
#         upstream_segs.append(np.searchsorted(idx_view, us_segs).tolist())
#
#     for ts in range(nsteps):
#         #breakpoint()
#         for local_idx, global_idx, us in zip(local_idxs, global_idxs, upstream_segs):
#             qup_quc[:] = 0
#             for x in us:
#                 qup_quc[0] += flowdepthvel[x, 0]
#                 qup_quc[1] += flowdepthvel[x, 4]
#             mc_reach.compute_reach(zip(local_idx, global_idx), qup_quc[0], qup_quc[1], flowdepthvel, data_values, assume_short_ts=assume_short_ts)
#     return idx_view, flowdepthvel


def translate_reach_to_index(reaches, index):
    """
    Translate the reach value to an index into data
    Args:
        reaches (list): list of reaches
        data (DataFrame): data

    Returns:

    """
    rv = []
    for reach in reaches:
        rv.append(np.stack([np.searchsorted(index, reach), reach], axis=1))
    return rv

def translate_network_to_index(connections, index):
    """
    Annotate a connections dictionary with indexes for each node into index
    Args:
        connections:
        index:

    Returns:

    """
    rv = {}
    for n, children in connections.items():
        r = {}
        if children:
            r['children'] = np.stack([np.searchsorted(index, children), children], axis=1)
        else:
            r['children'] = None
        try:
            r["index"] = index.get_loc(n)
        except KeyError:
            r["index"] = None
        rv[n] = r
    return rv

def main():


    global data_values

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo')

    # TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    #supernetwork = 'Brazos_LowerColorado_ge5'
    supernetwork = 'Pocono_TEST2'
    """##NHD CONUS order 5 and greater"""
    #supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    # supernetwork = 'Mainstems_CONUS'
    # supernetwork = 'CONUS_FULL_RES_v20'
    # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
    # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches

    if verbose: print('creating supernetwork connections set')
    if showtiming: start_time = time.time()
    # STEP 1
    network_data = nnu.set_supernetwork_data(supernetwork=supernetwork,
                                             geo_input_folder=geo_input_folder)

    cols = [v for c, v in network_data.items() if c.endswith("_col")]
    data = nhd_io.read(network_data['geo_file_path'])
    data = data[cols]
    data = data.set_index(network_data['key_col'])

    if 'mask_file_path' in network_data:
        data_mask = nhd_io.read_mask(network_data['mask_file_path'],
                                     layer_string=network_data['mask_layer_string'])
        data = data.filter(data_mask.iloc[:, network_data["mask_key"]], axis=0)

    data = data.sort_index()
    connections = nhd_network.extract_network(data, network_data["downstream_col"])
    rconn = nhd_network.reverse_network(connections)
    #rconn_annoated = translate_network_to_index(rconn, data.index)
    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into reaches ...')
    #reaches = nhd_network.dfs_decomposition(connections)
    #reaches_i = translate_reach_to_index(reaches, data.index)

    subreachable = nhd_network.reachable(rconn)
    subnets_ = nhd_network.reachable_network(rconn)
    subreaches = {}
    for tw, net in subnets_.items():
        path_func = partial(nhd_network.split_at_junction, net)
        reach = nhd_network.dfs_decomposition(nhd_network.reverse_network(net), path_func)
        subreaches[tw] = translate_reach_to_index(reach, data.index)
    subnets = {k: translate_network_to_index(v, data.index) for k, v in subnets_.items()}

    #networks = nru.compose_networks(
    #    supernetwork_values
    #    , verbose=False
    #    # , verbose = verbose
    #    , debuglevel=debuglevel
    #    , showtiming=showtiming
    #)
    if verbose: print('reach organization complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    if showtiming: start_time = time.time()

    flowdepthvel = np.zeros((len(data), 8), dtype='float32')
    #flowdepthvel = {node: {'flow': {'prev': 0, 'curr': 0}
    #    , 'depth': {'prev': 0, 'curr': 0}
    #    , 'vel': {'prev': 0, 'curr': 0}
    #    , 'qlat': {'prev': 0, 'curr': 0}} for node in nhd_network.nodes(connections)}

    parallelcompute = False

    # Data column ordering is very important as we directly lookup values.
    # The column order *must* be:
    # 0: bw, 1: tw, 2: twcc, 3: dx, 4: n_manning 5: n_manning_cc, 6: cs, 7: s0
    datasub = data[['BtmWdth', 'TopWdth', 'TopWdthCC', 'Length', 'n', 'nCC', 'ChSlp', 'So']]
    data_values = datasub.values.astype('float32')
    ts = 144
    compute_start = time.time()
    if parallelcompute:
        if verbose: print('executing computation on ordered reaches ...')
        with Parallel(n_jobs=5, pre_dispatch='all', backend='threading', verbose=5) as parallal:
            jobs = []
            for twi, (tw, reach) in enumerate(subreaches.items(), 1):
                jobs.append(delayed(mc_reach.compute_network)(ts, reach, subnets[tw], data_values))
            for findex, fdv in parallal(jobs):
                flowdepthvel[findex] = fdv


    else:
        for twi, (tw, reach) in enumerate(subreaches.items(), 1):
            findex, fdv = mc_reach.compute_network(ts, reach, subnets[tw], data_values)
            flowdepthvel[findex] = fdv

            if showtiming:
                print(f"... in {time.time()-compute_start} seconds ({twi}/{len(subreaches)})")

    print("Computation time: ", time.time() - compute_start)
    with np.printoptions(precision=5, suppress=True, linewidth=180, edgeitems=5):
        print(flowdepthvel)

        print(flowdepthvel.shape)
        #print(sorted(connections.keys()))

#if __name__ == '__main__':
#    main()
main()
