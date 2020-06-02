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
import pathlib
import sys
import time
import numpy as np
from operator import itemgetter
from cytoolz import pluck

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL:
    root = os.path.dirname(os.path.dirname(os.path.abspath('')))
    sys.path.append(r'../python_framework')
    sys.path.append(r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS')
    sys.setrecursionlimit(4000)

## Muskingum Cunge
COMPILE = True
if COMPILE:
    import subprocess
    print("Recompiling Fortran module")
    fortran_compile_call = ['f2py3', '-c', 'MCsingleSegStime_f2py_NOLOOP.f90', '-m', 'mc_sseg_stime_NOLOOP']
    subprocess.run(fortran_compile_call, cwd=r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
from mc_sseg_stime_NOLOOP import muskingcunge_module as mc
#import single_seg as mc
#import mc_reach

connections = None
networks = None
flowdepthvel = None
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

def compute_reach(
        reach,
        qup, quc,
        data,
        flowdepthvel,
        assume_short_ts=False,
):

    for current_segment in reach:
        dt = 60.0
        bw = data[current_segment, 0]
        tw = data[current_segment, 1]
        twcc = data[current_segment, 2]
        dx = data[current_segment, 3]
        n_manning = data[current_segment, 4]
        n_manning_cc = data[current_segment, 5]
        cs = data[current_segment, 6]
        s0 = data[current_segment, 7]

        flowdepthvel[current_segment, 7] = qlat = 10.0
        #qdp = flowdepthvel[current_segment, 0]
        #depthp = flowdepthvel[current_segment, 1]
        #velp = flowdepthvel[current_segment, 2]
        qdp, depthp, velp = flowdepthvel[current_segment, 0:3]

        flowdepthvel[current_segment, :4] = flowdepthvel[current_segment, 4:]

        if assume_short_ts:
            quc = qup

        # qdc, velc, depthc
        #args = tuple(map(c_float, (dt, qup, quc, qdp, qlat, dx, bw, tw, twcc, n_manning)))
        qdc, velc, depthc = mc.muskingcungenwm(
            dt,
            qup,
            quc,
            qdp,
            qlat,
            dx,
            bw,
            tw,
            twcc,
            n_manning,
            n_manning_cc,
            cs,
            s0,
            velp,
            depthp
        )

        flowdepthvel[current_segment, 4:7] = qdc, depthc, velc
        quc = qdc
        qup = qdp


# ### Psuedocode
# 



# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
        reach=None,
        connections=None,
        supernetwork_data=None,
        cols=None,
        data=None,
        ts=0,
        flowdepthvel=None,
        verbose=False,
        debuglevel=0,
        write_output=False,
        assume_short_ts=False,
):


    qup = 0.0
    quc = 0.0
    # import pdb; pdb.set_trace()
    for us in map(first, connections.get(reach[0][1], {}).get('children', ())):
        qup += flowdepthvel[us, 0]
        quc += flowdepthvel[us, 4]

    compute_reach(pluck(0, reach), qup, quc, data, flowdepthvel, assume_short_ts=assume_short_ts)


def singlesegment(
        dt # dt
        , qup = None # qup
        , quc = None # quc
        , qdp = None # qdp
        , qlat = None # ql
        , dx = None # dx
        , bw = None # bw
        , tw = None # tw
        , twcc = None # twcc
        , n_manning = None #
        , n_manning_cc = None # ncc
        , cs = None # cs
        , s0 = None # s0
        , velp = None # velocity at previous time step
        , depthp = None # depth at previous time step
    ):

    # call Fortran routine
    return mc.muskingcungenwm(
        dt, qup, quc, qdp, qlat, dx, bw, tw, twcc
        ,n_manning, n_manning_cc, cs, s0, velp, depthp
    )
    #return qdc, vel, depth

def flow_dict_to_arr(connections, fdv):
    flowdepthvel_arr = np.zeros((len(connections), 8))
    for i, c in enumerate(sorted(connections.keys())):
        d = fdv[c]
        flowdepthvel_arr[i] = (d['flow']['prev'], d['depth']['prev'], d['vel']['prev'], d['qlat']['prev'],
                               d['flow']['curr'], d['depth']['curr'], d['vel']['curr'], d['qlat']['curr'])
    return flowdepthvel_arr


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
        rv.append(tuple(zip(np.searchsorted(index, reach), reach)))
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
        rv[n] = {"index": index.get_loc(n), "children": tuple(zip(np.searchsorted(index, children), children))}
    return rv

def main():


    global connections
    global networks
    global flowdepthvel

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo')

    # TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    #supernetwork = 'Brazos_LowerColorado_ge5'
    #supernetwork = 'Pocono_TEST1'
    """##NHD CONUS order 5 and greater"""
    #supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    supernetwork = 'Mainstems_CONUS'
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
    rconn_annoated = translate_network_to_index(rconn, data.index)
    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into reaches ...')
    reaches = nhd_network.dfs_decomposition(connections)
    reaches_i = translate_reach_to_index(reaches, data.index)
    del connections
    del rconn
    del reaches

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
    if not parallelcompute:
        if verbose: print('executing computation on ordered reaches ...')

        data_values= data.values.astype('float32')
        compute_start = time.time()
        for ts in range(50):
            for i, reach in enumerate(reaches_i):
                compute_mc_reach_up2down(
                    reach=reach,
                    connections=rconn_annoated,
                    supernetwork_data=network_data,
                    cols=cols,
                    data=data_values,
                    # ts=ts,
                    flowdepthvel=flowdepthvel
                )
                print(str(i).ljust(20), end='\r')
            if showtiming: print("... in %s seconds." % (time.time() - start_time))
        print("Computation time: ", time.time() - compute_start)
        with np.printoptions(precision=5, suppress=True, linewidth=120):
            print(flowdepthvel)
            print(flowdepthvel.shape)
        #print(sorted(connections.keys()))
    else:
        if verbose: print(f'executing parallel computation on ordered reaches .... ')
        # for terminal_segment, network in networks.items():
        #    print(terminal_segment, network)
        # print(tuple(([x for x in networks.keys()][i], [x for x in networks.values()][i]) for i in range(len(networks))))
        nslist = ([terminal_segment
                      , network
                      , supernetwork_data  # TODO: This should probably be global...
                      , False
                      , debuglevel]
                  for terminal_segment, network in networks.items())
        with multiprocessing.Pool() as pool:
            results = pool.starmap(compute_network, nslist)

    if verbose: print('ordered reach computation complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

if __name__ == '__main__':
    main()

