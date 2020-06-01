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
    fortran_compile_call = []
    fortran_compile_call.append(r'f2py3')
    fortran_compile_call.append(r'-c')
    fortran_compile_call.append(r'MCsingleSegStime_f2py_NOLOOP.f90')
    fortran_compile_call.append(r'-m')
    fortran_compile_call.append(r'mc_sseg_stime_NOLOOP')
    subprocess.run(fortran_compile_call, cwd=r'../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
import mc_sseg_stime_NOLOOP as mc

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


def compute_network(
        reach,
        supernetwork_data,
        connections,
        data,
        flowdepthvel
):

    # = {connection:{'flow':{'prev':-999, 'curr':-999}
    #                            , 'depth':{'prev':-999, 'curr':-999}
    #                            , 'vel':{'prev':-999, 'curr':-999}} for connection in connections}

    # print(tuple(([x for x in network.keys()][i], [x for x in network.values()][i]) for i in range(len(network))))

    # if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_order']}")

    #ordered_reaches = {}
    #for head_segment, reach in network['reaches'].items():
    #    if reach['seqorder'] not in ordered_reaches:
    #        ordered_reaches.update({reach['seqorder']: []})  # TODO: Should this be a set/dictionary?
    #    ordered_reaches[reach['seqorder']].append([head_segment
    #                                                  , reach
    #                                               ])

    # initialize flowdepthvel dict
    #nts = 50  # one timestep
    #nts = 1440 # number fof timestep = 1140 * 60(model timestep) = 86400 = day

    reversed_conns = nhd_network.reverse_network(connections)
    #for ts in range(nts):
    compute_mc_reach_up2down(
        reach=reach,
        connections=reversed_conns,
        supernetwork_data=supernetwork_data,
        data=data,
        #ts=ts,
        flowdepthvel=flowdepthvel
    )

# ### Psuedocode
# 



# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
        reach=None,
        connections=None,
        supernetwork_data=None,
        data=None,
        ts=0,
        flowdepthvel=None,
        verbose=False,
        debuglevel=0,
        write_output=False,
        assume_short_ts=False,
):

    # if verbose: print(f"\nreach: {head_segment}")
    # if verbose: print(f"(reach: {reach})")
    # if verbose: print(f"(n_segs: {len(reach['segments'])})")

    #if write_output:
    #    filename = f"../../test/output/text/{head_segment}_{ts}.csv"
    #    file = open(filename, "w+")
    #    writeString = f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])}  isterminal: {reach['upstream_reaches'] == {supernetwork_data['terminal_code']}} )  reach tail: {reach['reach_tail']}  upstream seg : "

        # upstream flow per reach
    qup = 0.0
    quc = 0.0
    # import pdb; pdb.set_trace()
    for us in connections.get(reach[0], ()):
        if write_output:
            writeString = writeString + f"\n upstream seg : {us}"
        i = data.index.get_loc(us)
        qup += flowdepthvel[i, 0]
        quc += flowdepthvel[i, 4]
    #if write_output:
    #    writetoFile(file, writeString)

    #if write_output:
    #    writeString = (
    #            writeString
    #            + f" timestep: {ts} cur : {current_segment}  upstream flow: {qup}"
    #    )
    #    writetoFile(file, writeString)
    #    writeString = f"  , , , , , , "
    #    writetoFile(file, writeString)

    write_buffer = []


    for current_segment in reach:
        if current_segment not in data.index:
            break
        i = data.index.get_loc(current_segment)

        # for now treating as constant per reach
        dt = 60.0
        bw = data.loc[current_segment, supernetwork_data["bottomwidth_col"]]
        tw = data.loc[current_segment, supernetwork_data["topwidth_col"]]
        twcc = data.loc[current_segment, supernetwork_data["topwidthcc_col"]]
        dx = data.loc[current_segment, supernetwork_data["length_col"]]
        n_manning = data.loc[current_segment, supernetwork_data["manningn_col"]]
        n_manning_cc = data.loc[current_segment, supernetwork_data["manningncc_col"]]
        cs = data.loc[current_segment, supernetwork_data["ChSlp_col"]]
        s0 = data.loc[current_segment, supernetwork_data["slope_col"]]

        # add some flow
        flowdepthvel[i, 7] = qlat = 10.0
        #current_flow["qlat"][
        #    "curr"
        #] = qlat = 10.0  # (ts + 1) * 10.0  # lateral flow per segment

        qdp, depthp, velp = flowdepthvel[i, 0:3]

        flowdepthvel[i, :4] = flowdepthvel[i, 4:]

        if assume_short_ts:
            quc = qup

        # run M-C model

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
            depthp,
        )
        #print(f"ts={ts}", qdc, velc, depthc)
        # print(qdc_expected, velc_expected, depthc_expected)

        if write_output:
            write_buffer.append(
                ",".join(
                    map(
                        str,
                        (
                            current_segment,
                            qdp,
                            depthp,
                            velp,
                            qlat,
                            qup,
                            quc,
                            qdc,
                            depthc,
                            velc,
                        ),
                    )
                )
            )

        # for next segment qup / quc use the previous flow values
        flowdepthvel[i, 4:7] = qdc, depthc, velc

        quc = qdc
        qup = qdp

        #with np.printoptions(precision=5, suppress=True, linewidth=120):
        #    print('-'*40, f"ts={ts}", '-'*40)
        #    print(flowdepthvel)

        # end loop initialized the MC vars
    if write_output:
        writetoFile(file, "\n".join(write_buffer))
        file.close()


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
    data = nhd_io.read(network_data['geo_file_path'])
    data = data.set_index(network_data['key_col'])
    if 'mask_file_path' in network_data:
        data_mask = nhd_io.read_mask(network_data['mask_file_path'],
                                     layer_string=network_data['mask_layer_string'])
        data = data.filter(data_mask.iloc[:, network_data["mask_key_col"]], axis=0)

    connections = nhd_network.extract_network(data, network_data["downstream_col"])

    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    # STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into reaches ...')
    reaches = nhd_network.dfs_decomposition(connections)

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

    flowdepthvel = np.zeros((len(data), 8))
    #flowdepthvel = {node: {'flow': {'prev': 0, 'curr': 0}
    #    , 'depth': {'prev': 0, 'curr': 0}
    #    , 'vel': {'prev': 0, 'curr': 0}
    #    , 'qlat': {'prev': 0, 'curr': 0}} for node in nhd_network.nodes(connections)}

    parallelcompute = False
    if not parallelcompute:
        if verbose: print('executing computation on ordered reaches ...')

        compute_start = time.time()
        for ts in range(50):
            for i, reach in enumerate(reaches):
                compute_network(
                    reach,
                    network_data,
                    connections,
                    data,
                    flowdepthvel
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

