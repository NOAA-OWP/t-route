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
if ENV_IS_CL:
    root = "/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/"
elif not ENV_IS_CL:
    root = os.path.dirname(os.path.dirname((os.path.abspath(""))))

sys.setrecursionlimit(4000)
sys.path.append(os.path.join(root, r"src", r"python_framework_v01"))
sys.path.append(os.path.join(root, r"src", r"python_routing_v01"))
fortran_source_dir = os.path.join(
    root, r"src", r"fortran_routing", r"mc_pylink_v00", r"MC_singleSeg_singleTS"
)
sys.path.append(fortran_source_dir)

## Muskingum Cunge
COMPILE = True
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r"f2py3")
        fortran_compile_call.append(r"-c")
        fortran_compile_call.append(r"varPrecision.f90")
        fortran_compile_call.append(r"MCsingleSegStime_f2py_NOLOOP.f90")
        fortran_compile_call.append(r"-m")
        fortran_compile_call.append(r"mc_sseg_stime_NOLOOP")
        subprocess.run(
            fortran_compile_call,
            cwd=r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        from mc_sseg_stime_NOLOOP import muskingcunge_module as mc
    except Exception as e:
        print(e)
else:
    from mc_sseg_stime_NOLOOP import muskingcunge_module as mc

connections = None
networks = None
flowdepthvel = None
WRITE_OUTPUT = True

## network and reach utilities
import nhd_network_utilities_v01 as nnu
import nhd_reach_utilities as nru
import nhd_network


def writetoFile(file, writeString):
    file.write(writeString)
    file.write("\n")


def compute_network(
    network,
    conn_data,
    supernetwork_data,
    connections,
    flowdepthvel,
    verbose=False,
    debuglevel=0,
):

    ordered_reaches = [None] * (network["maximum_reach_seqorder"] + 1)
    for head_segment, reach in network["reaches"].items():
        seqorder = reach["seqorder"]
        ordered_reaches[seqorder] = reach

    # initialize flowdepthvel dict
    nts = 50  # one timestep
    # nts = 1440 # number fof timestep = 1140 * 60(model timestep) = 86400 = day

    for ts in range(0, nts):
        # print(f'timestep: {ts}\n')

        for reach in reversed(ordered_reaches):
            # print(f'{{{head_segment}}}:{reach}')

            compute_mc_reach_up2down(
                reach["reach_head"],
                reach,
                connections,
                flowdepthvel,
                conn_data,
                supernetwork_data,
                ts=ts,
                verbose=verbose,
                debuglevel=debuglevel,
            )
            # print(f'{head_segment} {flowdepthvel[head_segment]}')


# ### Psuedocode
#


# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
    head_segment,
    reach,
    connections,
    flowdepthvel,
    conn_data,
    supernetwork_data,
    ts=0,
    verbose=False,
    debuglevel=0,
):
    if verbose:
        print(
            f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments_list'])})"
        )

    if WRITE_OUTPUT:
        filename = f"../../test/output/text/{head_segment}_{ts}.csv"
        file = open(filename, "w+")
        writeString = f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments_list'])}  isterminal: {reach['upstream_reaches'] == {0}} )  reach tail: {reach['reach_tail']}  upstream seg : "

    inds = np.searchsorted(conn_data[:, 0], reach["segments_list"][::-1])
    # upstream flow per reach
    qup_tmp = 0
    if (
        reach["upstream_reaches"] != {0} and len(reach["upstream_reaches"]) > 1
    ):  # Not Headwaters
        # for us in reach['upstream_reaches']:
        qup_tmp = flowdepthvel[
            np.searchsorted(conn_data[:, 0], tuple(reach["upstream_reaches"])), 4
        ].sum()
    if WRITE_OUTPUT:
        writetoFile(file, writeString)

    qup = qup_tmp
    quc = qup

    current_segment = reach["reach_head"]

    if WRITE_OUTPUT:
        writeString = (
            writeString
            + f" timestep: {ts} cur : {current_segment}  upstream flow: {qup_tmp}"
        )
        writetoFile(file, writeString)
        writeString = f"  , , , , , , "
        writetoFile(file, writeString)

    dt = 60.0

    flowdepthvel[inds, 7] = 10

    for i in inds:
        flowdepthvel[i, :4] = flowdepthvel[i, 4:]
        qdc, velc, depthc = singlesegment(
            dt=dt,
            qup=qup,
            quc=quc,
            qdp=flowdepthvel[i, 0],
            qlat=flowdepthvel[i, 7],
            dx=conn_data[i, 1],
            bw=conn_data[i, 2],
            tw=conn_data[i, 3],
            twcc=conn_data[i, 7],
            n_manning=conn_data[i, 4],
            n_manning_cc=conn_data[i, 8],
            cs=conn_data[i, 6],
            s0=conn_data[i, 7],
            velp=flowdepthvel[i, 2],
            depthp=flowdepthvel[i, 1],
        )

        # update flowdepthvel
        flowdepthvel[i, 4:7] = qdc, velc, depthc
        qup = quc = qdc

    if WRITE_OUTPUT:
        file.close()

def process_edge(N, RN, flowdepthvel, conn_data):
    _key = conn_data[:,0]
    i = np.searchsorted(_key, N)

    dt = 60.0
    dx = conn_data[i, 1]
    bw = conn_data[i, 2]
    tw = conn_data[i, 3]
    twcc = conn_data[i, 7]
    n_manning = conn_data[i, 4]
    n_manning_cc = conn_data[i, 8]
    cs = conn_data[i, 6]
    s0 = conn_data[i, 7]

    qup = 0.0
    if N in RN and len(RN[N]) > 0:
        qup = flowdepthvel[np.searchsorted(_key, tuple(RN[N])), 4].sum()

    quc = qup

    flowdepthvel[i, 7] = 10
    flowdepthvel[i, :4] = flowdepthvel[i, 4:]

    # Current flow, depth, vel, qlat
    qdp, depthp, velp, qlat = flowdepthvel[i, 4:]

    qdc, velc, depthc = singlesegment(
        dt,
        qup,
        quc,
        qdp,
        qlat,
        dx, bw, tw, twcc, n_manning, n_manning_cc, cs, s0,
        velp,
        depthp
    )

    # update flowdepthvel
    flowdepthvel[i, 4:7] = qdc, depthc, velc

def singlesegment(
    dt,  # dt
    qup,  # qup
    quc,  # quc
    qdp,  # qdp
    qlat,  # ql
    dx,  # dx
    bw,  # bw
    tw,  # tw
    twcc,  # twcc
    n_manning,  #
    n_manning_cc,  # ncc
    cs,  # cs
    s0,  # s0
    velp,  # velocity at previous time step
    depthp,  # depth at previous time step
):

    # call Fortran routine
    return mc.muskingcungenwm(
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
    # return qdc, vel, depth


def main():

    global connections
    global networks
    global flowdepthvel

    verbose = True
    debuglevel = 0
    showtiming = True

    test_folder = os.path.join(root, r"test")
    geo_input_folder = os.path.join(test_folder, r"input", r"geo")

    # TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    supernetwork = "Pocono_TEST1"
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
    supernetwork_data, supernetwork_values = nnu.set_networks(
        supernetwork=supernetwork,
        geo_input_folder=geo_input_folder,
        verbose=False
        # , verbose = verbose
        ,
        debuglevel=debuglevel,
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
    networks = nru.compose_networks(
        supernetwork_values,
        verbose=False
        # , verbose = verbose
        ,
        debuglevel=debuglevel,
        showtiming=showtiming,
    )
    if verbose:
        print("reach organization complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))

    if showtiming:
        start_time = time.time()
    connections = supernetwork_values[0]

    connections, conn_data = separate_data(
        connections,
        {
            "key": supernetwork_data["key_col"],
            "length": supernetwork_data["length_col"],
            "bottomwidth": supernetwork_data["bottomwidth_col"],
            "topwidth": supernetwork_data["topwidth_col"],
            "manningn": supernetwork_data["manningn_col"],
            "ChSlp": supernetwork_data["ChSlp_col"],
            "slope": supernetwork_data["slope_col"],
            "topwidthcc": supernetwork_data["topwidthcc_col"],
            "manningncc": supernetwork_data["manningncc_col"],
        },
    )
    conn_data = conn_data[conn_data[:, 0].argsort()]

    # Index Map:
    # flow_prev, depth_prev, vel_prev, qlat_prev
    # flow_curr, depth_curr, vel_curr, qlat_curr
    flowdepthvel = np.zeros((len(connections), 8))
    RN = dict(nhd_network.reverse_network(connections))

    for ts in range(1440):
        for n in nhd_network.kahn_toposort(connections):
            process_edge(n, RN, flowdepthvel, conn_data)

    with np.printoptions(precision=5, suppress=True, linewidth=120):
        print(flowdepthvel)
    sorted_conns = sorted(connections.keys())
    print(sorted_conns, all(conn_data[:,0] == sorted_conns))

    # parallelcompute = False
    # if not parallelcompute:
    #     if verbose:
    #         print("executing computation on ordered reaches ...")
    #
    #     for terminal_segment, network in networks.items():
    #         compute_network(
    #             network,
    #             conn_data,
    #             supernetwork_data,
    #             connections,
    #             flowdepthvel,
    #             verbose=False,
    #             debuglevel=debuglevel,
    #         )
    #         print(f"{terminal_segment}")
    #         if showtiming:
    #             print("... in %s seconds." % (time.time() - start_time))
    #
    # else:
    #     if verbose:
    #         print(f"executing parallel computation on ordered reaches .... ")
    #     # for terminal_segment, network in networks.items():
    #     #    print(terminal_segment, network)
    #     # print(tuple(([x for x in networks.keys()][i], [x for x in networks.values()][i]) for i in range(len(networks))))
    #     nslist = (
    #         [
    #             network,
    #             conn_data,  # TODO: This should probably be global...
    #             connections,
    #             flowdepthvel,
    #             False,
    #             debuglevel,
    #         ]
    #         for terminal_segment, network in networks.items()
    #     )
    #     with multiprocessing.Pool() as pool:
    #         results = pool.starmap(compute_network, nslist)

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))


def separate_data(c, col_mapping):
    from operator import itemgetter

    _connections = {}
    _data = []

    col_labels, col_items = tuple(zip(*col_mapping.items()))
    col_getter = itemgetter(*col_items)
    for k, v in c.items():
        cv = v.copy()
        _data.append(col_getter(cv.pop("data")))
        _connections[k] = [cv.pop("downstream")]
    return _connections, np.array(_data)


if __name__ == "__main__":
    main()
