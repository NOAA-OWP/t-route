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
    root = os.path.dirname(os.path.dirname(os.path.abspath("")))
    sys.path.append(r"../python_framework")
    sys.setrecursionlimit(4000)

connections = None
networks = None
flowdepthvel = None
WRITE_OUTPUT = False
ASSUME_SHORT_TS = False

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru


if WRITE_OUTPUT:
    # write to file
    def writetoFile(file, writeString):
        file.write(writeString + "\n")
        file.flush()
        os.fsync(file.fileno())


def compute_network(
    terminal_segment=None,
    network=None,
    supernetwork_data=None
    # , connections = None
    ,
    verbose=False,
    debuglevel=0,
):
    global connections
    global flowdepthvel

    ordered_reaches = {}
    for head_segment, reach in network["reaches"].items():
        if reach["seqorder"] not in ordered_reaches:
            ordered_reaches.update(
                {reach["seqorder"]: []}
            )  # TODO: Should this be a set/dictionary?
        ordered_reaches[reach["seqorder"]].append([head_segment, reach])

    # initialize flowdepthvel dict
    nts = 50  # one timestep
    # nts = 1440 # number fof timestep = 1140 * 60(model timestep) = 86400 = day

    for ts in range(0, nts):
        # print(f'timestep: {ts}\n')

        for x in range(network["maximum_reach_seqorder"], -1, -1):
            for head_segment, reach in ordered_reaches[x]:
                # print(f'{{{head_segment}}}:{reach}')

                compute_mc_reach_up2down(
                    head_segment=head_segment,
                    reach=reach
                    # , connections = connections
                    ,
                    supernetwork_data=supernetwork_data,
                    ts=ts,
                    verbose=verbose,
                    debuglevel=debuglevel,
                )
                # print(f'{head_segment} {flowdepthvel[head_segment]}')


# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
    head_segment=None,
    reach=None
    # , connections = None
    ,
    supernetwork_data=None,
    ts=0,
    verbose=False,
    debuglevel=0,
):
    global connections
    global flowdepthvel
    # global network

    # if verbose: print(f"\nreach: {head_segment}")
    # if verbose: print(f"(reach: {reach})")
    # if verbose: print(f"(n_segs: {len(reach['segments'])})")
    if verbose:
        print(
            f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})"
        )

    if WRITE_OUTPUT:
        filename = f"../../test/output/text/{head_segment}_{ts}.csv"
        file = open(filename, "w+")
        writeString = f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])}  isterminal: {reach['upstream_reaches'] == {supernetwork_data['terminal_code']}} )  reach tail: {reach['reach_tail']}  upstream seg : "
    current_segment = reach["reach_head"]
    next_segment = connections[current_segment]["downstream"]
    # print(f'{current_segment}==>{next_segment} conections:{ncomp} timestep:{ts}')
    # input flow to upstream reach of current reach

    # upstream flow per reach
    qup_tmp = 0
    quc_tmp = 0
    # import pdb; pdb.set_trace()
    if reach["upstream_reaches"] == {supernetwork_data["terminal_code"]}:  # Headwaters
        qup_tmp = 0.0  # no flows in a head channel, only laterals
        quc_tmp = qup_tmp  # no flows in a head channel, only laterals
    else:  # Loop over upstream reaches
        for us in connections[reach["reach_head"]]["upstreams"]:
            if WRITE_OUTPUT:
                writeString = writeString + f" {us} "
            qup_tmp += flowdepthvel[us]["flow"]["curr"]
            if ASSUME_SHORT_TS:
                quc_tmp = qup_tmp
            else:
                quc_tmp += flowdepthvel[us]["flow"]["curr"]

    if WRITE_OUTPUT:
        writeString = (
            writeString
            + f" timestep: {ts} cur : {current_segment}  upstream flow: {qup_tmp}"
        )
        writetoFile(file, writeString)
        writeString = f"  , , , , , , "
        writetoFile(file, writeString)

    while True:
        if ASSUME_SHORT_TS:
            quc_tmp = qup_tmp

        qup = qup_tmp
        quc = quc_tmp

        # for now treating as constant per reach
        dt = 60.0
        # import pdb; pdb.set_trace()
        bw = connections[current_segment]["data"][supernetwork_data["bottomwidth_col"]]
        tw = connections[current_segment]["data"][supernetwork_data["topwidth_col"]]
        twcc = connections[current_segment]["data"][supernetwork_data["topwidthcc_col"]]
        dx = connections[current_segment]["data"][supernetwork_data["length_col"]]
        bw = connections[current_segment]["data"][supernetwork_data["bottomwidth_col"]]
        n_manning = connections[current_segment]["data"][
            supernetwork_data["manningn_col"]
        ]
        n_manning_cc = connections[current_segment]["data"][
            supernetwork_data["manningncc_col"]
        ]
        cs = connections[current_segment]["data"][supernetwork_data["ChSlp_col"]]
        s0 = connections[current_segment]["data"][supernetwork_data["slope_col"]]

        # add some flow
        flowdepthvel[current_segment]["qlat"][
            "curr"
        ] = 10.0  # (ts + 1) * 10.0  # lateral flow per segment

        flowdepthvel[current_segment]["flow"]["prev"] = flowdepthvel[current_segment][
            "flow"
        ]["curr"]
        flowdepthvel[current_segment]["depth"]["prev"] = flowdepthvel[current_segment][
            "depth"
        ]["curr"]
        flowdepthvel[current_segment]["vel"]["prev"] = flowdepthvel[current_segment][
            "vel"
        ]["curr"]
        flowdepthvel[current_segment]["qlat"]["prev"] = flowdepthvel[current_segment][
            "qlat"
        ]["curr"]

        # print (f'counter = {i}')
        # if current_segment == 5559368 or i == 100:
        #    import pdb; pdb.set_trace()

        qlat = flowdepthvel[current_segment]["qlat"]["curr"]  # temporary assigned qlat
        qdp = flowdepthvel[current_segment]["flow"]["prev"]  # temporary assigned qd
        velp = flowdepthvel[current_segment]["vel"]["prev"]
        depthp = flowdepthvel[current_segment]["depth"]["prev"]

        if WRITE_OUTPUT:
            # writeString = f'timestep: {ts} parameters : {current_segment}  {dx} {bw} {tw} {n_manning} {cs} {s0} {dt}'
            # writetoFile(file, writeString)
            writeString = (
                f"{current_segment} , {flowdepthvel[current_segment]['flow']['prev']} "
            )
            writeString = (
                writeString + f", {flowdepthvel[current_segment]['depth']['prev']} "
            )
            writeString = (
                writeString + f", {flowdepthvel[current_segment]['vel']['prev']} "
            )
            writeString = (
                writeString + f", {flowdepthvel[current_segment]['qlat']['curr']} "
            )
            writeString = writeString + f", {qup} "
            writeString = writeString + f", {quc}"
            # writetoFile(file, writeString)

        # run M-C model
        qdc, velc, depthc = singlesegment(
            dt=dt,
            qup=qup,
            quc=quc,
            qdp=qdp,
            qlat=qlat,
            dx=dx,
            bw=bw,
            tw=tw,
            twcc=twcc,
            n_manning=n_manning,
            n_manning_cc=n_manning_cc,
            cs=cs,
            s0=s0,
            velp=velp,
            depthp=depthp,
        )
        # print(qdc, velc, depthc)
        # print(qdc_expected, velc_expected, depthc_expected)

        if WRITE_OUTPUT:
            writeString = writeString + f",  {qdc},  {depthc},  {velc} "
            writetoFile(file, writeString)
        flowdepthvel[current_segment]["flow"]["curr"] = qdc
        flowdepthvel[current_segment]["depth"]["curr"] = depthc
        flowdepthvel[current_segment]["vel"]["curr"] = velc

        qup_tmp = qdp
        if ASSUME_SHORT_TS:
            quc_tmp = qdp
        else:
            quc_tmp = qdc

        if current_segment == reach["reach_tail"]:
            if verbose:
                print(f"{current_segment} (tail)")
            break
        if verbose:
            print(f"{current_segment} --> {next_segment}\n")
        current_segment = next_segment
        next_segment = connections[current_segment]["downstream"]
        # end loop initialized the MC vars


def singlesegment(
    dt,  # dt
    qup=None,  # qup
    quc=None,  # quc
    qdp=None,  # qdp
    qlat=None,  # ql
    dx=None,  # dx
    bw=None,  # bw
    tw=None,  # tw
    twcc=None,  # twcc
    n_manning=None,  #
    n_manning_cc=None,  # ncc
    cs=None,  # cs
    s0=None,  # s0
    velp=None,  # velocity at previous time step
    depthp=None,  # depth at previous time step
):

    # global model # Not necessary -- the tf object is apparently global by default.
    # call ML routine
    input_arr = []
    input_arr.append(
        [
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
        ]
    )
    input_arr = np.array(input_arr)
    vel = 0.5
    depth = 0.5
    qdc = model.predict_on_batch((input_arr))
    # Predicting with the more efficient batch function instead of regular model.predict()
    # See: https://github.com/tensorflow/tensorflow/issues/35835
    # See: https://github.com/tensorflow/tensorflow/issues/33009
    # While much faster, but still not nearly fast enough to compete with the raw SingleSegment
    # Muskingum Cunge as we've written it here.
    # TODO: explore methods to utilize a many-at-a-time call or
    # TODO: Consider finding a way to make the model persist in memory
    return qdc[0][0], vel, depth


def main(parallelcompute=True, num_cores=8, sort="natural"):

    global connections
    global networks
    global flowdepthvel

    verbose = True
    debuglevel = 0
    showtiming = True

    if sort is "natural":
        sorted_networks = list(networks.items())
    elif sort is "sort_rev":
        sorted_networks = list(networks.items())
        sorted_networks.sort(
            reverse=True,
            key=lambda networks_tuple: networks_tuple[1]["total_segment_count"],
        )
    elif sort is "sort":
        sorted_networks = list(networks.items())
        sorted_networks.sort(
            reverse=False,
            key=lambda networks_tuple: networks_tuple[1]["total_segment_count"],
        )

    if showtiming:
        compute_start_time = time.time()
    # parallelcompute = True
    if not parallelcompute:
        if verbose:
            print("executing computation on ordered reaches ...")

        for terminal_segment, network in sorted_networks:
            # for terminal_segment, network in sorted_networks.items():
            # for terminal_segment, network in sorted(list(networks.items())
            # , reverse = True
            # , key = lambda item: item[1]['total_segment_count']):
            if showtiming:
                n_start_time = time.time()
            compute_network(
                terminal_segment=terminal_segment,
                network=network,
                supernetwork_data=supernetwork_data
                # , connections = connections
                ,
                verbose=False
                # , verbose = verbose
                ,
                debuglevel=debuglevel,
            )
            if terminal_segment == 22811611:
                print(f"{terminal_segment}")
                if showtiming:
                    print("... in %s seconds." % (time.time() - n_start_time))

    else:
        if showtiming:
            start_time = time.time()
        if verbose:
            print(f"executing parallel computation on ordered reaches .... ")
        # for terminal_segment, network in networks.items():
        #    print(terminal_segment, network)
        # print(tuple(([x for x in networks.keys()][i], [x for x in networks.values()][i]) for i in range(len(networks))))
        nslist = (
            [
                terminal_segment,
                network,
                supernetwork_data,  # TODO: This should probably be global...
                False,
                debuglevel,
            ]
            for terminal_segment, network in sorted_networks
        )
        # for terminal_segment, network in networks.items())
        # for terminal_segment, network in sorted(list(networks.items())
        # , reverse = False
        # , key = lambda item: item[1]['total_segment_count']))
        if showtiming:
            p_start_time = time.time()
        with multiprocessing.Pool(num_cores) as pool:
            results = pool.starmap(compute_network, nslist)
        if showtiming:
            print("... in %s seconds." % (time.time() - p_start_time))

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - compute_start_time))


if __name__ == "__main__":

    # global connections
    # global networks
    # global flowdepthvel

    verbose = True
    debuglevel = -2
    showtiming = True

    test_folder = os.path.join(root, r"test")
    geo_input_folder = os.path.join(test_folder, r"input", r"geo")

    # TODO: Make these commandline args
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_Named_Streams'
    # supernetwork = 'Pocono_TEST1'
    """##NHD CONUS order 5 and greater"""
    # supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    supernetwork = "Mainstems_CONUS"
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

    import tensorflow as tf

    # Tensorflow is messy -- it needs to be loaded after reading the data
    # or it will mess up the libraries used by xarray to read the netcdf
    model = tf.keras.models.load_model("ML_MC_PRES6")

    connections = supernetwork_values[0]

    flowdepthvel = {
        connection: {
            "flow": {"prev": 0, "curr": 0},
            "depth": {"prev": 0, "curr": 0},
            "vel": {"prev": 0, "curr": 0},
            "qlat": {"prev": 0, "curr": 0},
        }
        for connection in connections
    }

    parallels = [False, True]
    sorts = ["sort_rev", "sort", "natural"]
    cores = range(1, 21)
    for sort in sorts:
        for parallel in parallels:
            if parallel:
                for core in cores:
                    print(f"parallel: {parallel} core:{core} sort:{sort}")
                    main(parallel, core, sort)
            else:
                print(f"parallel: {parallel} core:SERIAL sort:{sort}")
                main(parallel, None, sort)
