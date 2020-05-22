#!/usr/bin/env python
# coding: utf-8
# example usage: python compute_nhd_routing_SingleSeg.py -v -t -w -n Mainstems_CONUS


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
import argparse


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
    )
    parser.add_argument(
        "-v, --verbose",
        help="Verbose output (leave blank for quiet output)",
        dest="verbose",
        action="store_true",
    )
    parser.add_argument(
        "--assume_short_ts",
        help="Use the previous timestep value for upstream flow",
        dest="assume_short_ts",
        action="store_true",
    )
    parser.add_argument(
        "-o, --write_output",
        help="Write output files (leave blank for no writing)",
        dest="write_output",
        action="store_true",
    )
    parser.add_argument(
        "-t, --showtiming",
        help="Set the showtiming (leave blank for no timing information)",
        dest="showtiming",
        action="store_true",
    )
    parser.add_argument(
        "-w, --break_at_waterbodies",
        help="Use the waterbodies in the route-link dataset to divide the computation (leave blank for no splitting)",
        dest="break_network_at_waterbodies",
        action="store_true",
    )
    parser.add_argument(
        "-n, --supernetwork",
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

    return parser.parse_args()


ENV_IS_CL = False
if ENV_IS_CL:
    root = "/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/"
elif not ENV_IS_CL:
    root = os.path.dirname(os.path.dirname(os.path.abspath("")))
    sys.path.append(r"../python_framework")
    sys.path.append(r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS")
    sys.setrecursionlimit(4000)

## Muskingum Cunge
COMPILE = True
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r"f2py3")
        fortran_compile_call.append(r"-c")
        fortran_compile_call.append(r"MCsingleSegStime_f2py_NOLOOP.f90")
        fortran_compile_call.append(r"-m")
        fortran_compile_call.append(r"mc_sseg_stime")
        subprocess.run(
            fortran_compile_call,
            cwd=r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        import mc_sseg_stime as mc
    except Exception as e:
        print(e)
else:
    import mc_sseg_stime as mc

connections = None
networks = None
flowdepthvel = None

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru


def writetoFile(file, writeString):
    file.write(writeString)
    file.write("\n")


def compute_network(
    terminal_segment=None,
    network=None,
    supernetwork_data=None,
    verbose=False,
    debuglevel=0,
    write_output=False,
    assume_short_ts=False,
):
    global connections
    global flowdepthvel

    # = {connection:{'flow':{'prev':-999, 'curr':-999}
    #                            , 'depth':{'prev':-999, 'curr':-999}
    #                            , 'vel':{'prev':-999, 'curr':-999}} for connection in connections}

    # print(tuple(([x for x in network.keys()][i], [x for x in network.values()][i]) for i in range(len(network))))

    # if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_order']}")

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
                    reach=reach,
                    supernetwork_data=supernetwork_data,
                    ts=ts,
                    verbose=verbose,
                    debuglevel=debuglevel,
                    write_output=write_output,
                    assume_short_ts=assume_short_ts,
                )
                # print(f'{head_segment} {flowdepthvel[head_segment]}')


# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
    head_segment=None,
    reach=None,
    supernetwork_data=None,
    ts=0,
    verbose=False,
    debuglevel=0,
    write_output=False,
    assume_short_ts=False,
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

    if write_output:
        filename = f"../../test/output/text/{head_segment}_{ts}.csv"
        file = open(filename, "w+")
        writeString = f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])}  isterminal: {reach['upstream_reaches'] == {supernetwork_data['terminal_code']}} )  reach tail: {reach['reach_tail']}  upstream seg : "

    # upstream flow per reach
    qup_tmp = 0
    quc_tmp = 0
    # import pdb; pdb.set_trace()
    if reach["upstream_reaches"] == {supernetwork_data["terminal_code"]}:  # Headwaters
        qup_tmp = 0.0  # no flows in a head channel, only laterals
        quc_tmp = qup_tmp  # no flows in a head channel, only laterals
    else:  # Loop over upstream reaches
        for us in connections[reach["reach_head"]]["upstreams"]:
            if write_output:
                writeString = writeString + f"\n upstream seg : {us}"
            qup_tmp += flowdepthvel[us]["flow"]["curr"]
            if assume_short_ts:
                quc_tmp = qup_tmp
            else:
                quc_tmp += flowdepthvel[us]["flow"]["curr"]
    if write_output:
        writetoFile(file, writeString)

    current_segment = reach["reach_head"]
    next_segment = connections[current_segment]["downstream"]

    if write_output:
        writeString = (
            writeString
            + f" timestep: {ts} cur : {current_segment}  upstream flow: {qup_tmp}"
        )
        writetoFile(file, writeString)
        writeString = f"  , , , , , , "
        writetoFile(file, writeString)

    while True:
        if assume_short_ts:
            quc_tmp = qup_tmp

        data = connections[current_segment]["data"]
        current_flow = flowdepthvel[current_segment]

        qup = qup_tmp
        quc = quc_tmp

        # for now treating as constant per reach
        dt = 60.0
        # import pdb; pdb.set_trace()
        bw = data[supernetwork_data["bottomwidth_col"]]
        tw = data[supernetwork_data["topwidth_col"]]
        twcc = data[supernetwork_data["topwidthcc_col"]]
        dx = data[supernetwork_data["length_col"]]
        bw = data[supernetwork_data["bottomwidth_col"]]
        n_manning = data[supernetwork_data["manningn_col"]]
        n_manning_cc = data[supernetwork_data["manningncc_col"]]
        cs = data[supernetwork_data["ChSlp_col"]]
        s0 = data[supernetwork_data["slope_col"]]

        # add some flow
        current_flow["qlat"][
            "curr"
        ] = 10.0  # (ts + 1) * 10.0  # lateral flow per segment

        current_flow["flow"]["prev"] = qdp = current_flow["flow"]["curr"]
        current_flow["depth"]["prev"] = depthp = current_flow["depth"]["curr"]
        current_flow["vel"]["prev"] = velp = current_flow["vel"]["curr"]
        current_flow["qlat"]["prev"] = qlat = current_flow["qlat"]["curr"]

        # print (f'counter = {i}')
        # if current_segment == 5559368 or i == 100:
        #    import pdb; pdb.set_trace()

        if write_output:
            # writeString = f'timestep: {ts} parameters : {current_segment}  {dx} {bw} {tw} {n_manning} {cs} {s0} {dt}'
            # writetoFile(file, writeString)
            writeString = f"{current_segment} , {qdp} "
            writeString = writeString + f", {depthp} "
            writeString = writeString + f", {velp} "
            writeString = writeString + f", {qlat} "
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

        if write_output:
            writeString = writeString + f",  {qdc},  {depthc},  {velc} "
            writetoFile(file, writeString)
        # for next segment qup / quc use the previous flow values
        current_flow["flow"]["curr"] = qdc
        current_flow["depth"]["curr"] = depthc
        current_flow["vel"]["curr"] = velc

        qup_tmp = qdp
        if assume_short_ts:
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
    if write_output:
        file.close()


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

    args = _handle_args()

    global connections
    global networks
    global flowdepthvel

    debuglevel = -1 * int(args.debuglevel)
    verbose = args.verbose
    showtiming = args.showtiming
    supernetwork = args.supernetwork
    break_network_at_waterbodies = args.break_network_at_waterbodies
    write_output = args.write_output
    assume_short_ts = args.assume_short_ts

    test_folder = os.path.join(root, r"test")
    geo_input_folder = os.path.join(test_folder, r"input", r"geo")

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
        break_network_at_waterbodies=break_network_at_waterbodies,
        verbose=False,
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

    flowdepthvel = {
        connection: {
            "flow": {"prev": 0, "curr": 0},
            "depth": {"prev": 0, "curr": 0},
            "vel": {"prev": 0, "curr": 0},
            "qlat": {"prev": 0, "curr": 0},
        }
        for connection in connections
    }

    parallelcompute = True
    if not parallelcompute:
        if verbose:
            print("executing computation on ordered reaches ...")

        for terminal_segment, network in networks.items():
            compute_network(
                terminal_segment=terminal_segment,
                network=network,
                supernetwork_data=supernetwork_data,
                verbose=False,
                debuglevel=debuglevel,
                write_output=write_output,
                assume_short_ts=assume_short_ts,
            )
            print(f"{terminal_segment}")
            if showtiming:
                print("... in %s seconds." % (time.time() - start_time))

    else:
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
                write_output,
                assume_short_ts,
            ]
            for terminal_segment, network in networks.items()
        )
        with multiprocessing.Pool() as pool:
            results = pool.starmap(compute_network, nslist)

    if verbose:
        print("ordered reach computation complete")
    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))


if __name__ == "__main__":
    main()
