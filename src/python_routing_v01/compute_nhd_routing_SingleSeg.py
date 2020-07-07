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
import csv
import netCDF4
import numpy as np
import argparse
from datetime import datetime

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
        "-v",
        "--verbose",
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
        "-o",
        "--write_output",
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
        "--break_at_waterbodies",
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

    return parser.parse_args()


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
        fortran_compile_call.append(r"mc_sseg_stime")
        subprocess.run(
            fortran_compile_call,
            cwd=r"../fortran_routing/mc_pylink_v00/MC_singleSeg_singleTS",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        from mc_sseg_stime import muskingcunge_module as mc
    except Exception as e:
        print(e)
else:
    from mc_sseg_stime import muskingcunge_module as mc

connections = None
networks = None
flowdepthvel = None
#WRITE_OUTPUT = False  # True

## network and reach utilities
import nhd_network_utilities_v01 as nnu
import nhd_reach_utilities as nru


def writetoFile(file, writeString):
    file.write(writeString + '\n')
    file.flush()
    os.fsync(file.fileno())



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


    # print(tuple(([x for x in network.keys()][i], [x for x in network.values()][i]) for i in range(len(network))))

    # if verbose: print(f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_order']}")

    ordered_reaches = {}
    for head_segment, reach in network["reaches"].items():
        if reach["seqorder"] not in ordered_reaches:
            ordered_reaches.update(
                {reach["seqorder"]: []}
            )  
        ordered_reaches[reach["seqorder"]].append([head_segment, reach])

    # initialize flowdepthvel dict
    # nts = 50  # one timestep
    nts = 143  # test with dt =10
    dt = 300  # in seconds
    # nts = 1440 # number of  timestep = 1140 * 60(model timestep) = 86400 = day

    # initialize write to files variable
    writeToCSV = True
    writeToNETCDF = True
    pathToOutputFile = os.path.join(root, "test", "output", "text") 

    for ts in range(0, nts):
        for x in range(network["maximum_reach_seqorder"], -1, -1):
            for head_segment, reach in ordered_reaches[x]:

                compute_mc_reach_up2down(
                    head_segment=head_segment,
                    reach=reach,
                    supernetwork_data=supernetwork_data,
                    ts=ts,
                    dt=dt, 
                    verbose=verbose,
                    debuglevel=debuglevel,
                    write_output=write_output,
                    assume_short_ts=assume_short_ts,
                )

    if (writeToCSV):
        for x in range(network['maximum_reach_seqorder'], -1, -1):
            for head_segment, reach in ordered_reaches[x]:
                printarray(reach=reach
                           , verbose=verbose
                           , debuglevel=debuglevel
                           , pathToOutputFile=pathToOutputFile
                           )

    if (writeToNETCDF):
        writeArraytoNC(network=network
                       , nts=nts
                       , dt=dt
                       , verbose=verbose
                       , debuglevel=debuglevel
                       , pathToOutputFile=pathToOutputFile
                       )



# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
    head_segment=None,
    reach=None,
    supernetwork_data=None,
    ts=0,
    dt=60, 
    verbose=False,
    debuglevel=0,
    write_output=False,
    assume_short_ts=False,
):
    global connections
    global flowdepthvel

    if verbose: print(f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})")

    # upstream flow per reach
    qup = 0.0
    quc = 0.0
    # import pdb; pdb.set_trace()
    if reach["upstream_reaches"] != {
        supernetwork_data["terminal_code"]
    }:  # Not Headwaters
        for us in connections[reach["reach_head"]]["upstreams"]:
            qup += flowdepthvel[us]['flowval'][-1]
            quc += flowdepthvel[us]['flowval'][0]

    if assume_short_ts:
        quc = qup

    current_segment = reach["reach_head"]
    # next_segment = connections[current_segment]["downstream"]


    while True:
        data = connections[current_segment]["data"]
        #current_flow = flowdepthvel[current_segment]

        # for now treating as constant per reach
        bw = data[supernetwork_data["bottomwidth_col"]]
        tw = 0.01 * bw # data[supernetwork_data["topwidth_col"]]
        twcc = tw # data[supernetwork_data["topwidthcc_col"]]
        dx = data[supernetwork_data["length_col"]]
        bw = data[supernetwork_data["bottomwidth_col"]]
        n_manning = data[supernetwork_data["manningn_col"]]
        n_manning_cc = n_manning # data[supernetwork_data["manningncc_col"]]
        cs = data[supernetwork_data["ChSlp_col"]]
        s0 = data[supernetwork_data["slope_col"]]

        qlat = flowdepthvel[current_segment]['qlatval'][ts]


        if ts > 0:
            qdp = flowdepthvel[current_segment]['flowval'][-1]
            velp = flowdepthvel[current_segment]['velval'][-1]
            depthp = flowdepthvel[current_segment]['depthval'][-1]
        else:
            qdp = 0
            velp = 0
            depthp = 0


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


        # for next segment qup / quc use the previous flow values
        if ts > 0:
            qup = flowdepthvel[current_segment]['flowval'][-1]  # input for next segment
        else:
            qup = 0

        quc = qdc  # input for next segment
        if assume_short_ts:
            quc = qup

        flowdepthvel[current_segment]['flowval'].append(qdc)
        flowdepthvel[current_segment]['depthval'].append(depthc)
        flowdepthvel[current_segment]['velval'].append(velc)
        flowdepthvel[current_segment]['time'].append(ts * dt)

        if current_segment == reach["reach_tail"]:
            if verbose:
                print(f"{current_segment} (tail)")
            break
        if verbose:
            print(f"{current_segment} --> {next_segment}\n")
        #current_segment = next_segment
        #next_segment = connections[current_segment]["downstream"]
        next_segment = connections[current_segment]['downstream']
        current_segment = next_segment 
       # end loop initialized the MC vars
    # end while loop


# ### Psuedocode
# Write Array  to CSV file
# arguments reach , pathToOutputFile
# using global connections and flowdepthvel.
def printarray(reach=None
               , verbose=False
               , debuglevel=0
               , pathToOutputFile="../../test/output/text"
               ):
    global connections
    global flowdepthvel

    # define CSV file Header
    header = [['time', 'qlat', 'q', 'd', 'v']]

    # Loop over reach segments
    current_segment = reach['reach_head']
    next_segment = connections[current_segment]['downstream']

    while True:
        filename = f'{pathToOutputFile}/{current_segment}.csv'  #
        if verbose: print(f'printing to --> {filename} \n')
        with open(filename, 'w+') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_ALL)
            csvwriter.writerows(header)
            csvwriter.writerows(zip(flowdepthvel[current_segment]['time'],
                                    flowdepthvel[current_segment]['qlatval'],
                                    flowdepthvel[current_segment]['flowval'],
                                    flowdepthvel[current_segment]['depthval'],
                                    flowdepthvel[current_segment]['velval']))

        if current_segment == reach['reach_tail']:
            if verbose: print(f'{current_segment} (tail)')
            break
        if verbose: print(f'{current_segment} --> {next_segment}\n')
        current_segment = next_segment
        next_segment = connections[current_segment]['downstream']


# ### Psuedocode
# Write Array  to Arrays for Netcdf file and then call writeNC function to write output data to NC file
# arguments network, number of timesteps (nts),  timestep in seconds  (dt),  pathToOutputFile
# using global connections and flowdepthvel.
def writeArraytoNC(network=None
                   , nts=0
                   , dt=60
                   , verbose=False
                   , debuglevel=0
                   , pathToOutputFile="../../test/output/text/"
                   ):
    global connections
    global flowdepthvel
    # create  array variables to copy from python "flowdepthvel" which is global
    flowdepthvel_data = {'segment': []
        , 'time': []
        , 'qlatval': []
        , 'flowval': []
        , 'depthval': []
        , 'velval': []}

    ordered_reaches = {}
    for head_segment, reach in network['reaches'].items():
        if reach['seqorder'] not in ordered_reaches:
            ordered_reaches.update({reach['seqorder']: []})  # TODO: Should this be a set/dictionary?
        ordered_reaches[reach['seqorder']].append([head_segment, reach])

    # get data into array - preparation step
    TIME_WRITTEN = False
    write_segment = None
    for x in range(network['maximum_reach_seqorder'], -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            current_segment = reach['reach_head']
            while True:
                # appending data from each segments to a single list  "flowdepthvel_data"
                # preserving ordering same as segment in a reach
                flowdepthvel_data['qlatval'].append(flowdepthvel[current_segment]['qlatval'])
                flowdepthvel_data['flowval'].append(flowdepthvel[current_segment]['flowval'])
                flowdepthvel_data['depthval'].append(flowdepthvel[current_segment]['depthval'])
                flowdepthvel_data['velval'].append(flowdepthvel[current_segment]['velval'])
                # write segment flowdepthvel_data['segment']
                flowdepthvel_data['segment'].append(current_segment)
                if not TIME_WRITTEN:
                    # write time only once - for first segment
                    flowdepthvel_data['time'].append(flowdepthvel[current_segment]['time'])
                    TIME_WRITTEN = True

                if current_segment == reach['reach_tail']:
                    write_segment = current_segment
                    if verbose: print(f'{current_segment} (tail)')
                    break
                next_segment = connections[current_segment]['downstream']
                if verbose: print(f'{current_segment} --> {next_segment}\n')
                current_segment = next_segment

    # check number of timesteps should match the time the data is written
    if (int(len(flowdepthvel_data['time'][0])) != int(nts)):
        print(f"Number of timesteps  {nts} does not match data timesteps {len(flowdepthvel_data['time'][0])}\n")
        return


    writeNC(flowdepthvel_data=flowdepthvel_data
            , segment_count=network['total_segment_count']
            , terminal_segment=write_segment
            , nts=nts
            , dt=dt
            , pathToOutputFile=pathToOutputFile
            , verbose=verbose
            , debuglevel=debuglevel)

# ### Psuedocode
# Write  netcdf  file
# arguments flowdepthvel , segment_count, terminal_segment (for filename and as identifier of reach)
#           number of timesteps (nts),  timestep in seconds  (dt),  verbose , debuglevel , pathToOutputFile
def writeNC(flowdepthvel_data=None
            , segment_count=0
            , terminal_segment=None
            , nts=0
            , dt=60
            , pathToOutputFile="../../test/output/text"
            , verbose=False
            , debuglevel=0
            ):
    # start writing data to nc file
    filename = f"{pathToOutputFile}/{terminal_segment}.nc"  # ncfile'
    ncfile = netCDF4.Dataset(filename, mode="w", format="NETCDF4")
    # segcount = total segments for the current reach
    segcount = ncfile.createDimension('stations', segment_count)  # segment
    # timecount = number of timesteps
    timecount = ncfile.createDimension('time', nts)  # unlimited axis (can be appended to).
    # analysis time =  current time , hence count =1
    analysistimecount = ncfile.createDimension('analysistime', 1)  # unlimited axis (can be appended to).
    ncfile.title = f'Result of MC for Reach with terminal segment {terminal_segment}'
    ncfile.subtitle = "MC Python Module"
    ncfile.anything = f"streamflow , depth, velocity and lateral flow for Reach with terminal segment {terminal_segment}"
    # write streamflow
    flow = ncfile.createVariable('flow', np.float64, ('time', 'stations'))  # note: unlimited dimension is leftmost
    flow[:, :] = np.transpose(np.array(flowdepthvel_data['flowval'], dtype=float))
    flow.units = 'cu ft/s'
    flow.standard_name = 'streamflow'  # this is a CF standard name
    # write depth
    depth = ncfile.createVariable('depth', np.float64,
                                  ('time', 'stations'))  # note: unlimited dimension is leftmost
    depth[:, :] = np.transpose(np.array(flowdepthvel_data['depthval'], dtype=float))
    depth.units = 'ft'  #
    depth.standard_name = 'depth'  # this is a CF standard name
    # write velocity
    velocity = ncfile.createVariable('velocity', np.float64,
                                     ('time', 'stations'))  # note: unlimited dimension is leftmost
    velocity.units = 'ft/s'  #
    velocity[:, :] = np.transpose(np.array(flowdepthvel_data['velval'], dtype=float))
    velocity.standard_name = 'velocity'  # this is a CF standard name
    # write  lateral flow (input from NWM)
    lateralflow = ncfile.createVariable('lateralflow', np.float64,
                                        ('time', 'stations'))  # note: unlimited dimension is leftmost
    lateralflow[:, :] = np.transpose(np.array(flowdepthvel_data['qlatval'], dtype=float))
    lateralflow.units = 'cu ft/s'  #
    lateralflow.standard_name = 'lateralflow'  # this is a CF standard name
    # write time in seconds since  TODO get time from lateral flow from NWM
    time = ncfile.createVariable('time', np.float64, 'time')
    time.units = 'seconds since 2011-08-27 00:00:00'  ## TODO get time fron NWM as argument to this function
    time.long_name = 'time'
    time[:] = flowdepthvel_data['time']
    # write segment ids
    segment = ncfile.createVariable('station_id', np.int, ('stations'))
    segment.long_name = 'feature id'
    segment[:] = flowdepthvel_data['segment']
    # write analysis time = current time = time MC model is run
    analysistime = ncfile.createVariable('analysis_time', np.double, 'analysistime')
    analysistime.long_name = 'forecast_reference_time'
    analysistime.units = ' minutes since ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    analysistime[:] = 0
    #
    ncfile.close()
    if verbose: print(f'{filename} closed!')

## call to singlesegment MC Fortran Module
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


# Main Routine
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
    # supernetwork = 'Pocono_TEST2'
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
        start_time = time.time()

    connections = supernetwork_values[0]

    flowdepthvel = {connection: {'qlatval': []
        , 'time': []
        , 'flowval': []
        , 'depthval': []
        , 'velval': []} for connection in connections}


    # Lateral flow
    ## test 1. Take lateral flow from wrf-hydro output from Pocono Basin
    qlcol = 54
    qlrow = 144
    ql = np.zeros((qlrow, qlcol))
    ql_input_folder = os.path.join(root, r'./test/input/text/Pocono_ql_testsamp1_nwm_mc.txt')
    for j in range(0, qlcol):
        ql[0, j] = int(np.loadtxt(ql_input_folder, max_rows=1, usecols=(j + 2)))
        ql[1:, j] = np.loadtxt(ql_input_folder, skiprows=2, usecols=(j + 2))
    for j in range(0, qlcol):
        flowdepthvel[int(ql[0, j])]['qlatval'] = ql[1:, j].tolist()

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
