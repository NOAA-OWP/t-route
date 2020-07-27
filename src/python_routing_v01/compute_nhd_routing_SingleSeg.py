# command line args and usage with: python compute_nhd_routing_SingleSeg.py --help
# example to run test: python compute_nhd_routing_SingleSeg.py -v --test
# example usage: python compute_nhd_routing_SingleSeg.py -v -t -w -onc -n Mainstems_CONUS

# a notebook-based version of very similar code is found here:
# https://github.com/NOAA-OWP/t-route/blob/master/notebooks/compute_nhd_routing_v2_clean_with_lateral_inflow_data.ipynb

## Parallel execution
import os
import sys
import time
import traceback
import numpy as np
import argparse
import pathlib
import pandas as pd
import netCDF4
import csv
from datetime import datetime
import multiprocessing
import sys


def _handle_args():
    # TODO: Convert to global argparser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--debuglevel",
        help="Set the debuglevel",
        dest="debuglevel",
        type=int,
        choices=[0, 1, 2, 3],
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
        "--test" "--run-pocono-test-example",
        help="Use the data values stored in the repository for a test of the Pocono network",
        dest="run_pocono_test",
        action="store_true",
    )
    parser.add_argument(
        "--sts" "--assume-short-ts",
        help="Use the previous timestep value for upstream flow",
        dest="assume_short_ts",
        action="store_true",
    )
    parser.add_argument(
        "-ocsv",
        "--write-output-csv",
        help="Write csv output files (omit flag for no csv writing)",
        dest="write_csv_output",
        action="store_true",
    )
    parser.add_argument(
        "-onc",
        "--write-output-nc",
        help="Write netcdf output files (omit flag for no netcdf writing)",
        dest="write_nc_output",
        action="store_true",
    )
    parser.add_argument(
        "--dt",
        "--time-step",
        help="Set the default timestep length",
        dest="dt",
        default=300,
    )
    parser.add_argument(
        "--nts",
        "--number-of-timesteps",
        help="Set the number of timesteps to execute",
        dest="nts",
        default=144,
    )
    parser.add_argument(
        "--ql",
        "--constant_qlateral",
        help="Set the number of timesteps to execute",
        dest="qlat_const",
        default=10,
    )
    parser.add_argument(
        "-t",
        "--showtiming",
        help="Set the showtiming (omit flag for no timing information)",
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
    parser.add_argument(
        "--parallel",
        help="Use the parallel computation engine (omit flag for serial computation)",
        dest="parallel_compute",
        action="store_true",
    )

    return parser.parse_args()


root = pathlib.Path("../..").resolve()

sys.setrecursionlimit(4000)
sys.path.append(os.path.join(root, r"src", r"python_framework_v01"))
sys.path.append(os.path.join(root, r"src", r"python_routing_v01"))
fortran_routing_dir = os.path.join(
    root, r"src", r"fortran_routing", r"mc_pylink_v00", r"MC_singleSeg_singleTS"
)
fortran_reservoir_dir = os.path.join(
    root, r"src", r"fortran_routing", r"mc_pylink_v00", r"Reservoir_singleTS"
)
sys.path.append(fortran_routing_dir)

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
            cwd=fortran_routing_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        from mc_sseg_stime import muskingcunge_module as mc
    except Exception as e:
        print(e)
        traceback.print_exc()
else:
    from mc_sseg_stime import muskingcunge_module as mc

connections = None
networks = None
flowveldepth = None
waterbodies_df = None

# WRITE_OUTPUT = False  # True

## network and reach utilities
import nhd_network_utilities_v01 as nnu
import nhd_reach_utilities as nru

import pdb; pdb.set_trace()
sys.path.append(fortran_reservoir_dir)
if COMPILE:
    try:
        import subprocess

        fortran_compile_call = []
        fortran_compile_call.append(r"f2py3")
        fortran_compile_call.append(r"-c")
        fortran_compile_call.append(r"module_levelpool.f90")
        # fortran_compile_call.append(r"varPrecision.f90")
        # fortran_compile_call.append(r"Reservoir.f90")
        fortran_compile_call.append(r"-m")
        fortran_compile_call.append(r"pymodule_levelpool")
        subprocess.run(
            fortran_compile_call,
            cwd=fortran_reservoir_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        from pymodule_levelpool import module_levelpool as rc
    except Exception as e:
        print(e)
        traceback.print_exc()
else:
    from pymodule_levelpool import module_levelpool as rc


def compute_network(
    terminal_segment=None,
    network=None,
    supernetwork_data=None,
    waterbody=None,
    nts=0,
    dt=0,
    verbose=False,
    debuglevel=0,
    write_csv_output=False,
    write_nc_output=False,
    assume_short_ts=False,
):
    global connections
    global flowveldepth
    if debuglevel <= -1:
        print(
            f"\nExecuting simulation on network {terminal_segment} beginning with streams of order {network['maximum_reach_seqorder']}"
        )

    ordered_reaches = {}
    for head_segment, reach in network["reaches"].items():
        if reach["seqorder"] not in ordered_reaches:
            ordered_reaches.update({reach["seqorder"]: []})
        ordered_reaches[reach["seqorder"]].append([head_segment, reach])

    # initialize write to files variable
    writeToCSV = write_csv_output
    writeToNETCDF = write_nc_output
    pathToOutputFile = os.path.join(root, "test", "output", "text")

    for ts in range(0, nts):
        for x in range(network["maximum_reach_seqorder"], -1, -1):
            for head_segment, reach in ordered_reaches[x]:
                # print(f'timestep: {ts}\n')
                # for loops should contain both waterbodies and mc_reach as it loops entire network
                qup_reach, quc_reach = compute_reach_upstream_flows(
                    head_segment=head_segment,
                    reach=reach,
                    supernetwork_data=supernetwork_data,
                    ts=ts,
                    dt=dt,
                    verbose=verbose,
                    debuglevel=debuglevel,
                    assume_short_ts=assume_short_ts,
                )
                if waterbody:
                    compute_level_pool_reach_up2down(
                        qup_reach=qup_reach,
                        quc_reach=quc_reach,
                        head_segment=head_segment,
                        reach=reach,
                        supernetwork_data=supernetwork_data,
                        waterbody=waterbody,
                        ts=ts,
                        dt=dt,
                        verbose=verbose,
                        debuglevel=debuglevel,
                        assume_short_ts=assume_short_ts,
                    )

                else:
                    compute_mc_reach_up2down(
                        qup_reach=qup_reach,
                        quc_reach=quc_reach,
                        head_segment=head_segment,
                        reach=reach,
                        supernetwork_data=supernetwork_data,
                        ts=ts,
                        dt=dt,
                        verbose=verbose,
                        debuglevel=debuglevel,
                        assume_short_ts=assume_short_ts,
                    )

    if writeToCSV:
        for x in range(network["maximum_reach_seqorder"], -1, -1):
            for head_segment, reach in ordered_reaches[x]:
                writeArraytoCSV(
                    connections=connections,
                    flowveldepth=flowveldepth,
                    reach=reach,
                    verbose=verbose,
                    debuglevel=debuglevel,
                    pathToOutputFile=pathToOutputFile,
                )

    if writeToNETCDF:
        writeArraytoNC(
            connections=connections,
            flowveldepth=flowveldepth,
            network=network,
            nts=nts,
            dt=dt,
            verbose=verbose,
            debuglevel=debuglevel,
            pathToOutputFile=pathToOutputFile,
        )


# TODO: generalize with a direction flag
def compute_reach_upstream_flows(
    head_segment=None,
    reach=None,
    supernetwork_data=None,
    ts=0,
    dt=60,
    verbose=False,
    debuglevel=0,
    assume_short_ts=False,
):
    global connections
    global flowveldepth

    if debuglevel <= -2:
        print(
            f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})"
        )

    # upstream flow per reach
    qup = 0.0
    quc = 0.0
    ########################
    if reach["upstream_reaches"] != {
        supernetwork_data["terminal_code"]
    }:  # Not Headwaters
        for us in connections[reach["reach_head"]]["upstreams"]:
            quc += flowveldepth[us]["flowval"][-1]
            if ts > 0:
                qup += flowveldepth[us]["flowval"][-2]

    if assume_short_ts:
        quc = qup

    return quc, qup


# TODO: generalize with a direction flag
def compute_mc_reach_up2down(
    qup_reach,
    quc_reach,
    head_segment=None,
    reach=None,
    supernetwork_data=None,
    ts=0,
    dt=60,
    verbose=False,
    debuglevel=0,
    assume_short_ts=False,
):
    global connections
    global flowveldepth

    if debuglevel <= -2:
        print(
            f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})"
        )

    qup = qup_reach
    quc = quc_reach

    current_segment = reach["reach_head"]
    # next_segment = connections[current_segment]["downstream"]

    while True:
        data = connections[current_segment]["data"]
        # current_flow = flowveldepth[current_segment]

        # for now treating as constant per reach
        bw = data[supernetwork_data["bottomwidth_col"]]
        tw = data[supernetwork_data["topwidth_col"]]
        twcc = data[supernetwork_data["topwidthcc_col"]]
        dx = data[supernetwork_data["length_col"]]
        n_manning = data[supernetwork_data["manningn_col"]]
        n_manning_cc = data[supernetwork_data["manningncc_col"]]
        cs = data[supernetwork_data["ChSlp_col"]]
        s0 = data[supernetwork_data["slope_col"]]

        qlat = flowveldepth[current_segment]["qlatval"][ts]

        if ts > 0:
            qdp = flowveldepth[current_segment]["flowval"][-1]
            velp = flowveldepth[current_segment]["velval"][-1]
            depthp = flowveldepth[current_segment]["depthval"][-1]
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
            qup = flowveldepth[current_segment]["flowval"][-1]  # input for next segment
        else:
            qup = 0

        quc = qdc  # input for next segment
        if assume_short_ts:
            quc = qup
        # need to append with waterbody calculations instead
        flowveldepth[current_segment]["flowval"].append(qdc)
        flowveldepth[current_segment]["depthval"].append(depthc)
        flowveldepth[current_segment]["velval"].append(velc)
        flowveldepth[current_segment]["time"].append(ts * dt)

        next_segment = connections[current_segment]["downstream"]
        if current_segment == reach["reach_tail"]:
            if debuglevel <= -2:
                print(f"{current_segment} (tail)")
            break
        if debuglevel <= -3:
            print(f"{current_segment} --> {next_segment}\n")
        # current_segment = next_segment
        # next_segment = connections[current_segment]["downstream"]
        current_segment = next_segment
    # end loop initialized the MC vars
    # end while loop


def compute_level_pool_reach_up2down(
    qup_reach,
    quc_reach,
    head_segment=None,
    reach=None,
    supernetwork_data=None,
    waterbody=None,
    ts=0,
    dt=60,
    verbose=False,
    debuglevel=0,
    assume_short_ts=False,
):
    global connections
    global flowveldepth
    global waterbodies_df

    if debuglevel <= -2:
        print(
            f"\nreach: {head_segment} (order: {reach['seqorder']} n_segs: {len(reach['segments'])})"
        )

    qup = qup_reach
    quc = quc_reach

    current_segment = reach["reach_tail"]
    if ts > 0:
        depthp = flowveldepth[current_segment]["depthval"][-1]
    else:
        depthp = 0

    qlat = flowveldepth[current_segment]["qlatval"][ts]
    if debuglevel <= -2:
        print(f"executing reservoir computation on waterbody: {waterbody}")

    ln = waterbody
    qi0 = qup
    qi1 = quc
    ql = qlat
    dt = dt  # current timestep
    ar = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_waterbody_area"]]
    we = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_weir_elevation"]]
    maxh = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_waterbody_max_elevation"]]
    wc = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_outfall_weir_coefficient"]]
    wl = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_outfall_weir_length"]]
    dl = 10 * wl  # waterbodies_df.loc[waterbody][supernetwork_data["level_pool_overall_dam_length"]]
    oe = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_orifice_elevation"]]
    oc = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_orifice_coefficient"]]
    oa = waterbodies_df.loc[waterbody][supernetwork_data["level_pool_orifice_area"]]

    import pdb; pdb.set_trace()
    qdc, depthc = rc.levelpool_physics(
        ln, qi0, qi1, ql, dt, depthp, ar, we, maxh, wc, wl, dl, oe, oc, oa,
    )

    flowveldepth[current_segment]["flowval"].append(qdc)
    flowveldepth[current_segment]["depthval"].append(depthc)
    flowveldepth[current_segment]["velval"].append(0)
    flowveldepth[current_segment]["time"].append(ts * dt)


# ### Psuedocode
# Write Array  to CSV file
# arguments reach , pathToOutputFile
# using global connections and flowveldepth.
def writeArraytoCSV(
    flowveldepth=None,
    connections=None,
    reach=None,
    verbose=False,
    debuglevel=0,
    pathToOutputFile="../../test/output/text",
):

    # define CSV file Header
    header = ["time", "qlat", "q", "v", "d"]

    # Loop over reach segments
    current_segment = reach["reach_head"]
    next_segment = connections[current_segment]["downstream"]

    while True:
        filename = f"{pathToOutputFile}/{current_segment}.csv"  #
        if verbose:
            print(f"writing segment output to --> {filename}")
        with open(filename, "w+") as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=",", quoting=csv.QUOTE_ALL)
            csvwriter.writerow(header)
            csvwriter.writerows(
                zip(
                    flowveldepth[current_segment]["time"],
                    flowveldepth[current_segment]["qlatval"],
                    flowveldepth[current_segment]["flowval"],
                    flowveldepth[current_segment]["velval"],
                    flowveldepth[current_segment]["depthval"],
                )
            )

        if current_segment == reach["reach_tail"]:
            if debuglevel <= -2:
                print(f"{current_segment} (tail)")
            break
        if debuglevel <= -3:
            print(f"{current_segment} --> {next_segment}\n")
        current_segment = next_segment
        next_segment = connections[current_segment]["downstream"]


# ### Psuedocode
# Write Array  to Arrays for Netcdf file and then call writeNC function to write output data to NC file
# arguments network, number of timesteps (nts),  timestep in seconds  (dt),  pathToOutputFile
# using global connections and flowveldepth.
def writeArraytoNC(
    flowveldepth=None,
    connections=None,
    network=None,
    nts=0,
    dt=60,
    verbose=False,
    debuglevel=0,
    pathToOutputFile="../../test/output/text/",
):
    # create  array variables to copy from python "flowveldepth" which is global
    flowveldepth_data = {
        "segment": [],
        "time": [],
        "qlatval": [],
        "flowval": [],
        "depthval": [],
        "velval": [],
    }

    ordered_reaches = {}
    for head_segment, reach in network["reaches"].items():
        if reach["seqorder"] not in ordered_reaches:
            ordered_reaches.update(
                {reach["seqorder"]: []}
            )  # TODO: Should this be a set/dictionary?
        ordered_reaches[reach["seqorder"]].append([head_segment, reach])

    # get data into array - preparation step
    TIME_WRITTEN = False
    write_segment = None
    for x in range(network["maximum_reach_seqorder"], -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            current_segment = reach["reach_head"]
            while True:
                # appending data from each segments to a single list  "flowveldepth_data"
                # preserving ordering same as segment in a reach
                flowveldepth_data["qlatval"].append(
                    flowveldepth[current_segment]["qlatval"]
                )
                flowveldepth_data["flowval"].append(
                    flowveldepth[current_segment]["flowval"]
                )
                flowveldepth_data["depthval"].append(
                    flowveldepth[current_segment]["depthval"]
                )
                flowveldepth_data["velval"].append(
                    flowveldepth[current_segment]["velval"]
                )
                # write segment flowveldepth_data['segment']
                flowveldepth_data["segment"].append(current_segment)
                if not TIME_WRITTEN:
                    # write time only once - for first segment
                    flowveldepth_data["time"].append(
                        flowveldepth[current_segment]["time"]
                    )
                    TIME_WRITTEN = True

                if current_segment == reach["reach_tail"]:
                    write_segment = current_segment
                    if debuglevel <= -2:
                        print(f"{current_segment} (tail)")
                    break
                next_segment = connections[current_segment]["downstream"]
                if debuglevel <= -3:
                    print(f"{current_segment} --> {next_segment}\n")
                current_segment = next_segment

    # check number of timesteps should match the time the data is written
    if len(flowveldepth_data["time"][0]) != nts:
        print(
            f"Number of timesteps  {nts} does not match data timesteps {len(flowveldepth_data['time'][0])}\n"
        )
        return

    writeNC(
        flowveldepth_data=flowveldepth_data,
        segment_count=network["total_segment_count"],
        terminal_segment=write_segment,
        nts=nts,
        dt=dt,
        pathToOutputFile=pathToOutputFile,
        verbose=verbose,
        debuglevel=debuglevel,
    )


# ### Psuedocode
# Write  netcdf  file
# arguments flowveldepth , segment_count, terminal_segment (for filename and as identifier of reach)
#           number of timesteps (nts),  timestep in seconds  (dt),  verbose , debuglevel , pathToOutputFile
def writeNC(
    flowveldepth_data=None,
    segment_count=0,
    terminal_segment=None,
    nts=0,
    dt=60,
    pathToOutputFile="../../test/output/text",
    verbose=False,
    debuglevel=0,
):
    # start writing data to nc file
    filename = f"{pathToOutputFile}/{terminal_segment}.nc"  # ncfile'
    if verbose:
        print(f"writing netcdf output to --> {filename}")
    ncfile = netCDF4.Dataset(filename, mode="w", format="NETCDF4")
    # segcount = total segments for the current reach
    segcount = ncfile.createDimension("stations", segment_count)  # segment
    # timecount = number of timesteps
    timecount = ncfile.createDimension(
        "time", nts
    )  # unlimited axis (can be appended to).
    # analysis time =  current time , hence count =1
    analysistimecount = ncfile.createDimension(
        "analysistime", 1
    )  # unlimited axis (can be appended to).
    ncfile.title = f"Result of MC for Reach with terminal segment {terminal_segment}"
    ncfile.subtitle = "MC Python Module"
    ncfile.anything = f"streamflow , depth, velocity and lateral flow for Reach with terminal segment {terminal_segment}"
    # write streamflow
    flow = ncfile.createVariable(
        "flow", np.float64, ("time", "stations")
    )  # note: unlimited dimension is leftmost
    flow[:, :] = np.transpose(np.array(flowveldepth_data["flowval"], dtype=float))
    flow.units = "cu ft/s"
    flow.standard_name = "streamflow"  # this is a CF standard name
    # write depth
    depth = ncfile.createVariable(
        "depth", np.float64, ("time", "stations")
    )  # note: unlimited dimension is leftmost
    depth[:, :] = np.transpose(np.array(flowveldepth_data["depthval"], dtype=float))
    depth.units = "ft"  #
    depth.standard_name = "depth"  # this is a CF standard name
    # write velocity
    velocity = ncfile.createVariable(
        "velocity", np.float64, ("time", "stations")
    )  # note: unlimited dimension is leftmost
    velocity.units = "ft/s"  #
    velocity[:, :] = np.transpose(np.array(flowveldepth_data["velval"], dtype=float))
    velocity.standard_name = "velocity"  # this is a CF standard name
    # write  lateral flow (input from NWM)
    lateralflow = ncfile.createVariable(
        "lateralflow", np.float64, ("time", "stations")
    )  # note: unlimited dimension is leftmost
    lateralflow[:, :] = np.transpose(
        np.array(flowveldepth_data["qlatval"], dtype=float)
    )
    lateralflow.units = "cu ft/s"  #
    lateralflow.standard_name = "lateralflow"  # this is a CF standard name
    # write time in seconds since  TODO get time from lateral flow from NWM
    time = ncfile.createVariable("time", np.float64, "time")
    time.units = "seconds since 2011-08-27 00:00:00"  ## TODO get time fron NWM as argument to this function
    time.long_name = "time"
    time[:] = flowveldepth_data["time"]
    # write segment ids
    segment = ncfile.createVariable("station_id", np.int, ("stations"))
    segment.long_name = "feature id"
    segment[:] = flowveldepth_data["segment"]
    # write analysis time = current time = time MC model is run
    analysistime = ncfile.createVariable("analysis_time", np.double, "analysistime")
    analysistime.long_name = "forecast_reference_time"
    analysistime.units = " minutes since " + datetime.now().strftime(
        "%Y-%m-%d %H:%M:%S"
    )
    analysistime[:] = 0
    #
    ncfile.close()
    if debuglevel <= -1:
        print(f"{filename} closed!")


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
    global flowveldepth
    global waterbodies_df

    supernetwork = args.supernetwork
    break_network_at_waterbodies = args.break_network_at_waterbodies

    dt = float(args.dt)
    nts = int(args.nts)
    qlat_const = float(args.qlat_const)

    debuglevel = -1 * int(args.debuglevel)
    verbose = args.verbose
    showtiming = args.showtiming
    write_csv_output = args.write_csv_output
    write_nc_output = args.write_nc_output
    assume_short_ts = args.assume_short_ts
    parallel_compute = args.parallel_compute

    run_pocono_test = args.run_pocono_test

    if run_pocono_test:
        if verbose:
            print("running test case for Pocono_TEST2 domain")
        # Overwrite the following test defaults
        supernetwork = "Pocono_TEST2"
        break_network_at_waterbodies = False
        dt = 300
        nts = 144
        write_csv_output = True
        write_nc_output = True

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

    # initialize flowveldepth dict
    flowveldepth = {
        connection: {
            "qlatval": [],
            "time": [],
            "flowval": [],
            "velval": [],
            "depthval": [],
        }
        for connection in connections
    }

    # Lateral flow
    if (
        run_pocono_test
    ):  # test 1. Take lateral flow from wrf-hydro output from Pocono Basin
        ql_input_folder = os.path.join(
            root, r"test/input/geo/PoconoSampleData2/Pocono_ql_testsamp1_nwm_mc.csv"
        )
        ql = pd.read_csv(ql_input_folder, index_col=0)

    else:
        ql = pd.DataFrame(
            qlat_const, index=connections.keys(), columns=range(nts), dtype="float32"
        )

    for index, row in ql.iterrows():
        flowveldepth[index]["qlatval"] = row.tolist()
    ######################
    waterbodies_values = supernetwork_values[12]
    waterbodies_segments = supernetwork_values[13]

    waterbodies_df = nnu.read_waterbody_df(
        supernetwork_data["level_pool_waterbody_parameter_file_path"],
        waterbodies_values,
    )

    if not parallel_compute:
        if verbose:
            print("executing computation on ordered reaches ...")

        #########

        for terminal_segment, network in networks.items():

            # is_reservoir = terminal_segment in waterbodies_segments
            waterbody = waterbodies_segments.get(terminal_segment)

            compute_network(
                terminal_segment=terminal_segment,
                network=network,
                supernetwork_data=supernetwork_data,
                waterbody=waterbody,
                nts=nts,
                dt=dt,
                verbose=verbose,
                debuglevel=debuglevel,
                write_csv_output=write_csv_output,
                write_nc_output=write_nc_output,
                assume_short_ts=assume_short_ts,
            )
            if showtiming:
                print("... in %s seconds." % (time.time() - start_time))
    else:  # serial execution
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
                nts,
                dt,
                verbose,
                debuglevel,
                write_csv_output,
                write_nc_output,
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
