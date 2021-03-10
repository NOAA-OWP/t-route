## testing: python nhd_network_augment.py CapeFear_FULL_RES -return_original -p -s
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from functools import partial
from itertools import chain
import os
import sys
import pathlib
import argparse
import json
from tqdm import tqdm
import time
from memory_profiler import profile

root = pathlib.Path("../../../").resolve()
sys.path.append(os.path.join(root, "src", "python_framework_v01"))

import nhd_network_utilities_v02 as nnu
import nhd_network
import nhd_io
import network_dl


def _handle_args():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-t",
        "--threshold_length",
        help="threshold segment length (meters)",
        dest="threshold",
        default=500,
        type=int,
    )
    parser.add_argument(
        "--network",
        help="Choose from among the pre-programmed supernetworks (Florence_FULL_RES, CONUS_FULL_RES_v20, CapeFear_FULL_RES)",
        dest="supernetwork",
        choices=[
            "CONUS_FULL_RES_v20",
            "CapeFear_FULL_RES",
            "Florence_FULL_RES",
        ],
        default="Florence_FULL_RES",
    )
    parser.add_argument(
        "-p",
        "--prune",
        help="prune short headwater reaches",
        dest="prune",
        action="store_true",
    )
    parser.add_argument(
        "-s",
        "--snap",
        help="snap junctions adjacent to short reaches",
        dest="snap",
        action="store_true",
    )
    parser.add_argument(
        "-return_original",
        "--return_original",
        help="return an unmodified RouteLink.nc file for the specified domain",
        dest="return_original",
        action="store_true",
    )
    parser.add_argument(
        "-write_output",
        "--write_output",
        help="write resulting augmented routelink data",
        dest="write_output",
        action="store_true",
    )

    return parser.parse_args()


def get_network_data(network_name):

    # Create directory path variable for test/input/geo, where NHD data and masks are stored
    test_folder = os.path.join(root, r"test")
    geo_input_folder = os.path.join(test_folder, r"input", r"geo")

    # Load network meta data 
    supernetwork = network_name
    network_data = nnu.set_supernetwork_parameters(
        supernetwork=supernetwork, geo_input_folder=geo_input_folder
    )

    # if the NHDPlus RouteLink file does not exist, download it.
    if not os.path.exists(network_data["geo_file_path"]):
        filename = os.path.basename(network_data["geo_file_path"])
        network_dl.download(network_data["geo_file_path"], network_data["data_link"])

    # read-in NHD network data, as-is, unmodified
    RouteLink = nhd_io.read(network_data["geo_file_path"])

    # select only the necessary parameter columns
    cols = list(network_data["columns"].values())
    # GET THE STRAHLER ORDER DATA TOO!
    cols.append("order")

    data = nhd_io.read(network_data["geo_file_path"])    
    data = data[cols]
    data = data.set_index(network_data["columns"]["key"])

    # mask to supernetwork  extent
    if "mask_file_path" in network_data:
        data_mask = nhd_io.read_mask(
            network_data["mask_file_path"],
            layer_string=network_data["mask_layer_string"],
        )
        data = data.filter(data_mask.iloc[:, network_data["mask_key"]], axis=0)

    # sort index
    data = data.sort_index()

    # replace downstreams
    data = nhd_io.replace_downstreams(data, network_data["columns"]["downstream"], 0)

    return data, RouteLink, network_data

def network_connections(data, network_data):

    """
    Extract upstream and downstream connections between segments in network
    Args:
        data (DataFrame): Network parameter dataset, prepared
        network_data (dict): network metadata
    Returns:
        conn (dict): downstream connections
        rconn (dict): upstream connections 
    """

    # extract downstream connections
    conn = nhd_network.extract_connections(data, network_data["columns"]["downstream"])

    # extract upstream connections
    rconn = nhd_network.reverse_network(conn)

    return conn, rconn


def build_reaches(rconn):
    
    """
    Construct reaches from network connections
    Args:
        rconn (dict): upstream connections
    Returns:
        sureachable (dict): lists of all segments (values) connected to each network tailwater (keys)
        subreaches (dict): list of reaches (linearly connected segments) connected to each network tailwater (keys)
    """

    # isolate independent subnetworks
    subnets = nhd_network.reachable_network(rconn)

    # identify the segments in each subnetwork
    subreachable = nhd_network.reachable(rconn)

    # break each subnetwork into reaches
    subreaches = {}
    for tw, net in subnets.items():
        path_func = partial(nhd_network.split_at_junction, net)
        subreaches[tw] = nhd_network.dfs_decomposition(net, path_func)
    
    return subreachable, subreaches

def find_headwater_reaches(subreaches, hws):
    
    """
    Construct a list of headwater reaches
    Args:
        subreaches (dict): reaches (linearly connected segments) connected to each network tailwater (keys)
        hws (list): headwater segments
    Returns:
        hw_reaches (list): lists of headwater reaches 
    """
    
    hw_reaches = []
    for sublists in list(subreaches.values()):
        for rch in sublists:
            for val in rch:
                if val in hws:
                    hw_reaches.append(rch)
                    pass
    
    return hw_reaches

def find_short_headwater_reaches(hw_reaches, idx_array, length_array, threshold):

    """
    Construct a list of short headwater reaches
    Args:
        hw_reaches (list): lists of headwater reaches 
        idx_array (numpy array): ordered network segment indices (link IDs)
        length_array (numpy array): network segment lengths, ordered by idx_array
        threshold (int): threshold length that defines short
    Returns:
        short_hw_reaches (list): lists of SHORT headwater reaches 
    """
    
    short_hw_reaches = []
    for rch in hw_reaches:

        # find row locations of reach segments
        a = np.searchsorted(idx_array, rch)

        # check if reach is shorter than threshold
        if np.sum(length_array[a]) < threshold:

            # note reach as being a short headwater
            short_hw_reaches.append(rch)
            
    return short_hw_reaches

def headwater_junctions(short_hw_reaches, conn):
    
    """
    Build a dictionary of headwater junctions (keys) and contributing reaches (values)
    Args:
        short_hw_reaches (list): lists of short headwater reaches 
        conn (dict): downstream connections
    Returns:
        hw_junctions (dict): headwater reaches (values) and the junction nodes they drain to (keys)
    """
        
    hw_junctions = {}
    for rch in short_hw_reaches:
        jun = conn[rch[-1]]
        hw = rch

        if len(jun) == 0:
            # no downstream connection
            hw_junctions[rch[-1]] = [hw]

        elif jun[0] in hw_junctions.keys():
            # junction already eists as dict key
            hw_junctions[jun[0]].append(hw)
        else:
            # new junction
            hw_junctions[jun[0]] = [hw]
            
    return hw_junctions

def prune_headwaters(data, threshold, network_data):
    
    """
    Remove short headwater reaches from the network. Short is defined as having a length less that the threshold
    Args:
        data (pandas DataFrame): network parameters
        threshold (int): threshold length that defines short
        network_data (dict): network metadata
    Returns:
        data_pruned (pandas DataFrame) : parameters for network with short headwater reaches removed
    """

    # initialize list to store pruned reaches
    hw_prune_list = []
    
    # numpy arrays for quick indexing
    idx_array = data.index.to_numpy()
    length_array = data.Length.to_numpy()
    
    iter_count = 1
    # continue so long as short headwaters reaches exist in the network
    while 1 == 1:
        
        t1 = time.time()
        # STEP 1: Find headwater reaches:
        # --------------------------------------#
        
        # flatten headwater prune list
        flat_list = [item for sublist in hw_prune_list for item in sublist]
        flatter_list = [item for sublist in flat_list for item in sublist]
        
        # build connections and reverse connections
        connections, rconn = network_connections(
            data.drop(flatter_list), # network data, with pruned headwater reaches removed
            network_data
        )

        # identify headwater segments
        hws = nhd_network.headwaters(connections)

        # build reaches
        subreachable, subreaches = build_reaches(rconn)

        # find headwater reaches
        hw_reaches = find_headwater_reaches(subreaches, hws)

        # STEP 2: identify short headwater reaches
        # --------------------------------------#
    
        # find headwater reaches shorter than threshold
        short_hw_reaches = find_short_headwater_reaches(hw_reaches, idx_array, length_array, threshold)
        
        # CHECK: Are there any short headwter reaches?
        # --------------------------------------#
        if len(short_hw_reaches) == 0:
            print("no more short headwaters to prune")

            # if no more reaches, exit while loop
            break

        # STEP 3: trim short headwater reaches
        # --------------------------------------
        
        # build a dictionary of headwater junction segments (keys) and contributing reaches (values) 
        hw_junctions = headwater_junctions(short_hw_reaches, connections)
            
        for jun, hw in hw_junctions.items():
            
            # do two short headwaters converge at this junction?
            if len(hw) > 1:
                # if so, prune the shorter of the two
            
                # length of reach 1
                a = np.searchsorted(idx_array, hw[0])
                r1_len = np.sum(length_array[a])
                
                # length of reach 2
                a = np.searchsorted(idx_array, hw[1])
                r2_len = np.sum(length_array[a])
                
                if r1_len < r2_len:
                    # prune reach 1
                    hw_prune_list.append([hw[0]])
                else:
                    # prune reach 2
                    hw_prune_list.append([hw[1]])
                    
            else:
                # if not, just prune the headwater reach
                hw_prune_list.append(hw)
        
        t2 = time.time()
        print("completed", iter_count, "iterations of headwater pruning in", np.round(t2 - t1, 3),)
        iter_count += 1

    # remove pruned headwater segments from the network
    flat_list = [item for sublist in hw_prune_list for item in sublist]
    flatter_list = [item for sublist in flat_list for item in sublist]
    data_pruned = data.drop(flatter_list).copy()
    
    return data_pruned

def find_short_nonheadwater_reaches(subreaches,rconn,idx_array,length_array,threshold):
    
    """
    Construct a list of short non-headwater reaches
    Args:
        subreaches (list): network reaches 
        rconn(dict): upstream connections
        idx_array (numpy array): ordered network segment indices (link IDs)
        length_array (numpy array): network segment lengths, ordered by idx_array
        threshold (int): threshold length that defines short
    Returns:
        short_reaches (tuple): short non-headwater reaches
    """
    
    short_reaches = ()
    for e, sublists in enumerate(list(subreaches.values())):
        for rch in sublists:
            if rconn[rch[0]]:
                # the reach is not a headwater
                # check reach length
                a = np.searchsorted(idx_array, rch)
                if np.sum(length_array[a]) < threshold:
                    short_reaches += tuple(rch)
                
    return(short_reaches)

def snap_junctions(data, threshold, network_data):

    """
    This function snaps junctions on oposite ends of a short reach
    by forcing the lowest order upstream tributary to drain to the reach tail. 
    Short reaches are defined by a user-appointed threshold length.

    For example, consider a short reach (*) stranded between two junctions:

          \  /
           \/ | 4th
       2nd  \ |  
             \| / 
             *|/
              |
              | 

    The algoritm would select the lowest order tributary segment to the reach head, 
    which in this case is a second-order drainage, and change the downstream connection
    to the reach tail. This produces the following "snapped" network config:

          \  / 
           \/ |  
            \ | / 
             \|/
              |
              | 

    Inputs
    --------------
    - data (DataFrame): NHD RouteLink data, must contain Strahler order attribute
        - this input can be either the native NHD data OR the data with pruned headwaters

    Outputs
    --------------
    - snapped_data (DataFrame): NHD RouteLink data containing no short reaches

    """
    # create a copy of the native data for snapping
    data_snapped = data.copy()

    # snap junctions
    iter_num = 1
    while 1 == 1:

        t1 = time.time()
        
        # evaluate connections
        connections, rconn = network_connections(data_snapped, network_data)

        # build reaches
        subreachable, subreaches = build_reaches(rconn)

        # build list of headwater segments
        hws = nhd_network.headwaters(connections)
        
        # find short reaches stranded between junctions
        idx_array = data_snapped.index.to_numpy()
        length_array = data_snapped.Length.to_numpy()
        short_reaches = find_short_nonheadwater_reaches(subreaches,rconn,idx_array,length_array,threshold)

        print(
            "Prior to itteration",
            iter_num,
            ",",
            len(short_reaches),
            "short, non-headwater, reaches exist in the network",
        )
        
        # check that short reaches exist, if none - terminate process
        if len(short_reaches) == 0:
            break

        # for each short reach, snap lower order upstream trib to downstream drainage destination
        order_array = data_snapped.order.to_numpy()
        to_array = data_snapped.to.to_numpy()
        rch_to_move = []
        snap_destination = []
        for i, rch in enumerate(short_reaches):

            # identify reach tail (downstream-most) and head (upstream-most) segments
            if type(rch) is tuple:
                tail = rch[-1]
                head = rch[0]
            else:
                tail = rch
                head = rch

            # identify segments that drain to reach head
            us_conn = rconn[head]

            # select the upstream segment to snap
            a = np.searchsorted(idx_array, us_conn)
            o = order_array[a]
            l = length_array[a]            
            if min(o) != max(o):
                rch_to_move.extend([idx_array[a[np.argmin(o)]]])
            else:
                rch_to_move.extend([idx_array[a[np.argmin(l)]]])

            # if reach tail doesn't have a downstrem connection (it is a network tailwater)
            # then the snap desination is the ocean
            if len(connections[tail]) == 0:
                a = np.searchsorted(idx_array, tail)
                snap_destination.extend([to_array[a]])
            else:
                snap_destination.extend(connections[tail])
        
        data_snapped.loc[rch_to_move, "to"] = snap_destination
        t2 = time.time()
        print("completed", iter_num, "iterations of junction snapping in", np.round(t2 - t1, 3),)
        iter_num += 1

    return data_snapped

def qlat_destination_compute(
    data_native, data_merged, merged_segments, pruned_segments, network_data
):
    
    idx_array = data_merged.index.to_numpy()
    
    # build a list of all segments that need crosswalking
    if bool(list(pruned_segments)):
        segments = list(merged_segments) + list(pruned_segments)

    else:
        segments = list(merged_segments)

    # compute connections using native network data
    conn = nhd_network.extract_connections(
        data_native, network_data["columns"]["downstream"]
    )
    rconn = nhd_network.reverse_network(conn)

    # initialize a dictionary to store qlat destination nodes for pruned/merged segments
    qlat_destinations = {}

    for idx in segments:

        # find the segment to recieve qlats from the pruned or merged segment
        if conn[idx]:
            ds_idx = conn[idx]
            while bool(ds_idx[0] in idx_array) == False:
                ds_idx = conn[ds_idx[0]]

        elif rconn[idx]:
            ds_idx = rconn[idx]
            while bool(ds_idx[0] in idx_array) == False:
                us_idx = conn[ds_idx[0]]

        else:
            ds_idx = []

        # update the qlat destination dict
        qlat_destinations[str(idx)] = str(ds_idx)

    return qlat_destinations

def segment_merge(data_native, data, network_data, thresh, pruned_segments):

    # create a copy of the pruned network dataset, which will be updated with merged data
    data_merged = data.copy()

    # build connections and reverse connections
    conn, rconn = network_connections(data, network_data)

    # organize network into reaches
    subreachable, subreaches = build_reaches(rconn)
    
    # numpy arrays of indices and segment lengths
    idx_array = data.index.to_numpy()
    len_array = data.Length.to_numpy()

    # loop through each reach in the network
    print("evaluating segment merge orientation")
    t1 = time.time()
    merges = {}
    for twi, (tw, rchs) in enumerate(subreaches.items(), 1):

        for rch in rchs:
            
            # numpy arrays of index and length for this reach
            a = np.searchsorted(idx_array, rch)
            idx_rch = idx_array[a]
            len_rch = len_array[a]
            
            rch_len = sum(len_rch)
            L = data.loc[rch]
            
            if min(len_rch) < thresh:

                ##################################################
                # orphaned short single segment reaches
                ##################################################
                # if reach length is shorter than threshold and composed of a single segment
                if rch_len < thresh and len(idx_rch) == 1:
                    continue

                ##################################################
                # multi segment reaches - combine into a single segment reach
                ##################################################
                # if reach length is shorter than threshold and composed of more than one segment
                elif rch_len < thresh and len(idx_rch) > 1:

                    eater = idx_rch[-1]
                    eaten = idx_rch[:-1]
                    merges[eater] = eaten

                ##################################################
                # multi segment reaches longer than threshold with some segments shorter than threshold
                ##################################################
                # if reach length is longer than threshold and smallest segment length is less than threshold
                else:
                    
                    idx_merge = idx_rch
                    len_merge = len_rch
                    
                    while min(len_merge) < thresh:

                        # find the shortest segment in the reach
                        idx_short = idx_merge[np.argmin(len_merge)]    
                        iloc_short = np.argmin(len_merge)
                            
                        # ----- Determine segment merge order ------
                        # if shortest segment is reach tail, short segment eats upstream segment
                        if idx_short == idx_merge[-1]:
                            eater = idx_short
                            eaten = idx_merge[-2]

                        # if shortest segment is not reach tail, downstream segment eats short segment
                        else:
                            eater = idx_merge[iloc_short + 1]
                            eaten = idx_short
                            
                        # adjust length of eater segment
                        ieater = np.where(idx_merge == eater)
                        ieaten = np.where(idx_merge == eaten)
                        
                        len_merge[ieater] = len_merge[ieater] + len_merge[ieaten]
                        len_merge = np.delete(len_merge, ieaten)
                        idx_merge = np.delete(idx_merge, ieaten)

                        # ---- Append the merges dictionary - a running record of what is merged into what -----
                        # if eater and eaten are already listed as eaters
                        if eater in merges.keys() and eaten in merges.keys():
                            merges[eater] = merges[eater] + merges[eaten] + [eaten]
                            del merges[eaten] 
                        
                        # if eater is already listed as an eater
                        elif eater in merges.keys():
                            merges[eater] = merges[eater] + [eaten]

                        # if eaten is already listed as eater
                        elif eaten in merges.keys():
                            merges[eater] = merges[eaten] + [eaten]
                            del merges[eaten]
                            
                        else:
                            # update merge dictionary
                            merges[eater] = [eaten]

    t2 = time.time()
    print("Done evaluating merge orientation in", np.round(t2 - t1, 3),)
    
    # build lists of all eaten segments and their corresponding eaters
    # build a cross walk dictionary to help re parameterize downstream connections ('to')
    print("Creating 'to' cross walk dict and all_eaters/eaten lists")
    t1 = time.time()
    all_eaten = []
    all_eaters = []
    to_cw = {}
    for eater, eaten in merges.items():
        for e in eaten:
            to_cw[e] = eater
            all_eaten.append(e)
            all_eaters.append(eater)
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))

    # build a dataframe of eaten segments, with a new column for their eaters
    print("Building data_eaten dataframe")
    t1 = time.time()
    data_eaten = data.loc[all_eaten]
    data_eaten["eater"] = all_eaters
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # build a dataframe of eater segments, with a new eater column for the sake of grouping
    print("Building data_eaters dataframe")
    t1 = time.time()
    data_eaters = data.loc[np.unique(all_eaters)]
    data_eaters["eater"] = data_eaters.index
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # append eater and eaten dataframes together
    print("merging eater and eaten dataframes")
    t1 = time.time()
    data_to_merge = data_eaten.append(data_eaters)
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # group by eater segmment IDs and modify routing parameters
    print("Grouping by eater and changing parameter values")
    t1 = time.time()
    idx_array = data_to_merge.index.to_numpy()
    length_array = data_to_merge.Length.to_numpy()
    sort = np.argsort(idx_array)
    data_parameter_merge = data_to_merge.groupby("eater").agg({'to': 'min',
                                               'Length':'sum', 
                                               'n': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'nCC': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'So': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'BtmWdth': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'TopWdth': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'TopWdthCC': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'NHDWaterbodyComID': 'min',
                                               'MusK': 'min',
                                               'MusX': 'min',
                                               'ChSlp': lambda x: np.average(x, weights = length_array[sort[np.searchsorted(idx_array, x.index.to_numpy(), sorter = sort)]]),
                                               'order': 'min'})
    
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # update 'to' variable in data_merged
    print("updating 'to' variable")
    t1 = time.time()
    data_parameter_merge['to'] = data_to_merge.loc[data_parameter_merge.index]['to']
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # replace
    print("dropping and replacing merged segments")
    t1 = time.time()
    data_merged = data_merged.drop(all_eaten)  
    data_merged.loc[data_parameter_merge.index] = data_parameter_merge.loc[data_parameter_merge.index]
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # adjust 'to' variable
    print("re-wiring connections")
    t1 = time.time()
    to_array = data_merged.to.to_numpy()
    idx_array = data_merged.index.to_numpy()
    idx = []
    to_replace = []
    for i, (eaten, eater) in enumerate(to_cw.items()):
        if eaten in to_array:

            a = np.where(to_array == eaten)
            idx.extend(list(idx_array[a]))
            to_replace.extend([eater]*len(list(idx_array[a])))

    data_merged.loc[idx,'to'] = to_replace
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    # create a qlateral destinations dictionary
    print("building qlateral crosswalk")
    t1 = time.time()
    qlat_destinations = qlat_destination_compute(
        data_native, data_merged, all_eaten, pruned_segments, network_data
    ) 
    t2 = time.time()
    print("Done. That took:",np.round(t2-t1,3))
    
    return data_merged, qlat_destinations

def main():

    # unpack command line arguments
    args = _handle_args()
    supernetwork = args.supernetwork
    threshold = args.threshold
    prune = args.prune
    snap = args.snap
    return_original = args.return_original
    write_output = args.write_output

    # get network data
    print("Extracting and organizing supernetwork data")
    data, RouteLink, network_data = get_network_data(supernetwork)
    RouteLink = RouteLink.set_index(network_data["columns"]["key"])

    pruned_segs = set()
    
    #-------------------------------
    # PRUNE - SNAP - MERGE
    #-------------------------------
    t1 = time.time()
    if prune and snap:
        dirname = (
            "RouteLink_"
            + supernetwork
            + "_" + str(threshold)
            + "m_prune_snap_merge"
        )
        filename_cw = (
            "CrossWalk_"
            + supernetwork
            + "_"
            + str(threshold)
            + "m_prune_snap_merge.json"
        )

        print("------------- PRUNING HEADWATERS -------------")
        data_pruned = prune_headwaters(data, threshold, network_data)

        # identify pruned segments
        pruned_segs.update(set(np.setdiff1d(data.index, data_pruned.index)))
        
        print("------------- SNAPPING JUNCTIONS -------------")
        data_snapped = snap_junctions(data_pruned, threshold, network_data)

        print("------------- MERGING SEGMENTS -------------")
        data_merged, qlat_destinations = segment_merge(
            data, data_snapped, network_data, threshold, pruned_segs
        )
        
    print("------------- AUGMENTATION COMPLETE -------------")
    t2 = time.time()
    print("The total augmentation processing time:", np.round(t2-t1,3))
    
    print("------------- QA/QC TESTING -------------")
    # Check augmented network connections
    def connection_check(N1, N2, network_data):
        
        def count_tailwaters(N, network_data):
            connections, rconn = network_connections(N, network_data)
            subreachable, subreaches = build_reaches(rconn)
            return(len(subreachable.keys()))
        
        tws_N1 = count_tailwaters(N1, network_data)
        tws_N2 = count_tailwaters(N2, network_data)
        return(tws_N1 == tws_N2)
    
    check = connection_check(data_snapped, data_merged, network_data)
    if check:
        print("Passed connection test")
    else:
        print("Error! - The merging process has created additional tailwaters (i.e. it broke the network)")
        raise ValueError
        
    def crosswalk_check(data, data_merged):
        
        # make sure all segments removed from the original network have a qlateral destination in the augmented network
        
        # which segments no longer exist in the augmented network?
        idx_original = data.index.to_numpy()
        idx_augmented = data_merged.index.to_numpy()
        
        mask = np.in1d(idx_original, idx_augmented, invert = True)
        idx_dropped = idx_original[mask]
        
        # check that all dropped ind
        if idx_dropped 
        
        print(len(idx_dropped))
        print(len(data))
        print(len(qlat_destinations.keys()))
        raise ValueError
        
    a = crosswalk_check(data,data_merged)

    #-------------------------------
    #      SNAP - MERGE
    #-------------------------------
    if snap and not prune:
        dirname = (
            "RouteLink_"
            + supernetwork
            + "_"
            + str(threshold)
            + "m_snap_merge"
        )
        filename_cw = (
            "CrossWalk_" + 
            supernetwork + 
            "_" + 
            str(threshold) + 
            "m_snap_merge.json"
        )

        print("Snap and merge:")
        print("snapping junctions...")
        data_snapped = snap_junctions(data, threshold, network_data)

        print("merging segments...")
        data_merged, qlat_destinations = segment_merge(
            data, data_snapped, network_data, threshold, pruned_segs
        )

    #-------------------------------
    #      PRUNE - MERGE
    #-------------------------------
    if not snap and prune:
        dirname = (
            "RouteLink_"
            + supernetwork
            + "_"
            + str(threshold)
            + "m_prune_merge"
        )
        filename_cw = (
            "CrossWalk_"
            + supernetwork
            + "_"
            + str(threshold)
            + "m_prune_merge.json"
        )

        print("Prune and merge:")
        print("pruning headwaters...")
        data_pruned = prune_headwaters(data, threshold, network_data)

        # identify pruned segments
        pruned_segs.update(set(np.setdiff1d(data.index, data_pruned.index)))

        print("merging segments...")
        data_merged, qlat_destinations = segment_merge(
            data, data_snapped, network_data, threshold, pruned_segs
        )

    #-------------------------------
    #         MERGE
    #-------------------------------
    if not snap and not prune:
        dirname = (
            "RouteLink_"
            + supernetwork
            + "_"
            + str(threshold)
            + "m_merge"
                  )
        filename_cw = (
            "CrossWalk_"
            + supernetwork
            + "_"
            + str(threshold)
            + "m_merge.json"
        )
        
        print("Just merge:")
        print("merging segments...")
        data_merged, qlat_destinations = segment_merge(
            data, data, network_data, threshold, pruned_segs
        )

    
    if write_output:
        
        # write a new *edited* RouteLink DataFrame
        RouteLink_edit = RouteLink.loc[data_merged.index.values]
        # replace parameters
        for (columnName, columnData) in data_merged.iteritems():
            RouteLink_edit.loc[:, columnName] = columnData
        # replace "to"
        for idx in RouteLink_edit.index:
            if RouteLink_edit.loc[idx, "to"] < 0:
                RouteLink_edit.loc[idx, "to"] = 0

        # export merged data
        print(
            "exporting RouteLink DataFrame:",
            dirname + ".pkl",
            "to",
            os.path.join(root, "test", "input", "geo", "Channels", dirname),
        )

        dir_path = os.path.join(root, "test", "input", "geo", "Channels", dirname)
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

        # save RouteLink data as pickle
        pd.DataFrame(RouteLink_edit).to_pickle(os.path.join(root, "test", "input", "geo", "Channels", dirname, dirname + ".pkl"))

        # save cross walk as json
        print(
            "exporting CrossWalk file:",
            filename_cw,
            "to",
            os.path.join(root, "test", "input", "geo", "Channels", dirname),
        )
        with open(os.path.join(dir_path, filename_cw), "w") as outfile:
            json.dump(qlat_destinations, outfile)

        # export original data
        if return_original:
            RouteLink.loc[data.index.values].to_pickle(os.path.join(root, "test", "input", "geo", "Channels", dirname, supernetwork + ".pkl"))


if __name__ == "__main__":
    main()