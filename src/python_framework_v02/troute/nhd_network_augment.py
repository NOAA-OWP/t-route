import nhd_network
import numpy as np
from functools import partial


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
    conn = nhd_network.extract_connections(data, network_data["downstream_col"])
    
    # extract upstream connections
    rconn = nhd_network.reverse_network(conn)
    
    return conn, rconn

def build_reaches(rconn):
    
    # isolate independent subnetworks
    subnets = nhd_network.reachable_network(rconn)

    # identify the segments in each subnetwork
    subreachable = nhd_network.reachable(rconn)

    # break each subnetwork into reaches
    subreaches = {}
    for tw, net in subnets.items():
        path_func = partial(nhd_network.split_at_junction, net)
        subreaches[tw] = nhd_network.dfs_decomposition(net, path_func)
        
    return subreachable, subreaches, subnets

def headwater_connections(data, network_data):

    """
    Determine which segments are and are not headwaters. 
    Headwaters are defined as segments with no upstream connection, only downstream connections. 
    Non-headwaters are defined as segments with both upstream and downstream connections. 
    Args:
        data (DataFrame): Network parameter dataset, prepared
    Returns:
        hw_conn (dict): downstream connections, headwaters only
        non_hw_conn (dict): downstream connections, non headwaters only
    """

    # extract network connections
    conn, rconn = network_connections(data, network_data)
    
    hw = []  # store headwater segments
    non_hw = []  # store non-headwater segments

    for seg in rconn.keys():
        # if there is no upstream connection anda downstream connection, it is a headwater
        if bool(rconn[seg]) == False and bool(conn[seg]) == True:
            hw.append(seg)

        # if there is an upstream connection and a downstream connection, it is a non-headwater (midwater?)
        elif bool(rconn[seg]) == True and bool(conn[seg]) == True:
            non_hw.append(seg)

    # get segment key-value pairs from the connections dictionary
    hw_conn = {key: conn[key] for key in hw}
    non_hw_conn = {key: conn[key] for key in non_hw}

    return hw_conn, non_hw_conn

def prune_headwaters(hw_connections, non_hw_connections, data, thresh):

    # TO DO: Bug fix - why does the function return non-unique headwater IDs in the chop list? 
    #        some headwaters are being identified for chopping more than once?
    
    """
    Prune headwaters from the network
    Headwaters are pruned if they are shorter than the threshold length
    OR if they merge with a midwater that is less than the threshold length
    Args:
        hw_connections (dict): downstream connections, headwaters only
        non_hw_connections (dict): downstream connections, midwaters only
        data (DataFrame): Network to be pruned
        thresh (int): length threshold, segments below this length are slated for pruning
    Returns:
        chop (list): list of pruned headwater indices, updated
        data_pruned (DataFrame): pruned network
    """

    # initialize a list of pruned headwater segments
    chop = []
    
    # loop through keys and values in headwater connections dictionary
    for hw_k, srch_v in hw_connections.items():

        """
        Trim headwater-headwater junctions
        
        - Headwater-headwater junctions are junctions where two (or more) headwaters merge. 
        - Once headwater-headwater junctions are identified, the shortest headwater below the
          threshold length is pruned. 
        """

        # create a list of headwaters draining to the specified downstream connections, srch_v
        G = [k for k, v in hw_connections.items() if v == srch_v]

        # if more than one headwater drains into the segment, it is a headwater-headwater junction
        if len(G) > 1:

            # get segment lengths from the test dataset
            hw_len = data.loc[G].Length

            # find the shortest segment
            hw_min_len = hw_len[hw_len.idxmin()]

            # check to see if headwater shorter than threshold length
            if hw_min_len < thresh:

                # update chopping block
                chop.append(hw_len.idxmin())

        ############################################################
        #  Trim headwater-midwater junctions
        ############################################################

        """
        Trim headwater-midwater junctions
        
        - Headwater-midwater junctions are junctions a headwater merges with a non-headwater (i.e. midwater). 
        - If the headwater is shorter than the threshold, it is pruned.
        - If the headwater is longer than the threshold, but joins with or connects to a midwater that is
          shorter than the threshold, the headwater is pruned. 
              - This is done to minimize the number of small midwater segments stranded between headwaters.
                Else, these short, stranded segments cannot be merged upstream or downstream. 
        """

        # midwaters draining into a junction that a headwater does too
        H = [k for k, v in non_hw_connections.items() if v == srch_v]

        # if there is a corresponding midwater, then the headwater is a candidate for trimming
        if bool(H) == True:

            # get segment lengths from the dataset
            hw_len = data.loc[hw_k].Length

            # get midwater (in) length from the test dataset
            mw_in_len = data.loc[H[0]].Length

            # get midwater (out) length from the test dataset
            mw_out_len = data.loc[srch_v[0]].Length

            # check to see if headwater or midwaters (in/out) are shorter than the threshold
            if hw_len < thresh or mw_in_len < thresh or mw_out_len < thresh:

                # update chopping block
                chop.append(hw_k)

    # chop segments from the network
    data_pruned = data.drop(np.unique(chop))

    return np.unique(chop), data_pruned

def len_weighted_av(df, var, weight):

    """
    Calculate a weighted average
    Args:
        df (DataFrame): DataFrame containing variables to be averaged and used as weights
        var (str): name of the variable to be averaged
        weight (str): name of the variable to be used as a weight
    Returns:
        x (float32): weighted average
    """

    x = (df[weight] * df[var]).sum() / df[weight].sum()

    return x

def merge_parameters(to_merge):

    """
    length-weighted averaging of channel routing parameters across merged segments
    Args:
        to_merge (DataFrame): DataFrame containing routing parameters for segments to be merged together
    Returns:
        replace (DataFrame): weighted average
    """

    data_replace = to_merge.tail(1)

    idx = to_merge.tail(1).index
    data_replace.loc[idx, "Length"] = to_merge.Length.sum()
    data_replace.loc[idx, "n"] = len_weighted_av(to_merge, "n", "Length")
    data_replace.loc[idx, "nCC"] = len_weighted_av(to_merge, "nCC", "Length")
    data_replace.loc[idx, "So"] = len_weighted_av(to_merge, "So", "Length")
    data_replace.loc[idx, "BtmWdth"] = len_weighted_av(to_merge, "BtmWdth", "Length")
    data_replace.loc[idx, "TopWdth"] = len_weighted_av(to_merge, "TopWdth", "Length")
    data_replace.loc[idx, "TopWdthCC"] = len_weighted_av(
        to_merge, "TopWdthCC", "Length"
    )
    data_replace.loc[idx, "MusK"] = len_weighted_av(to_merge, "MusK", "Length")
    data_replace.loc[idx, "MusX"] = len_weighted_av(to_merge, "MusX", "Length")
    data_replace.loc[idx, "ChSlp"] = len_weighted_av(to_merge, "ChSlp", "Length")

    return data_replace

def correct_reach_connections(data_merged):

    """
    Update downstream connections ("to") for segments in a merged reach.
    Only updates *in-reach* connections.
    Args:
        data_merged (DataFrame): Routing parameters for segments in merged reach
    Returns:
        data_merged (DataFrame): Routing parameters for segments in merged reach with updated donwstream connections
    """

    for i, idx in enumerate(data_merged.index.values[0:-1]):
        data_merged.loc[idx, "to"] = data_merged.index.values[i + 1]

    return data_merged

def upstream_merge(data_merged, chop):

    """
    Merge a short reach tail segment with upstream neighbor
    Args:
        data_merged (DataFrame): Routing parameters for segments in merged reach
        chop (list): list of merged-out segments
    Returns:
        data_merged (DataFrame): Routing parameters for segments in merged reach with updated donwstream connections
        chop (list): updated list of merged-out segments
    """

    # grab the two segments that need to be merged - simply the last two segments of the reach
    to_merge = data_merged.tail(2)

    # calculate new parameter values
    data_replace = merge_parameters(to_merge)

    # paste new parameters in to data_merged
    data_merged.loc[to_merge.tail(1).index] = data_replace

    # remove merged segments from data_merged
    data_merged = data_merged.drop(to_merge.head(1).index)

    # update "chop" list with merged-out segment IDs
    chop.append(to_merge.head(1).index.values[0])

    return data_merged, chop

def downstream_merge(data_merged, chop, thresh):

    """
    Merge short segments with their downstream neighbors
    Args:
        data_merged (DataFrame): Routing parameters for segments in merged reach
        chop (list): list of merged-out segments
        thresh (int): theshold reach length (meters)
    Returns:
        data_merged (DataFrame): Routing parameters for segments in merged reach with updated donwstream connections
        chop (list): updated list of merged-out segments
    """

    # find the upstream-most short segment and it's downstream connection
    idx_us = data_merged.loc[data_merged.Length < thresh].head(1).index.values[0]

    pos_idx_us = data_merged.index.get_loc(idx_us)
    idx_to = data_merged.iloc[pos_idx_us + 1].name

    # grab segments to be merged
    to_merge = data_merged.loc[[idx_us, idx_to]]

    # calculate new parameter values
    data_replace = merge_parameters(to_merge)

    # paste new parameters in to data_merged
    data_merged.loc[to_merge.tail(1).index] = data_replace

    # remove merged segments from data_merged
    data_merged = data_merged.drop(to_merge.head(1).index)

    # update "chop" list with merged-out segment IDs
    chop.append(to_merge.head(1).index.values[0])

    return data_merged, chop

def merge_all(rch, data, chop):

    """
    Merge all segments in a reach
    Args:
        rch (list): Segment indices in the reach to be merged
        data (DataFrame): Routing parameters for network containing the reach to be merged
        chop (list): list of merged-out segments
    Returns:
        data_merged (DataFrame): Routing parameters for segments in merged reach with updated donwstream connections
        chop (list): updated list of merged-out segments
    """

    # subset the model parameter data for this reach
    data_merged = data.loc[rch].copy()

    # grab the two segments that need to be merged - in this case, merge all segments!
    to_merge = data_merged.copy()

    # calculate new parameter values
    data_replace = merge_parameters(to_merge)

    # paste new parameters in to data_merged
    data_merged.loc[to_merge.tail(1).index] = data_replace

    # remove merged segments from data_merged - in this case, all but the last
    data_merged = data_merged.drop(data_merged.iloc[:-1, :].index)

    # update "chop" list with merged-out segment IDs
    chop.extend(list(to_merge.iloc[:-1, :].index))

    return data_merged, chop

def update_network_data(data, rch, data_merged, chop, rconn):

    """
    Update the network routing parameter data with merged segment data
    Args:
        data (DataFrame): Routing parameters for network to be updated
        rch (list): Segment indices in the reach to be merged
        data_merged (DataFrame): Routing parameters for merged reach
    Returns:
        data (DataFrame): Updated network routing parameters
    """

    # make a copy of the data to be replaced with merged
    data_old = data.loc[rch].copy()

    # drop the segments that disapeared with merger
    data = data.drop(chop)

    # adjust the segment data for those that remain
    data.loc[data_merged.index] = data_merged

    # update out of reach connections - these will change in the first segment was merged out
    upstreams = rconn[
        data_old.head(1).index.values[0]
    ]  # upstream connection of the OLD reach head

    if bool(upstreams):
        
        data.loc[upstreams, "to"] = data_merged.head(1).index.values[
            0
        ]  # index of NEW reach head

    return data

def qlat_destination_compute(data_native, data_merged, merged_segments, pruned_segments, network_data):
    
    if bool(list(pruned_segments)):
        
        segments = merged_segments + list(pruned_segments)
        
    else:
        
        segments = merged_segments

    # compute connections
    conn = nhd_network.extract_connections(data_native, network_data["downstream_col"])
    
    # initialize a library to store qlat desinations for pruned/merged segments
    qlat_destinations = {x: {} for x in segments}

    for idx in segments:

        # get the downstream connection of this segment in the original network
        ds_idx = conn[idx]

        # find the nearest downstream segment remaining in the pruned/merged network 
        while bool(ds_idx[0] in data_merged.index) == False:
            ds_idx = conn[ds_idx[0]]

        # update the qlat destination dict
        qlat_destinations[idx] = ds_idx
        
    return qlat_destinations

def segment_merge(data_native, data, network_data, thresh, pruned_segments):
    
    # create a copy of the pruned network dataset, which will be updated with merged data
    data_merged = data.copy()

    # initialize list to store merged segment IDs
    merged_segments = []
    
    # organize network into reaches
    conn = nhd_network.extract_connections(data, network_data["downstream_col"])
    rconn = nhd_network.reverse_network(conn)
    subnets = nhd_network.reachable_network(rconn)
    subreachable = nhd_network.reachable(rconn)
    
    subreaches = {}
    for tw, net in subnets.items():
        path_func = partial(nhd_network.split_at_junction, net)
        subreaches[tw] = nhd_network.dfs_decomposition(net, path_func)
    
    # loop through each reach in the network
    for twi, (tw, rchs) in enumerate(subreaches.items(), 1):

        for rch in rchs:

            # calculate reach length
            rch_len = data.loc[rch].Length.sum()

            ##################################################
            # orphaned short single segment reaches
            ##################################################
            # if reach length is shorter than threshold and composed of a single segment
            if rch_len < thresh and len(data.loc[rch]) == 1:
                continue  # do nothing

            ##################################################
            # multi segment reaches - combine into a single segment reach
            ##################################################
            # if reach length is shorter than threshold and composed more than one segment
            if rch_len < thresh and len(data.loc[rch]) > 1:

                # merge ALL reach segments into one
                chop = []
                reach_merged, chop = merge_all(rch, data, chop)

                # update network with merged reach data
                data_merged = update_network_data(
                    data_merged, 
                    rch, 
                    reach_merged, 
                    chop, 
                    rconn
                )

                # update merged_segmetns list with merged-out segments
                merged_segments.extend(chop)

            ##################################################
            # multi segment reaches longer than threshold with some segments shorter than threshold
            ##################################################
            # if reach length is longer than threshold and smallest segment length is less than threshold
            if rch_len > thresh and data.loc[rch].Length.min() < thresh:

                # initialize data_merged - this DataFrame will be subsequently revised
                reach_merged = data.loc[rch]

                # initialize list of segments chopped from this reach
                chop_reach = []

                # so long as the shortest segment is shorter than the threshold...
                while reach_merged.Length.min() < thresh:

                    # if shortest segment is the last segment in the reach - conduct an upstream merge.
                    if (
                        reach_merged.Length.idxmin() == reach_merged.tail(1).index.values[0]
                        and reach_merged.Length.min() < thresh
                    ):

                        # upstream merge
                        chop = []
                        reach_merged, chop = upstream_merge(reach_merged, chop)

                        # update chop_reach list with merged-out segments
                        chop_reach.extend(chop)

                    # if shortest segment is NOT the last segment in the reach - conduct a downstream merge
                    if (
                        reach_merged.Length.idxmin() != reach_merged.tail(1).index.values[0]
                        and reach_merged.Length.min() < thresh
                    ):

                        # downstream merge
                        chop = []
                        reach_merged, chop = downstream_merge(reach_merged, chop, thresh)

                        # update chop_reach list with merged-out segments
                        chop_reach.extend(chop)

                # correct segment connections within reach
                reach_merged = correct_reach_connections(reach_merged)

                # update the greater network data set
                data_merged = update_network_data(
                    data_merged,
                    rch,
                    reach_merged,
                    chop_reach,
                    rconn
                )

                # update merged_segmetns list with merged-out segments
                merged_segments.extend(chop_reach)
                
    # create a qlateral destinations dictionary
    qlat_destinations = qlat_destination_compute(data_native,
                                                 data_merged,
                                                 merged_segments,
                                                 pruned_segments, 
                                                 network_data)
    
    return data_merged, qlat_destinations