import recursive_print
import logging 

LOG = logging.getLogger(__name__)

def get_down_connections(
    rows,
    key_col,
    downstream_col,
    length_col,
    mask_set=None,
    length_key=r"length",
    upstreams_key=r"upstreams",
    downstream_key=r"downstream",
    data_key=r"data",
    verbose=False,
    debuglevel=0,
    LOG=LOG,
):
    # TODO: Consider moving debug and verbose prints to the calling function
    if debuglevel <= -100:
        breakpoint()
    if verbose:
        LOG.info("down connections ...")

    connections = {
        row[key_col]: {
            downstream_key: row[downstream_col],
            length_key: row[length_col],
            data_key: list(row),
        }
        for row in rows
        if row[key_col] in mask_set
    }

    if debuglevel <= -1:
        LOG.debug(f"found {len(connections.keys())} segments")
    if debuglevel <= -3:
        if verbose:
            LOG.info(f"The complete 'connections' object is as follows:")
        LOG.critical(connections)
    if verbose:
        LOG.info("down_connections complete")

    return connections


def get_waterbody_segments(
    connections=None,
    terminal_code=-999,
    waterbody_col=3,
    waterbody_null_code=0,
    data_key=r"data",
    downstream_key=r"downstream",
    upstreams_key=r"upstreams",
    verbose=False,
    debuglevel=0,
    LOG=LOG,    
):

    if verbose:
        LOG.info("level_pool_waterbody_set ...")
    waterbody_dict = {}
    level_pool_waterbody_set = {
        con[data_key][waterbody_col] for key, con in connections.items()
    }
    level_pool_waterbody_set.discard(waterbody_null_code)
    waterbody_dict["level_pool"] = level_pool_waterbody_set
    if debuglevel <= -1:
        LOG.debug(f"found {len(level_pool_waterbody_set)} waterbodies")
    if debuglevel <= -3:
        LOG.critical(level_pool_waterbody_set)
    if verbose:
        LOG.info("level_pool_waterbody_set complete")

    if verbose:
        LOG.info("waterbody segments ...")
    waterbody_segments = {
        key: con[data_key][waterbody_col]
        for key, con in connections.items()
        if not (con[data_key][waterbody_col] == waterbody_null_code)
    }
    if debuglevel <= -1:
        LOG.debug(f"found {len(waterbody_segments)} segments that are part of a waterbody")
    if debuglevel <= -3:
        LOG.critical(waterbody_segments)
    if verbose:
        LOG.info("waterbody_segments complete")

    if verbose:
        LOG.info("waterbody_outlet_set ...")
    waterbody_outlet_set = set()
    for waterbody_segment in waterbody_segments:
        if connections[waterbody_segment][downstream_key] not in waterbody_segments:
            waterbody_outlet_set.add(waterbody_segment)
    if debuglevel <= -1:
        LOG.debug(
            f"found {len(waterbody_outlet_set)} segments that are outlets of a waterbody"
        )
    if debuglevel <= -3:
        LOG.critical(waterbody_outlet_set)
    if verbose:
        LOG.info("waterbody_outlet_set complete")

    if verbose:
        LOG.info("waterbody_downstream_set ...")
    waterbody_downstream_set = set()
    for outlet_segment in waterbody_outlet_set:
        waterbody_downstream_set.add(connections[outlet_segment][downstream_key])
    if debuglevel <= -1:
        LOG.debug(
            f"found {len(waterbody_downstream_set)} segments that are below outlets of a waterbody"
        )
    if debuglevel <= -3:
        LOG.critical(waterbody_downstream_set)
    if verbose:
        LOG.info("waterbody_downstream_set complete")

    if verbose:
        LOG.info("waterbody_upstreams_set ...")
    waterbody_upstreams_set = set()
    for waterbody_segment in waterbody_segments:
        for upstream in connections[waterbody_segment][upstreams_key]:
            if not upstream == terminal_code and not upstream in waterbody_segments:
                waterbody_upstreams_set.add(upstream)
    waterbody_upstreams_set.discard(
        terminal_code
    )  # TODO: Is this the best place for this filter -- check if ever used.
    if debuglevel <= -1:
        LOG.debug(
            f"found {len(waterbody_upstreams_set)} segments that are upstream of a waterbody"
        )
    if debuglevel <= -3:
        LOG.critical(waterbody_upstreams_set)
    if verbose:
        LOG.info("waterbody_upstreams_set complete")

    return (
        waterbody_dict,
        waterbody_segments,
        waterbody_outlet_set,
        waterbody_upstreams_set,
        waterbody_downstream_set,
    )


def determine_keys(
    connections
    # , key_col, downstream_col
    ,
    terminal_code,
    upstreams_key=r"upstreams",
    downstream_key=r"downstream",
    verbose=False,
    debuglevel=0,
    LOG=LOG,
):

    if verbose:
        LOG.info("ref_keys ...")
    ref_keys = {con[downstream_key] for key, con in connections.items()}
    if debuglevel <= -1:
        LOG.debug(f"found {len(ref_keys)} ref_keys")
    if debuglevel <= -3:
        LOG.critical(ref_keys)
    if verbose:
        LOG.info("ref_keys complete")

    if verbose:
        LOG.info("headwater_keys ...")
    headwater_keys = {x for x in connections.keys() if x not in ref_keys}
    if debuglevel <= -1:
        LOG.debug(f"found {len(headwater_keys)} headwater segments")
    if debuglevel <= -3:
        LOG.critical(headwater_keys)
    if verbose:
        LOG.info("headwater_keys complete")

    # Get the downstream terminating nodes
    if verbose:
        LOG.info("terminal_keys ...")
    # Find the pointing-to keys not found in the key dataset.
    terminal_ref_keys = {x for x in ref_keys if x not in connections.keys()}

    # Then collect the keys associated with those 'pointing-tos'
    terminal_keys = set()
    for key, con in connections.items():
        curr_term_ref_key = con[downstream_key]
        if curr_term_ref_key in terminal_ref_keys:
            if curr_term_ref_key != terminal_code:
                if debuglevel <= -2:
                    LOG.warning(
                        f"Non-standard terminal key {con[downstream_key]} found in segment {key}"
                    )
            elif curr_term_ref_key == terminal_code:
                if debuglevel <= -3:
                    print(
                        f"Standard terminal key {con[downstream_key]} found in segment {key}"
                    )
            terminal_keys.add(key)
    if debuglevel <= -1:
        LOG.debug(f"found {len(terminal_keys)} terminal segments")
    if debuglevel <= -1:
        LOG.debug(
            f"of those, {len([x for x in terminal_ref_keys if x != terminal_code])} had non-standard terminal keys"
        )
    if debuglevel <= -3:
        LOG.critical(terminal_keys)
    if verbose:
        LOG.info("terminal_keys complete")

    if verbose:
        LOG.info("circular_keys ...")
    circular_keys = set()
    for key, value in connections.items():
        try:
            # TODO: benchmark try/except vs. nested if statment on 'in' to handle terminal keys
            # e.g., "if key not in terminal_keys: ... etc.
            if connections[connections[key][downstream_key]][downstream_key] == key:
                circular_keys.add(key)
            elif (
                connections[
                    connections[connections[key][downstream_key]][downstream_key]
                ][downstream_key]
                == key
            ):
                circular_keys.add(key)
            elif (
                connections[
                    connections[
                        connections[connections[key][downstream_key]][downstream_key]
                    ][downstream_key]
                ][downstream_key]
                == key
            ):
                circular_keys.add(key)
            elif (
                connections[
                    connections[
                        connections[
                            connections[connections[key][downstream_key]][
                                downstream_key
                            ]
                        ][downstream_key]
                    ][downstream_key]
                ][downstream_key]
                == key
            ):
                circular_keys.add(key)
        except:
            pass

    if debuglevel <= -1:
        LOG.debug(
            f"identified at least {len(circular_keys)} segments with circular references testing to the fourth level"
        )
    if debuglevel <= -3:
        LOG.critical(circular_keys)
    if verbose:
        LOG.info("circular_keys complete")

    return (
        connections.keys(),
        ref_keys,
        headwater_keys,
        terminal_keys,
        terminal_ref_keys,
        circular_keys,
    )


def get_up_connections(
    connections,
    terminal_code,
    headwater_keys,
    terminal_keys,
    upstreams_key=r"upstreams",
    downstream_key=r"downstream",
    verbose=False,
    debuglevel=0,
    LOG=LOG,
):

    # Create inverse of connections looking upstream
    if verbose:
        LOG.info("identifying upstream connections and junctions ...")

    # Using Sets for Junction and Visited keys is REALLY, REALLY, REALLY, FAST!!!
    junction_keys = set()
    visited_keys = set()
    visited_terminal_keys = set()
    junction_count = 0
    for hkey in headwater_keys:
        # TODO: Create a dictionary key identifying relationship to the terminal segment.

        # Start with the headwater keys and label the upstream connections
        # with the terminal_code...
        connections[hkey].update({upstreams_key: {terminal_code}})
        visited_keys.add(hkey)
        # Then iterate through the list and search for the other values
        ukey = hkey
        # print(ukey, hkey)
        # print(visited_keys)
        # print(ukey not in terminal_keys)
        # print(ukey not in junction_keys)
        while True:
            dkey = connections[ukey][downstream_key]
            if (ukey in terminal_keys) or (ukey in junction_keys):
                # If we have hit the bottom (a terminal_key) or if
                # we have joined into an already explored branch, STOP.
                if ukey in terminal_keys:
                    visited_terminal_keys.add(ukey)
                break
            if (
                upstreams_key not in connections[dkey]
            ):  # Check for key in dictionary https://able.bio/rhett/check-if-a-key-exists-in-a-python-dictionary--73iajoz
                connections[dkey].update({upstreams_key: set()})
                connections[dkey][upstreams_key].add(ukey)
                visited_keys.add(dkey)
            else:
                if terminal_code in connections[dkey][upstreams_key]:
                    # If the downstream node here is labeled as a headwater (because it
                    # has an upstream set with the terminal code), it means
                    # that the network had a break and that the traversal has
                    # spanned the gap and the headwater is not actually not a terminating node.
                    # In that case, reset the node to be a blank list (or set, if using
                    # that method), then proceed downstream.
                    # TODO: THIS IS A DANGEROUS/FRAGILE STEP AND DESERVES ADDITIONAL REVIEW
                    # TODO: TO MAKE SURE IT IS DOING WHAT WE INTEND AS DESCRIBED ABOVE
                    # TODO: RESERVOIRS: For instance, this will probably break for subnetworks containing reservoirs

                    connections[dkey].update({upstreams_key: set()})

                connections[dkey][upstreams_key].add(ukey)
                visited_keys.add(dkey)
                # print(dkey, connections[dkey][upstreams_key], visited_keys)
                if len(connections[dkey][upstreams_key]) == 2:
                    if dkey not in junction_keys:
                        junction_keys.add(dkey)
                        junction_count += 1
                    if debuglevel <= -2:
                        LOG.warning(
                            f"Junction found above/into Segment {dkey} with upstream Segments {connections[dkey][upstreams_key]}"
                        )
                elif len(connections[dkey][upstreams_key]) > 2:
                    if dkey not in junction_keys:
                        # At this point, the logic does not allow for this to be a non-junction
                        # TODO: raise/handle error/warning
                        LOG.warning(
                            "key error -- junction analysis has an undetermined anomaly!"
                        )
                    #                         print(dkey in visited_keys)
                    #                         for temp_ukey in connections[dkey][upstreams_key]:
                    #                             print(temp_ukey, temp_ukey in visited_keys)
                    if debuglevel <= -2:
                        LOG.warning(
                            f"revisited Junction above/into Segment {dkey} now with upstream Segments {connections[dkey][upstreams_key]}"
                        )
                    junction_count += 1
            ukey = dkey

    if debuglevel <= -1:
        LOG.debug(f"visited {len(visited_keys)} segments")
    if debuglevel <= -1:
        LOG.debug(
            f"found {junction_count} junctions in {len(junction_keys)} junction nodes"
        )
    if debuglevel <= -3:
        LOG.critical(junction_keys)
    # if debuglevel <= -4:
    #     print(connections)
    if verbose:
        LOG.info("up_connections complete")
    if verbose:
        LOG.info("")

    if verbose:
        LOG.info("confluence segments ...")
    confluence_segment_set = {
        seg for seg, con in connections.items() if con[downstream_key] in junction_keys
    }
    if debuglevel <= -1:
        LOG.debug(f"found {len(confluence_segment_set)} confluence segments")
    if debuglevel <= -3:
        LOG.critical(confluence_segment_set)
    if verbose:
        LOG.info("confluence_segment_set complete")

    return (
        junction_keys,
        confluence_segment_set,
        visited_keys,
        visited_terminal_keys,
        junction_count,
    )


def main():
    """##TEST"""
    print("")
    print("Executing Test")
    # Test data
    debuglevel = -3
    verbose = True
    test_rows = [
        [50, 178, 51, 0],
        [51, 178, 50, 0],
        [60, 178, 61, 0],
        [61, 178, 62, 0],
        [62, 178, 60, 0],
        [70, 178, 71, 0],
        [71, 178, 72, 0],
        [72, 178, 73, 0],
        [73, 178, 70, 0],
        [80, 178, 81, 0],
        [81, 178, 82, 0],
        [82, 178, 83, 0],
        [83, 178, 84, 0],
        [84, 178, 80, 0],
        [0, 456, -999, 0],
        [1, 178, 4, 0],
        [2, 394, 0, 0],
        [3, 301, 2, 0],
        [4, 798, 0, 0],
        [5, 679, 4, 0],
        [6, 523, 0, 0],
        [7, 815, 2, 0],
        [8, 841, -999, 0],
        [9, 514, 8, 0],
        [10, 458, 9, 0],
        [11, 832, 10, 0],
        [12, 543, 11, 0],
        [13, 240, 12, 0],
        [14, 548, 13, 0],
        [15, 920, 14, 0],
        [16, 920, 15, 401],
        [17, 514, 16, 401],
        [18, 458, 17, 0],
        [180, 458, 17, 0],
        [181, 458, 180, 0],
        [19, 832, 18, 0],
        [20, 543, 19, 0],
        [21, 240, 16, 401],
        [22, 548, 21, 0],
        [23, 920, 22, 0],
        [24, 240, 23, 0],
        [25, 548, 12, 0],
        [26, 920, 25, 0],
        [27, 920, 26, 0],
        [28, 920, 27, 0],
    ]

    test_key_col = 0
    test_downstream_col = 2
    test_length_col = 1
    test_terminal_code = -999
    test_waterbody_col = 3
    test_waterbody_null_code = 0

    (test_connections) = get_down_connections(
        rows=test_rows,
        key_col=test_key_col,
        mask_set={row[test_key_col] for row in test_rows},
        downstream_col=test_downstream_col,
        length_col=test_length_col,
        verbose=verbose,
        debuglevel=debuglevel,
    )

    (
        test_all_keys,
        test_ref_keys,
        test_headwater_keys,
        test_terminal_keys,
        test_terminal_ref_keys,
        test_circular_keys,
    ) = determine_keys(
        connections=test_connections,
        terminal_code=test_terminal_code,
        verbose=verbose,
        debuglevel=debuglevel,
    )

    (
        test_junction_keys,
        test_confluence_segment_set,
        test_visited_keys,
        test_visited_terminal_keys,
        test_junction_count,
    ) = get_up_connections(
        connections=test_connections,
        terminal_code=test_terminal_code,
        headwater_keys=test_headwater_keys,
        terminal_keys=test_terminal_keys,
        verbose=verbose,
        debuglevel=debuglevel,
    )

    # TODO: Set/pass/identify a proper flag value
    if test_waterbody_col is not None:
        (
            test_waterbody_dict,
            test_waterbody_segments,
            test_waterbody_outlet_set,
            test_waterbody_upstreams_set,
            test_waterbody_downstream_set,
        ) = get_waterbody_segments(
            connections=test_connections,
            terminal_code=test_terminal_code,
            waterbody_col=test_waterbody_col,
            waterbody_null_code=test_waterbody_null_code,
            verbose=verbose,
            debuglevel=debuglevel,
        )

    recursive_print.print_connections(
        headwater_keys=test_headwater_keys,
        down_connections=test_connections,
        up_connections=test_connections,
        terminal_code=test_terminal_code,
        terminal_keys=test_terminal_keys,
        terminal_ref_keys=test_terminal_ref_keys,
        debuglevel=debuglevel,
    )

    recursive_print.print_basic_network_info(
        connections=test_connections,
        headwater_keys=test_headwater_keys,
        junction_keys=test_junction_keys,
        terminal_keys=test_terminal_keys,
        terminal_code=test_terminal_code,
        verbose=verbose,
        debuglevel=debuglevel,
    )


if __name__ == "__main__":
    main()
