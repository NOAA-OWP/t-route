def print_basic_network_info(connections, headwater_keys
                            , junction_keys = None
                            , terminal_keys = None
                            , terminal_code = None
                            , verbose = False, debuglevel = 0):
    i = 0
    for key, value in connections.items():
        try:
            i += (len(value['upstreams']) - 1)
            if debuglevel <= -2 and len(value['upstreams']) > 1:
                print (f"junction# {i} at key {key} has these upstream segments {value['upstreams']}")
        except:
            if debuglevel <= -1: print(f"No valid upstream computed for segment {key} -- possible circular reference")
            if debuglevel <= -2: print (f" key {key} points to {value['downstream']} which points to {connections[value['downstream']]['downstream']}")


    print(f'Total Segments {len(connections)}')
    print(f'Head Segments {len(headwater_keys)}')
    print(f'found {i} junctions')
    if junction_keys:
        print(f'...in {len(junction_keys)} junction nodes')
        print(f'Total Reaches ( = head_segments + junction_nodes ) {len(junction_keys) + len(headwater_keys)}')
    if terminal_keys:
        print(f'Total Independent Networks (estimated from number of terminal keys) {len(terminal_keys)}')
    #     if terminal_code is not None:
    #         print(f'(of those, {sum( 1 for key in terminal_keys if connections[key]['downstream'] != terminal_code)})')
    print('\n')

def rec_print_down(key, down_connections, terminal_ref_keys, debuglevel = 0):
    if key in terminal_ref_keys: return
    print(f"{key} with length {down_connections[key]['length']}")
    rec_print_down(down_connections[key]['downstream'], down_connections, terminal_ref_keys, debuglevel)

def rec_print_up(keys, tab_count, up_connections, down_connections
                , terminal_code, debuglevel = 0):
    if not isinstance(keys, set): keys = {keys}
    tab_count += 1
    for key in keys:
        if not key == terminal_code:
        # takes advantage of terminal code assigned as upstream to headwater keys...
            print(f"{'..' * (tab_count)}\\{key} with length {down_connections[key]['length']}\\")
            rec_print_up(up_connections[key]['upstreams'], tab_count
                                , up_connections, down_connections
                                , terminal_code
                                , debuglevel)

def print_connections(headwater_keys = None, terminal_keys = None
                    , down_connections = None, up_connections = None
                    , terminal_code = None, terminal_ref_keys = None
                    , debuglevel = 0):
    try:
        if headwater_keys:
            print("########################")
            print("Downstream Connections")
            print("########################")
            for key in headwater_keys:
                rec_print_down(key, down_connections, terminal_ref_keys, debuglevel)
                print("########################")

        if terminal_keys:
            print("########################")
            print("Upstream Connections")
            print("########################")
            for key in terminal_keys:
                rec_print_up({key}, -1, up_connections, down_connections
                                , terminal_code, debuglevel = debuglevel)
                print("########################")
    except:
        if debuglevel <= -1: print('''Error: Provide headwater_keys, down_connections, and a terminal code
to print Downstream Connections.

Provide terminal_keys, up_connections, down_connections, and a terminal code
to print Upstream Connections.''')
