import networkbuilder

def recursive_junction_read (
                             segments
                             , order_iter
                             , connections
                             , network
                             , terminal_code = 0
                             , verbose = False
                             , debuglevel = 0
                            ):

    for segment in segments:
        csegment = segment
        if 1 == 1:
        #try:
            reach = {}
            reach.update({'reach_tail':csegment})
            reach.update({'downstream_reach':connections[csegment]['downstream']})
            segmentset = set()
            segmentlist = [] # Ordered Segment List bottom to top
            usegments = connections[segment]['upstreams']
            while True: 
                if usegments == {terminal_code}: # HEADWATERS
                    if debuglevel <= -3: print(f"headwater found at {csegment}")
                    network['total_segment_count'] += 1
                    if debuglevel <= -3: print(f"segs at csegment {csegment}: {network['total_segment_count']}")
                    segmentset.add(csegment)
                    segmentlist.append(csegment) # Ordered Segment List bottom to top
                    reach.update({'reach_head':csegment})
                    reach.update({'seqorder':order_iter})
                    if order_iter == 0: network.update({'terminal_reach':csegment})#; import pdb; pdb.set_trace() #TODO: FIX THIS; SEEMS FRAGILE
                    network.update({'maximum_reach_seqorder':max(network['maximum_reach_seqorder'],order_iter)})
                    reach.update({'segments':segmentset})
                    reach.update({'segments_list':segmentlist}) # Ordered Segment List bottom to top
                    network['reaches'].update({csegment:reach})
                    network['headwater_reaches'].add(csegment)
                    break
                elif len(usegments) >= 2: # JUNCTIONS
                    if debuglevel <= -3: print(f"junction found at {csegment} with upstreams {usegments}")
                    network['total_segment_count'] += 1
                    if debuglevel <= -3: print(f"segs at csegment {csegment}: {network['total_segment_count']}")
                    segmentset.add(csegment)
                    segmentlist.append(csegment) # Ordered Segment List
                    reach.update({'reach_head':csegment})
                    reach.update({'seqorder':order_iter})
                    if order_iter == 0: network.update({'terminal_reach':csegment})#; import pdb; pdb.set_trace() #TODO: FIX THIS; SEEMS FRAGILE
                    network.update({'maximum_reach_seqorder':max(network['maximum_reach_seqorder'],order_iter)})
                    reach.update({'segments':segmentset})
                    reach.update({'segments_list':segmentlist}) # Ordered Segment List bottom to top
                    network['reaches'].update({csegment:reach})
                    network['total_junction_count'] += 1 #the Terminal Segment
                    network['junctions'].add(csegment)
                    recursive_junction_read (
                            usegments
                            , order_iter + 1
                            , connections
                            , network
                            , terminal_code = terminal_code
                            , verbose = verbose
                            , debuglevel = debuglevel) 
                    break
                network['total_segment_count'] += 1
                if debuglevel <= -3: print(f"segs at csegment {csegment}: {network['total_segment_count']}")
                # the terminal code will indicate a headwater
                if debuglevel <= -4: print(usegments)
                segmentset.add(csegment)
                segmentlist.append(csegment) # Ordered Segment List
                (csegment,) = usegments
                usegments = connections[csegment]['upstreams']

                # print(usegments)
        #except:
            #if debuglevel <= -2: 
                #print(f'There is a problem with connection: {segment}: {connections[segment]}')

#TODO: re-implement as parallel process on terminal segments
#TODO: NOTE: this will be a little more complicated now that the global `connections` 
#TODO: NOTE: object is not available from the calling function. Probably requires a class.
def network_trace(
        terminal_segment
        , order_iter
        , connections
        , terminal_code = 0
        , verbose= False
        , debuglevel = 0
    ):

    network = {}
    us_length_total = 0
    
    if debuglevel <= -2: print(f'\ntraversing upstream on network {terminal_segment}:')
    # try:
    if 1 == 1:
        network.update({'total_segment_count': 0}) 
        network.update({'total_junction_count': 0})
        network.update({'maximum_reach_seqorder':0})
        network.update({'junctions':set()})
        network.update({'headwater_reaches':set()})
        network.update({'reaches':{}}) 
        recursive_junction_read(
                  [terminal_segment]
                  , order_iter
                  , connections
                  , network
                  , verbose = verbose
                  , terminal_code = terminal_code
                  , debuglevel = debuglevel)
        if verbose: print(f"junctions: {network['total_junction_count']}")
        if verbose: print(f"segments: {network['total_segment_count']}")
    # except Exception as exc:
    #     print(exc)
    #TODO: compute upstream length as a surrogate for the routing computation
    return {terminal_segment: network, 'upstream_length': us_length_total}

def compose_networks(
        supernetwork_values = None
        , terminal_code = 0
        , debuglevel = 0
        , verbose = False
        , showtiming = False
    ):

    terminal_segments = supernetwork_values[4] 
    circular_segments = supernetwork_values[6]
    terminal_segments_super = terminal_segments - circular_segments
    connections = supernetwork_values[0]
        
    networks = {terminal_segment:{}
                      for terminal_segment in terminal_segments_super 
                     }  

    if verbose: print('verbose output')
    if verbose: print(f'number of Independent Networks to be analyzed is {len(networks)}')
    if verbose: print(f'Multi-processing will use {multiprocessing.cpu_count()} CPUs')
    if verbose: print(f'debuglevel is {debuglevel}')

    init_order = 0
    for terminal_segment, network in networks.items():
        network.update(network_trace(terminal_segment = terminal_segment
            , order_iter = init_order
            , connections = connections
            , terminal_code = terminal_code
            , verbose = verbose
            , debuglevel = debuglevel)[terminal_segment])
        up_reaches = networkbuilder.get_up_connections(
            network['reaches']
            , terminal_code
            , network['headwater_reaches']
            , {network['terminal_reach']}
            , r'upstream_reaches'
            , r'downstream_reach'
            , verbose = verbose
            , debuglevel = debuglevel
            )

        if debuglevel <= -2:
            if debuglevel <=-2: print(f'terminal_segment: {terminal_segment}')
            if debuglevel <=-3: 
                for k, v in network.items():
                    if type(v) is dict:
                        print (f'{k}:')
                        for k1, v1 in v.items():
                            print(f'\t{k1}: {v1}')
                    else: print(f'{k}: {v}')
    if debuglevel <= -1: print(f'Number of networks in the Supernetwork: {len(networks.items())}')

    return networks

