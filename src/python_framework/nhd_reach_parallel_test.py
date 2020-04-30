## Basic imports
import sys
import os
import time

connections = None

ENV_IS_CL = False
if ENV_IS_CL: root = '/content/wrf_hydro_nwm_public/trunk/NDHMS/dynamic_channel_routing/'
elif not ENV_IS_CL: 
    sys.setrecursionlimit(4000)
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.append(os.path.join(root, r'src', r'python_framework'))

## network and reach utilities
import nhd_network_utilities as nnu
import nhd_reach_utilities as nru

def compute_network_parallel_totaltreedepth(
        large_networks = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections

    overall_max = -1
    for terminal_segment, network in large_networks.items():
        overall_max = max(network['maximum_reach_seqorder'], overall_max)

    if debuglevel <= -1: print(f"Executing simulation for all reaches beginning with maximum order {overall_max}")

    #if 1 == 1:
    ordered_reaches = {}
    ordered_reaches.update({overall_max:[]})
    #import pdb; pdb.set_trace()
    for terminal_segment, network in large_networks.items():
        for head_segment, reach in network['reaches'].items():
            if reach['seqorder'] not in ordered_reaches:
                ordered_reaches.update({reach['seqorder']:[]}) #TODO: Should this be a set/dictionary?
            ordered_reaches[reach['seqorder']].append([head_segment
                      , reach
                      , supernetwork_data
                      , verbose
                      , debuglevel])
            #import pdb; pdb.set_trace()
            #print(ordered_reaches)
    

    if debuglevel <= -1:
        for x in range(overall_max,-1,-1):
            nslist3 = ordered_reaches[x]
            print(f'{x}: {len(nslist3)}')

def compute_network_parallel_totaltreedepth_wHEADS(
        large_networks = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections

    overall_max = -1
    for terminal_segment, network in large_networks.items():
        overall_max = max(network['maximum_reach_seqorder'], overall_max)

    if debuglevel <= -1: print(f"Executing simulation for all reaches beginning with maximum order {overall_max}")

    #if 1 == 1:
    ordered_reaches = {}
    ordered_reaches.update({overall_max:[]})
    #import pdb; pdb.set_trace()
    for terminal_segment, network in large_networks.items():
        for head_segment, reach in network['reaches'].items():
            if reach['seqorder'] not in ordered_reaches:
                ordered_reaches.update({reach['seqorder']:[]}) #TODO: Should this be a set/dictionary?
            if reach['upstream_reaches'] == {supernetwork_data['terminal_code']}:
                ordered_reaches[overall_max].append([head_segment
                          , reach
                          , supernetwork_data
                          , verbose
                          , debuglevel])
            else:
                ordered_reaches[reach['seqorder']].append([head_segment
                          , reach
                          , supernetwork_data
                          , verbose
                          , debuglevel])
            #import pdb; pdb.set_trace()
            #print(ordered_reaches)
    

    if debuglevel <= -1:
        for x in range(overall_max,-1,-1):
            nslist3 = ordered_reaches[x]
            print(f'{x}: {len(nslist3)}')

def compute_network_parallel_HEADS(
        large_networks = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections

    overall_max = -1
    for terminal_segment, network in large_networks.items():
        overall_max = max(network['maximum_reach_seqorder'], overall_max)

    if debuglevel <= -1: print(f"Executing simulation for all large network beginning with maximum order {overall_max}")

    #if 1 == 1:
    ordered_reaches = {}
    ordered_reaches.update({overall_max:[]})
    #import pdb; pdb.set_trace()
    for terminal_segment, network in large_networks.items():
        for head_segment, reach in network['reaches'].items():
            if reach['seqorder'] not in ordered_reaches:
                ordered_reaches.update({reach['seqorder']:[]}) #TODO: Should this be a set/dictionary?
            if reach['upstream_reaches'] == {supernetwork_data['terminal_code']}:
                print('headwater')
                ordered_reaches[overall_max].append([head_segment
                          , reach
                          , supernetwork_data
                          , verbose
                          , debuglevel])
            else:
                ordered_reaches[reach['seqorder']].append([head_segment
                          , reach
                          , supernetwork_data
                          , verbose
                          , debuglevel])
            #import pdb; pdb.set_trace()
            #print(ordered_reaches)
    

    if debuglevel <= -1:
        for x in range(overall_max,-1,-1):
            nslist3 = ordered_reaches[x]
            print(f'{x}: {len(nslist3)}')


def compute_network_parallel_opportunistic(
        large_networks = None
        , supernetwork_data = None
        # , connections = None
        , verbose = False
        , debuglevel = 0
        ):

    global connections

    overall_max = -1
    for terminal_segment, network in large_networks.items():
        overall_max = max(network['maximum_reach_seqorder'], overall_max)

    if debuglevel <= -1: print(f"Executing simulation all reaches beginning with maximum order {overall_max}")

    #if 1 == 1:
    ordered_reaches = {}
    #import pdb; pdb.set_trace()

    handled_reaches = set() # When this set contains the whole collection, we are done
    for terminal_segment, network in large_networks.items():
        network_ordered_reaches = {}
        network_handled_reaches = set() # When this set contains the whole collection, we are done
        network_tail_reaches = set() # When this set contains the whole collection, we are done
        parorder = 0
        if parorder not in network_ordered_reaches: network_ordered_reaches.update({parorder:set()})
        network_ordered_reaches[parorder].update(network['headwater_reaches'] | network['receiving_reaches'])
        network_tail_reaches.update(network_ordered_reaches[parorder])
        network_handled_reaches.update(network_ordered_reaches[parorder])

        parorder = 1
        while not network_handled_reaches == (network['receiving_reaches'] | network['headwater_reaches'] | network['junctions']):
            network_handled_reaches_this_round = set()
            tail_reaches_added_this_round = set()
            if parorder not in network_ordered_reaches: network_ordered_reaches.update({parorder:set()})
            for head_segment in network_tail_reaches:
                down = network['reaches'][head_segment]['downstream_reach']
                if not down == supernetwork_data['terminal_code']:
                    downups = network['reaches'][down]['upstream_reaches']
                    if downups.issubset(network_handled_reaches):
                        network_ordered_reaches[parorder].add(down)
                        network_handled_reaches_this_round.add(head_segment)
                        tail_reaches_added_this_round.add(down)

            network_tail_reaches.difference_update(network_handled_reaches_this_round) 
            network_tail_reaches.update(tail_reaches_added_this_round) 
            network_handled_reaches.update(network_ordered_reaches[parorder])
            parorder += 1

        handled_reaches.update(network_handled_reaches)
        for po, rs in network_ordered_reaches.items():
            if po not in ordered_reaches: ordered_reaches.update({po:set()})
            ordered_reaches[po].update(rs)

    if debuglevel <= -1:
        for po, l in ordered_reaches.items():
            print(f'{po}: {len(l)}')
    
def main():

    global connections

    print('This script demonstrates the parallel traversal of reaches developed from NHD datasets')

    verbose = True
    debuglevel = -2
    showtiming = True
    break_network_at_waterbodies = False

    test_folder = os.path.join(root, r'test')
    geo_input_folder = os.path.join(test_folder, r'input', r'geo')

    #TODO: Make these commandline args
    supernetwork = 'Pocono_TEST1'
    """##NHD Subset (Brazos/Lower Colorado)"""
    # supernetwork = 'Brazos_LowerColorado_ge5'
    """##NHD CONUS order 5 and greater"""
    # supernetwork = 'CONUS_ge5'
    """These are large -- be careful"""
    # supernetwork = 'Mainstems_CONUS'
    # supernetwork = 'CONUS_FULL_RES_v20'
    # supernetwork = 'CONUS_Named_Streams' #create a subset of the full resolution by reading the GNIS field
    # supernetwork = 'CONUS_Named_combined' #process the Named streams through the Full-Res paths to join the many hanging reaches

    if verbose: print('creating supernetwork connections set')
    if showtiming: start_time = time.time()
    #STEP 1
    supernetwork_data, supernetwork_values = nnu.set_networks(
        supernetwork = supernetwork
        , geo_input_folder = geo_input_folder
        , verbose = False
        # , verbose = verbose
        , debuglevel = debuglevel
        )
    if verbose: print('supernetwork connections set complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 2
    if showtiming: start_time = time.time()
    if verbose: print('organizing connections into networks and reaches ...')
    networks = nru.compose_networks(
        supernetwork_values
        , break_network_at_waterbodies = break_network_at_waterbodies
        , debuglevel = debuglevel
        # , verbose = False
        , verbose = verbose
        , showtiming = showtiming
        )
    if verbose: print('reach organization complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 3
    if verbose: print(f'Now computing the reaches in parallel')
    if verbose: print(f'This is just a DUMMY computation')
    if verbose: print(f'Only the number of potentially parallelizable reaches is shown for each order')

    connections = supernetwork_values[0]
    #initialize flowdepthvel dict
    parallel_split = -1 # -1 turns off the splitting and runs everything through the lumped execution

    ##STEP 3_ -- Small Networks
    # Another parallelization method is to simply execute a network or group of networks independent from the others. 
    ##TODO:Come back to this -- Essentially, we isolated each network and found that parallel speedup was minimal. 
    #TODO: add the Parallel Split to show network splitting
    ##Each of the subgroups could use one of the three parallel methods or simply execute 
    ##serially for itself (with other subgoups running at the same time).
    ##Probably with some effort, this could be effective.
    if parallel_split >= 0: print(r'DO NOT RUN WITH `parallel_split` >= 0')


    #STEP 3a -- Large Networks by total tree depth
    if showtiming: start_time = time.time()
    if verbose: print(f'executing computation on reaches ordered by distance from terminal node')

    large_networks = {terminal_segment: network \
                      for terminal_segment, network in networks.items() \
                      if network['maximum_reach_seqorder'] > parallel_split}
    # print(large_networks)
    compute_network_parallel_totaltreedepth(
        large_networks
        , supernetwork_data = supernetwork_data
        # , connections = connections
        # , verbose = False
        , verbose = verbose
        , debuglevel = debuglevel
    )
    if verbose: print(f'ordered reach traversal complete for reaches ordered by distance from terminal node')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    #STEP 3b -- Large Networks by total tree depth -- segregating the headwaters
    if showtiming: start_time = time.time()
    if verbose: print(f'executing computation for all headwaters, then by reaches ordered by distance from terminal node')

    large_networks = {terminal_segment: network \
                      for terminal_segment, network in networks.items() \
                      if network['maximum_reach_seqorder'] > parallel_split}
    # print(large_networks)
    compute_network_parallel_totaltreedepth_wHEADS(
        large_networks
        , supernetwork_data = supernetwork_data
        # , connections = connections
        # , verbose = False
        , verbose = verbose
        , debuglevel = debuglevel
    )
    if verbose: print(f'ordered reach computation complete for doing headwaters first')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

    ##STEP 3c -- Opportunistic Network Search
    # This method does as many as it can immediately in each round 
    # (which means that there are not many left in the later rounds 
    # of parallelization...)
    if showtiming: start_time = time.time()
    if verbose: print(f'executing computation on reaches opportunistically')
    if verbose: print(f'(this takes a little longer... we need to improve the method of assigning the opportunistic order...)')

    large_networks = {terminal_segment: network \
                      for terminal_segment, network in networks.items() \
                      if network['maximum_reach_seqorder'] > parallel_split}
    # print(large_networks)
    compute_network_parallel_opportunistic(
        large_networks
        , supernetwork_data = supernetwork_data
        # , connections = connections
        # , verbose = False
        , verbose = verbose
        , debuglevel = debuglevel
    )
    if verbose: print(f'opportunistic reach computation complete')
    if showtiming: print("... in %s seconds." % (time.time() - start_time))

if __name__ == '__main__':
    main()
