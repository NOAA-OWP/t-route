#This program reads in and modifies lateral flows from the next generation water model and feeds them
#into mc_reach.compute_network 
import os
import sys
import time
import numpy as np
import argparse
import pathlib
import pandas as pd
from functools import partial
from itertools import chain, islice
import next_gen_io
import nhd_network
import mc_reach

#The following 2 values are currently hard coded for this test domain
nts = 720  # number of timestep = 1140 * 60(model timestep) = 86400 = day
dt_mc = 300.0  # time interval for MC

#Currently tested on the Sugar Creek domain
ngen_network_df = next_gen_io.read_ngen_network_geojson("./sugar_creek_waterbody_data_subset.geojson")

#Create dictionary mapping each connection ID
ngen_network_dict = dict(zip(ngen_network_df.ID, ngen_network_df.toID))

waterbody_connections = {}

#Cycle through dictionary and extract only the ID integer values
for key in ngen_network_dict:
   to_waterbody_value = ngen_network_dict.get(key)
   key = int(key[4 :])
   to_waterbody_value = int(to_waterbody_value[4 :])
   waterbody_connections[key] = to_waterbody_value

#Convert dictionary connections to data frame and make ID column the index
waterbody_df = pd.DataFrame.from_dict(waterbody_connections, orient='index', columns=['to'])

#Sort ID index column
waterbody_df = waterbody_df.sort_index()

waterbody_df = next_gen_io.replace_downstreams(waterbody_df, "to", 0)

connections = nhd_network.extract_connections(waterbody_df, "to")

#At this moment, set to read lateral flows from the current directory.
#Convert catchment lateral flows to format that can be processed by compute_network
next_gen_io.read_catchment_lateral_flows(".")

#TODO: Figure out whether the outflows should be from nexus or waterbodies
qlats = next_gen_io.read_qlat("./sugar_creek_qlats.csv")

rconn = nhd_network.reverse_network(connections)

subnets = nhd_network.reachable_network(rconn, check_disjoint=False)

waterbody_df['dt'] = 300.0

#Setting all below to 1.0 until we can get the appropriate parameters
waterbody_df['bw'] = 1.0
waterbody_df['tw'] = 1.0
waterbody_df['twcc'] = 1.0
waterbody_df['dx'] = 1.0
waterbody_df['n'] = 1.0
waterbody_df['ncc'] = 1.0
waterbody_df['cs'] = 1.0
waterbody_df['s0'] = 1.0

#Set types as float32
waterbody_df = waterbody_df.astype({"dt": "float32", "bw": "float32", "tw": "float32", "twcc": "float32", "dx": "float32", "n": "float32", "ncc": "float32", "cs": "float32", "s0": "float32"})

subreaches = {}

for tw, net in subnets.items():
    path_func = partial(nhd_network.split_at_junction, net)
    subreaches[tw] = nhd_network.dfs_decomposition(net, path_func)


results = []
for twi, (tw, reach) in enumerate(subreaches.items(), 1):
    r = list(chain.from_iterable(reach))
    data_sub = waterbody_df.loc[r, ['dt', 'bw', 'tw', 'twcc', 'dx', 'n', 'ncc', 'cs', 's0']].sort_index()
    qlat_sub = qlats.loc[r].sort_index()
    results.append(mc_reach.compute_network(
        nts, reach, subnets[tw], data_sub.index.values, data_sub.columns.values, data_sub.values, qlat_sub.values
        )
    )

fdv_columns = pd.MultiIndex.from_product([range(nts), ['q', 'v', 'd']]).to_flat_index()
flowveldepth = pd.concat([pd.DataFrame(d, index=i, columns=fdv_columns) for i, d in results], copy=False)
flowveldepth = flowveldepth.sort_index()
flowveldepth.to_csv(f"Sugar_Creek_MC_Results.csv")
print(flowveldepth)
