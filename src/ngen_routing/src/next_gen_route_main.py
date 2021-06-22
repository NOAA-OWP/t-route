#!/usr/bin/env python
import os
import sys
import time
import numpy as np
import argparse
import pathlib
import pandas as pd
from functools import partial
from itertools import chain, islice
#import next_gen_io
import json

def set_paths(root_path):

  framework_path = "../../python_framework_v02"
  routing_path = "../../python_routing_v02"

  framework_path_full = os.path.join(root_path, framework_path)
  routing_path_full = os.path.join(root_path, routing_path)

  sys.path.append(framework_path_full)
  sys.path.append(routing_path_full)

  import troute.nhd_network as nhd_network
  #import compute_nhd_routing_v02
  #from troute.routing.compute import compute_nhd_routing_v02
  #import mc_reach

  import troute.nhd_io as nhd_io

  from preprocess import ngen_preprocess
  import next_gen_io



  return

def ngen_routing(number_of_timesteps, delta_time):


  return

def call_read_catchment_lateral_flows(path):
  import next_gen_io as next_gen_io

  nexus_flows = next_gen_io.read_catchment_lateral_flows(path)

  print (nexus_flows)

def main():
    from troute.routing.compute import compute_nhd_routing_v02
    import troute.nhd_io as nhd_io

    from .preprocess import ngen_preprocess



    """
    *****************************
        Params that should come from dynamic input
    *****************************
    """
    nts = 720
    dt = 300.0
    
    """
    (
        connections,
        param_df,
        wbodies,
        waterbodies_df,
        break_network_at_waterbodies,
        independent_networks,
        reaches_bytw,
        rconn,
    ) =
    """ 
    ngen_preprocess(
        #supernetwork_parameters,
        #waterbody_parameters,
        #showtiming=showtiming,
        #verbose=verbose,
        #debuglevel=debuglevel,
    )
    """
    results = compute_nhd_routing_v02(
        connections,
        rconn,
        wbodies,
        reaches_bytw,
        compute_func,
        run_parameters.get("parallel_compute_method", None),
        run_parameters.get("subnetwork_target_size", 1),
        # The default here might be the whole network or some percentage...
        run_parameters.get("cpu_pool", None),
        run_parameters.get("dt"),
        run_parameters.get("nts", 1),
        run_parameters.get("qts_subdivisions", 1),
        independent_networks,
        param_df,
        q0,
        qlats,
        usgs_df,
        last_obs_df,
        run_parameters.get("assume_short_ts", False),
        run_parameters.get("return_courant", False),
        waterbodies_df_reduced,
        diffusive_parameters,
    )
    """

if __name__ == "__main__":
    main()
