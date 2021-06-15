#!/usr/bin/env python
from troute.routing.compute import compute_nhd_routing_v02
import troute.nhd_io as nhd_io

from .preprocess import ngen_preprocess

def main():


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
