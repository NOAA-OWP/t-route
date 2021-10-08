import argparse
import time
from datetime import datetime
from collections import defaultdict
import pathlib
import pandas as pd

## network and reach utilities
#import troute.nhd_network as nhd_network
#import troute.nhd_io as nhd_io
#import troute.nhd_network_utilities_v02 as nnu
#import build_tests  # TODO: Determine whether and how to incorporate this into setup.py
import sys
import os

#https://github.com/googleapis/oauth2client/issues/642
if not hasattr(sys, 'argv'):
    sys.argv  = ['']

    #import troute.routing.diffusive_utils as diff_utils

    #from .input import _input_handler_v02, _input_handler_v03
    #from .preprocess import (
    #    nwm_network_preprocess,
    #    nwm_initial_warmstate_preprocess,
    #    nwm_forcing_preprocess,
    #)
    #from .output import nwm_output_generator

    #from troute.routing.compute import compute_nhd_routing_v02


def set_paths(root_path):
    '''
    framework_path = "../../../python_framework_v02"
    routing_path = "../../../python_routing_v02"
    network_utilities_path = "../../../python_framework_v02/troute"
    #nwm_routing_path =  "."


    framework_path_full = os.path.join(root_path, framework_path)
    routing_path_full = os.path.join(root_path, routing_path)
    network_utilities_path_full = os.path.join(root_path, network_utilities_path)
 
    sys.path.append(framework_path_full)
    sys.path.append(routing_path_full)
    sys.path.append(network_utilities_path_full)

    import troute.nhd_network as nhd_network
    #import compute_nhd_routing_v02
    from troute.routing.compute import compute_nhd_routing_v02
    #import mc_reach

    import troute.nhd_io as nhd_io

    import troute.nhd_network_utilities_v02 as nnu

    #from troute.routing.compute import compute_nhd_routing_v02

    #import nwm_routing.input as input1
 
    from nwm_routing.input import _input_handler_v02, _input_handler_v03
    from nwm_routing.preprocess import (
        nwm_network_preprocess,
        nwm_initial_warmstate_preprocess,
        nwm_forcing_preprocess,
    )
    from nwm_routing.output import nwm_output_generator

    #from nwm_routing.__main__ import main_v02
    #from nwm_routing.a__main__ import main_v02

    #import input._input_handler_v02 as _input_handler_v02
    #import nwm_routing.input._input_handler_v02 as _input_handler_v02

    #from input import _input_handler_v02, _input_handler_v03
    '''


    '''
    from nwm_routing.preprocess import (
        nwm_network_preprocess,
        nwm_initial_warmstate_preprocess,
        nwm_forcing_preprocess,
    )
    from nwm_routing.output import nwm_output_generator
    '''


    #from preprocess import ngen_preprocess
    #import next_gen_io

    return



def ngen_main(argv):
    from nwm_routing.__main__ import main_v02
    #from nwm_routing.a__main__ import main_v02
    main_v02(argv)
    print ("test2----------")


