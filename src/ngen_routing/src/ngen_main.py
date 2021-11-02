import argparse
import time
from datetime import datetime
from collections import defaultdict
import pathlib
import pandas as pd
import sys
import os

if not hasattr(sys, 'argv'):
    sys.argv  = ['']

def ngen_main(argv):
    from nwm_routing.__main__ import main_v02
    main_v02(argv)


