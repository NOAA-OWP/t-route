import argparse
import time
from datetime import datetime
from collections import defaultdict
import pathlib
import pandas as pd
import sys
import os

from functools import partial
import hashlib
non_crypto_md5 = partial(hashlib.md5, usedforsecurity=False)
hashlib.md5 = non_crypto_md5

if not hasattr(sys, 'argv'):
    sys.argv  = ['']

def ngen_main(argv):
    from nwm_routing.__main__ import main_v02
    main_v02(argv)

if __name__ == "__main__":
    ngen_main(sys.argv[1:])
