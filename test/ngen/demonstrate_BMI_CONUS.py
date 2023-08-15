import sys
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

sys.path.append("../t-route/src/")
import bmi_troute

#--------- Load geopackage data tables

file_path = '/home/dongha.kim/github/data/bmi_conus_hyfeature/conus.gpkg'
flowpaths = gpd.read_file(file_path, layer='flowpaths')
flowpath_attributes = gpd.read_file(file_path, layer='flowpath_attributes')
lakes = gpd.read_file(file_path, layer='lakes')
network = gpd.read_file(file_path, layer='network')

import pdb; pdb.set_trace()
# Initialize model with configuration file
troute = bmi_troute.bmi_troute()
troute.initialize(bmi_cfg_file='/home/sean.horvath/projects/t-route/test/ngen/test_conus_AnA.yaml')