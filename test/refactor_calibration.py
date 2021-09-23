#!/usr/bin/env python3

import argparse
import os
import geopandas as gpd
import pandas as pd
import json
import numpy as np
import xarray as xr
from pathlib import Path
from nwm_routing.__main__ import _run_everything_v02
from troute.nhd_io import read_custom_input
import time
import glob
import seaborn as sns
import matplotlib.pyplot as plt

huc_id = '01a'
nts = 288
testing_dir = '/media/sf_Inland_hydrofabric/testing'

testing_dir = Path(testing_dir)

inputs_dir = testing_dir / "inputs"
routlink_netcdf_filename = inputs_dir / "v2.1" / 'RouteLink_CONUS.nc'
nwm_v2_1_streams_filename = inputs_dir / "v2.1" / 'ngen-reference' / str("ngen_reference_" + str(huc_id) + ".gpkg")
refactored_streams_dir = inputs_dir / "refactored" / str(huc_id)
refactored_hydrofabrics = os.listdir(refactored_streams_dir)
 
        
outputs_dir = testing_dir / "outputs" / str(huc_id)
   
aggregate_cn_table_filename = outputs_dir / 'aggregate_cn_table.csv'
aggregate_cn_summary_table_filename = outputs_dir / 'aggregate_cn_summary_table.csv'
nwm_v2_1_streams_subset_filename = outputs_dir / str("routlink_" + str(huc_id) + ".gpkg")

diagnostic_dir = outputs_dir / 'diagnostic'
  
aggregate_cn_summary_table = pd.read_csv(aggregate_cn_summary_table_filename)


### Plots

sns.set(style="ticks", rc={"grid.linewidth": 0.1})
sns.set_context("paper", font_scale=2)


# objective function: equal weight -> sum(1-(reference_segments/refactored_segments) - 1 (reference_cn_violations/refactored_cn_violations))
total_reference_segments = aggregate_cn_summary_table.loc[aggregate_cn_summary_table.Network==str('reference_' + str(huc_id))]['Total Segments']
total_reference_violating_segments = aggregate_cn_summary_table.loc[aggregate_cn_summary_table.Network==str('reference_' + str(huc_id))]['Number of Segments Violating CN']



# Initialize a grid of plots with an Axes for each walk
# grid = sns.FacetGrid(aggregate_cn_summary_table, col="Number of Segments Violating CN", hue="Network", palette="tab20c",
#                      col_wrap=4, height=3, sharey=False)
plt.plot(aggregate_cn_summary_table["Network"][1:],aggregate_cn_summary_table["Number of Segments Violating CN"][1:])


plt.set_axis_labels(x_var="Network", y_var="Number of Segments Violating CN")
# Adjust the arrangement of the plots
grid.fig.tight_layout(w_pad=1)