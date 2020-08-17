#Utility I/O functions for the next generation water model network.
#The read_qlat and replace_downstreams functions are copies from nhd_io.py
#and are placed here for now, so that this program can run without importing that module.
import xarray as xr
import pandas as pd
import geopandas as gpd
import os


def read_ngen_network_csv(path):
    return pd.read_csv(path, index_col=0)


def read_ngen_network_geojson(path):
    return gpd.read_file(path)


def read_ngen_network_json(path):
    return pd.read_json(path)


#Convert catchment lateral flows to format that can be processed by compute_network
def read_catchment_lateral_flows(path):
    ql = []
    qlats_df = pd.DataFrame({'A' : []})

    past_first_loop = False

    catchment_id_list = []

    for file_name in os.listdir(path):
       if file_name.startswith("cat-"):
          file_name_str_list = file_name.split("_")
          catchment_id = int(file_name_str_list[0][4 :])
          catchment_id_list.append(catchment_id)
          catchment_qlats = pd.read_csv(file_name, names=["", catchment_id])
          catchment_qlats = catchment_qlats.iloc[:, [1]]
          catchment_qlats_transposed = catchment_qlats.transpose()
          ql.append(catchment_qlats_transposed)

    qlats_df = pd.concat(ql)
    qlats_df.to_csv('sugar_creek_qlats.csv')


def read_qlat(path):
    ql = pd.read_csv(path, index_col=0)
    ql.index = ql.index.astype(int)
    ql.columns = ql.columns.astype(int)
    ql = ql.sort_index(axis='index')
    return ql.astype('float32')


def replace_downstreams(data, downstream_col, terminal_code):
    ds0_mask = data[downstream_col] == terminal_code
    new_data = data.copy()
    new_data.loc[ds0_mask, downstream_col] = ds0_mask.index[ds0_mask]

    # Also set negative any nodes in downstream col not in data.index
    new_data.loc[~data[downstream_col].isin(data.index), downstream_col] *= -1
    return new_data

