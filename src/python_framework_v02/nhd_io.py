import zipfile

import xarray as xr
import pandas as pd
import geopandas as gpd

def read_netcdf(geo_file_path):
    return xr.open_dataset(geo_file_path).to_dataframe()


def read_csv(geo_file_path, header='infer', layer_string=None):
    if geo_file_path.endswith(".zip"):
        if layer_string is None:
            raise ValueError("layer_string is needed if reading from compressed csv")
        with zipfile.ZipFile(geo_file_path, "r") as zcsv:
            with zcsv.open(layer_string) as csv:
                return pd.read_csv(csv, header=header)
    else:
        return pd.read_csv(geo_file_path)


def read_geopandas(geo_file_path, layer_string=None, driver_string=None):
    return gpd.read_file(geo_file_path, driver=driver_string, layer_string=layer_string)


def read(geo_file_path, layer_string=None, driver_string=None):
    if geo_file_path.endswith(".nc"):
        return read_netcdf(geo_file_path)
    else:
        return read_geopandas(
            geo_file_path, layer_string=layer_string, driver_string=driver_string
        )


def read_mask(path, layer_string=None):
    return read_csv(path, header=None, layer_string=layer_string)


def read_qlat(path):
    ql = pd.read_csv(path, index_col=0)
    ql.index = ql.index.astype(int)
    ql.columns = ql.columns.astype(int)
    ql = ql.sort_index(axis='index')
    return ql.astype('float32')


def replace_downstreams(routelink_data, downstream_col, terminal_code):
    ds0_mask = routelink_data[downstream_col] == terminal_code
    new_data = routelink_data.copy()
    new_data.loc[ds0_mask, downstream_col] = ds0_mask.index[ds0_mask]

    # Also set negative any nodes in downstream col not in data.index
    new_data.loc[~routelink_data[downstream_col].isin(routelink_data.index), downstream_col] *= -1
    return new_data


def read_waterbody_df(parm_file, lake_id_mask=None):
    df1 = xr.open_dataset(parm_file).to_dataframe()
    df1 = df1.set_index("lake_id").sort_index()
    if lake_id_mask is None:
        return df1
    return df1.loc[lake_id_mask]
