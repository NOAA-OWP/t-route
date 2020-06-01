def read_netcdf(geo_file_path):
    import xarray as xr

    return xr.open_dataset(geo_file_path).to_dataframe()


def read_csv(geo_file_path, header='infer', layer_string=None):
    import pandas as pd
    import zipfile

    if geo_file_path.endswith(".zip"):
        if layer_string is None:
            raise ValueError("layer_string is needed if reading from compressed csv")
        with zipfile.ZipFile(geo_file_path, "r") as zcsv:
            with zcsv.open(layer_string) as csv:
                return pd.read_csv(csv, header=header)
    else:
        return pd.read_csv(geo_file_path)


def read_geopandas(geo_file_path, layer_string=None, driver_string=None):
    import geopandas as gpd

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
