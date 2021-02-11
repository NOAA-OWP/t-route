import zipfile

import xarray as xr
import pandas as pd
import geopandas as gpd
import json
import yaml
import numpy as np
import glob
from toolz import compose


def read_netcdf(geo_file_path):
    with xr.open_dataset(geo_file_path) as ds:
        return ds.to_dataframe()


def read_csv(geo_file_path, header="infer", layer_string=None):
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


def read_custom_input(custom_input_file):
    if custom_input_file[-4:] == "yaml":
        with open(custom_input_file) as custom_file:
            data = yaml.load(custom_file, Loader=yaml.SafeLoader)
    else:
        with open(custom_input_file) as custom_file:
            data = json.load(custom_file)
    supernetwork_parameters = data.get("supernetwork_parameters", None)
    waterbody_parameters = data.get("waterbody_parameters", {})
    forcing_parameters = data.get("forcing_parameters", {})
    restart_parameters = data.get("restart_parameters", {})
    output_parameters = data.get("output_parameters", {})
    run_parameters = data.get("run_parameters", {})
    parity_parameters = data.get("parity_parameters", {})
    data_assimilation_parameters = data.get("data_assimilation_parameters", {})
    # TODO: add error trapping for potentially missing files
    return (
        supernetwork_parameters,
        waterbody_parameters,
        forcing_parameters,
        restart_parameters,
        output_parameters,
        run_parameters,
        parity_parameters,
        data_assimilation_parameters,
    )


def replace_downstreams(data, downstream_col, terminal_code):
    ds0_mask = data[downstream_col] == terminal_code
    new_data = data.copy()
    new_data.loc[ds0_mask, downstream_col] = ds0_mask.index[ds0_mask]

    # Also set negative any nodes in downstream col not in data.index
    new_data.loc[~data[downstream_col].isin(data.index), downstream_col] *= -1
    return new_data


def read_waterbody_df(waterbody_parameters, waterbodies_values, wbtype="level_pool"):
    """
    General waterbody dataframe reader. At present, only level-pool
    capability exists.
    """
    if wbtype == "level_pool":
        wb_params = waterbody_parameters[wbtype]
        return read_level_pool_waterbody_df(
            wb_params["level_pool_waterbody_parameter_file_path"],
            wb_params["level_pool_waterbody_id"],
            waterbodies_values[wbtype],
        )


def read_level_pool_waterbody_df(
    parm_file, lake_index_field="lake_id", lake_id_mask=None
):
    """
    Reads LAKEPARM file and prepares a dataframe, filtered
    to the relevant reservoirs, to provide the parameters
    for level-pool reservoir computation.

    Completely replaces the read_waterbody_df function from prior versions
    of the v02 routing code.

    Prior version filtered the dataframe as opposed to the dataset as in this version.
    with xr.open_dataset(parm_file) as ds:
        df1 = ds.to_dataframe()
    df1 = df1.set_index(lake_index_field).sort_index(axis="index")
    if lake_id_mask is None:
        return df1
    else:
        return df1.loc[lake_id_mask]
    """

    # TODO: avoid or parameterize "feature_id" or ... return to name-blind dataframe version
    with xr.open_dataset(parm_file) as ds:
        ds = ds.swap_dims({"feature_id": lake_index_field})
        df1 = ds.sel({lake_index_field: list(lake_id_mask)}).to_dataframe()
    df1 = df1.sort_index(axis="index")
    return df1


def get_ql_from_csv(qlat_input_file, index_col=0):
    """
    qlat_input_file: comma delimted file with header giving timesteps, rows for each segment
    index_col = 0: column/field in the input file with the segment/link id
    """
    ql = pd.read_csv(qlat_input_file, index_col=index_col)
    ql.index = ql.index.astype(int)
    ql.columns = ql.columns.astype(int)
    ql = ql.sort_index(axis="index")
    return ql.astype("float32")


def read_qlat(path):
    """
    retained for backwards compatibility with early v02 files
    """
    return get_ql_from_csv(path)


def get_ql_from_wrf_hydro_mf(qlat_files, index_col="feature_id", value_col="q_lateral"):
    """
    qlat_files: globbed list of CHRTOUT files containing desired lateral inflows
    index_col: column/field in the CHRTOUT files with the segment/link id
    value_col: column/field in the CHRTOUT files with the lateral inflow value

    In general the CHRTOUT files contain one value per time step. At present, there is
    no capability for handling non-uniform timesteps in the qlaterals.

    The qlateral may also be input using comma delimited file -- see
    `get_ql_from_csv`


    Note/Todo:
    For later needs, filtering for specific features or times may
    be accomplished with one of:
        ds.loc[{selectors}]
        ds.sel({selectors})
        ds.isel({selectors})

    Returns from these selection functions are sub-datasets.

    For example:
    ```
    (Pdb) ds.sel({"feature_id":[4186117, 4186169],"time":ds.time.values[:2]})['q_lateral'].to_dataframe()
                                     latitude  longitude  q_lateral
    time                feature_id
    2018-01-01 13:00:00 4186117     41.233807 -75.413895   0.006496
    2018-01-02 00:00:00 4186117     41.233807 -75.413895   0.006460
    ```

    or...
    ```
    (Pdb) ds.sel({"feature_id":[4186117, 4186169],"time":[np.datetime64('2018-01-01T13:00:00')]})['q_lateral'].to_dataframe()
                                     latitude  longitude  q_lateral
    time                feature_id
    2018-01-01 13:00:00 4186117     41.233807 -75.413895   0.006496
    ```
    """
    filter_list = None

    with xr.open_mfdataset(
        qlat_files,
        combine="by_coords",
        # combine="nested",
        # concat_dim="time",
        # data_vars="minimal",
        # coords="minimal",
        # compat="override",
        preprocess=drop_all_coords,
        # parallel=True,
    ) as ds:
        ql = pd.DataFrame(
            ds[value_col].values.T,
            index=ds[index_col].values[0],
            columns=ds.time.values,
            # dtype=float,
        )

    return ql


def drop_all_coords(ds):
    return ds.reset_coords(drop=True)


def get_ql_from_wrf_hydro(qlat_files, index_col="station_id", value_col="q_lateral"):
    """
    qlat_files: globbed list of CHRTOUT files containing desired lateral inflows
    index_col: column/field in the CHRTOUT files with the segment/link id
    value_col: column/field in the CHRTOUT files with the lateral inflow value

    In general the CHRTOUT files contain one value per time step. At present, there is
    no capability for handling non-uniform timesteps in the qlaterals.

    The qlateral may also be input using comma delimited file -- see
    `get_ql_from_csv`
    """

    li = []

    for filename in qlat_files:
        with xr.open_dataset(filename) as ds:
            df1 = ds[["time", value_col]].to_dataframe()

        li.append(df1)

    frame = pd.concat(li, axis=0, ignore_index=False)
    mod = frame.reset_index()
    ql = mod.pivot(index=index_col, columns="time", values=value_col)

    return ql


def read_netcdfs(files, dim, transform_func=None):
    def process_one_path(path):
        with xr.open_dataset(path) as ds:
            if transform_func is not None:
                ds = transform_func(ds)
            ds.load()
            return ds

    paths = sorted(glob.glob(files))
    datasets = [process_one_path(p) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined


def preprocess_time_station_index(xd):
    stationId_da_mask = list(
        map(compose(bytes.isdigit, bytes.strip), xd.stationId.values)
    )
    stationId = xd.stationId[stationId_da_mask].values.astype(int)

    unique_times_str = np.unique(xd.time.values).tolist()

    unique_times = np.array(unique_times_str, dtype="str")

    data_var_dict = {}
    data_vars = ("discharge", "discharge_quality")
    for v in data_vars:
        data_var_dict[v] = (["stationId"], xd[v].values[stationId_da_mask])
    return xr.Dataset(
        data_vars=data_var_dict, coords={"stationId": stationId, "time": unique_times}
    )


def get_usgs_from_time_slices(
    routelink_subset_file, usgs_timeslices_folder, data_assimilation_filter,usgs_csv
):
    if usgs_csv:
        df2 = pd.read_csv(usgs_csv,index_col=0)
    else:
        usgs_files = glob.glob(usgs_timeslices_folder + data_assimilation_filter)

        with read_netcdfs(
            usgs_timeslices_folder + data_assimilation_filter,
            "time",
            preprocess_time_station_index,
        ) as ds2:
            df2 = pd.DataFrame(
                ds2["discharge"].values.T,
                index=ds2["stationId"].values,
                columns=ds2.time.values,
            )


    with xr.open_dataset(routelink_subset_file) as ds:
        gage_list = list(map(bytes.strip, ds.gages.values))
        gage_mask = list(map(bytes.isdigit, gage_list))

        gage_da = ds.gages[gage_mask].values.astype(int)

        data_var_dict = {}
        data_vars = ("link", "to", "ascendingIndex")
        for v in data_vars:
            data_var_dict[v] = (["gages"], ds[v].values[gage_mask])
        ds = xr.Dataset(data_vars=data_var_dict, coords={"gages": gage_da})
    df = ds.to_dataframe()

    usgs_df = df.join(df2)
    usgs_df = usgs_df.reset_index()
    usgs_df = usgs_df.rename(columns={"index": "gages"})
    usgs_df = usgs_df.set_index("link")
    usgs_df = usgs_df.drop(["gages", "ascendingIndex", "to"], axis=1)
    columns_list = usgs_df.columns
    
    for i in range(0, (len(columns_list) * 3) - 12, 12):
        original_string = usgs_df.columns[i]
        original_string_shortened = original_string[:-5]
        temp_name1 = original_string_shortened + str("05:00")
        temp_name2 = original_string_shortened + str("10:00")
        temp_name3 = original_string_shortened + str("20:00")
        temp_name4 = original_string_shortened + str("25:00")
        temp_name5 = original_string_shortened + str("35:00")
        temp_name6 = original_string_shortened + str("40:00")
        temp_name7 = original_string_shortened + str("50:00")
        temp_name8 = original_string_shortened + str("55:00")
        usgs_df.insert(i + 1, temp_name1, np.nan)
        usgs_df.insert(i + 2, temp_name2, np.nan)
        usgs_df.insert(i + 4, temp_name3, np.nan)
        usgs_df.insert(i + 5, temp_name4, np.nan)
        usgs_df.insert(i + 7, temp_name5, np.nan)
        usgs_df.insert(i + 8, temp_name6, np.nan)
        usgs_df.insert(i + 10, temp_name7, np.nan)
        usgs_df.insert(i + 11, temp_name8, np.nan)

    usgs_df = usgs_df.interpolate(method="linear", axis=1)
    usgs_df.drop(usgs_df[usgs_df.iloc[:,0] == -999999.000000].index , inplace=True)

    return usgs_df


def get_stream_restart_from_wrf_hydro(
    channel_initial_states_file,
    crosswalk_file,
    channel_ID_column,
    us_flow_column="qlink1",
    ds_flow_column="qlink2",
    depth_column="hlink",
    default_us_flow_column="qu0",
    default_ds_flow_column="qd0",
    default_depth_column="h0",
):
    """
    channel_initial_states_file: WRF-HYDRO standard restart file
    crosswalk_file: File containing channel IDs IN THE ORDER of the Restart File
    channel_ID_column: field in the crosswalk file to assign as the index of the restart values
    us_flow_column: column in the restart file to use for upstream flow initial state
    ds_flow_column: column in the restart file to use for downstream flow initial state
    depth_column: column in the restart file to use for depth initial state
    default_us_flow_column: name used in remainder of program to refer to this column of the dataset
    default_ds_flow_column: name used in remainder of program to refer to this column of the dataset
    default_depth_column: name used in remainder of program to refer to this column of the dataset

    The Restart file gives hlink, qlink1, and qlink2 values for channels --
    the order is simply the same as that found in the Route-Link files.
    *Subnote 1*: The order of these values is NOT the order found in the CHRTOUT files,
    though the number of outputs is the same as in those files.
    *Subnote 2*: If there is no mask applied, then providing the path to the RouteLink file should
    be sufficient for the crosswalk file. If there is a mask applied (so that
    not all segments in the RouteLink file are used in the routing computations) then
    a pre-processing step will be need to provide only the relevant segments in the
    crosswalk file.
    """

    with xr.open_dataset(crosswalk_file) as xds:
        xdf = xds[channel_ID_column].to_dataframe()
    xdf = xdf.reset_index()
    xdf = xdf[[channel_ID_column]]
    with xr.open_dataset(channel_initial_states_file) as qds:
        if depth_column in qds:
            qdf2 = qds[[us_flow_column, ds_flow_column, depth_column]].to_dataframe()
        else:
            qdf2 = qds[[us_flow_column, ds_flow_column]].to_dataframe()
            qdf2[depth_column] = 0
    qdf2 = qdf2.reset_index()
    qdf2 = qdf2[[us_flow_column, ds_flow_column, depth_column]]
    qdf2.rename(
        columns={
            us_flow_column: default_us_flow_column,
            ds_flow_column: default_ds_flow_column,
            depth_column: default_depth_column,
        },
        inplace=True,
    )
    qdf2[channel_ID_column] = xdf
    qdf2 = qdf2.reset_index().set_index([channel_ID_column])

    q_initial_states = qdf2

    q_initial_states = q_initial_states.drop(columns="index")

    return q_initial_states


def get_reservoir_restart_from_wrf_hydro(
    waterbody_intial_states_file,
    crosswalk_file,
    waterbody_ID_field,
    crosswalk_filter_file=None,
    crosswalk_filter_file_field=None,
    waterbody_flow_column="qlakeo",
    waterbody_depth_column="resht",
    default_waterbody_flow_column="qd0",
    default_waterbody_depth_column="h0",
):
    """
    waterbody_intial_states_file: WRF-HYDRO standard restart file
    crosswalk_file: File containing reservoir IDs IN THE ORDER of the Restart File
    waterbody_ID_field: field in the crosswalk file to assign as the index of the restart values
    crosswalk_filter_file: file containing only reservoir ids needed in the crosswalk file (see notes below)
    crosswalk_filter_file_field: waterbody id field in crosswalk filter file
    waterbody_flow_column: column in the restart file to use for upstream flow initial state
    waterbody_depth_column: column in the restart file to use for downstream flow initial state
    default_waterbody_flow_column: name used in remainder of program to refer to this column of the dataset
    default_waterbody_depth_column: name used in remainder of program to refer to this column of the dataset

    The Restart file gives qlakeo and resht values for waterbodies.
    The order of values in the file is the same as the order in the LAKEPARM file from WRF-Hydro.
    However, there are many instances where only a subset of waterbodies described in the lakeparm
    file are used. In these cases, a filter file must be provided which specifies
    which of the reservoirs in the crosswalk file are to be used.
    """

    with xr.open_dataset(crosswalk_file) as xds:
        X = xds[waterbody_ID_field]

        if crosswalk_filter_file:
            with xr.open_dataset(crosswalk_filter_file) as fds:
                xdf = X.loc[X.isin(fds[crosswalk_filter_file_field])].to_dataframe()
        else:
            xdf = X.to_dataframe()

    xdf = xdf.reset_index()[waterbody_ID_field]

    # read initial states from r&r output
    with xr.open_dataset(waterbody_intial_states_file) as resds:
        resdf = resds[[waterbody_flow_column, waterbody_depth_column]].to_dataframe()
    resdf = resdf.reset_index()
    resdf = resdf[[waterbody_flow_column, waterbody_depth_column]]
    resdf.rename(
        columns={
            waterbody_flow_column: default_waterbody_flow_column,
            waterbody_depth_column: default_waterbody_depth_column,
        },
        inplace=True,
    )

    mod = resdf.join(xdf).reset_index().set_index(waterbody_ID_field)
    init_waterbody_states = mod

    return init_waterbody_states
