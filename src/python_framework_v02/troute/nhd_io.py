import zipfile
import xarray as xr
import pandas as pd
import geopandas as gpd
import json
import yaml
import numpy as np
from toolz import compose
import dask.array as da
import sys
import math
from datetime import *
import pathlib

from troute.nhd_network import reverse_dict


def read_netcdf(geo_file_path):
    with xr.open_dataset(geo_file_path) as ds:
        return ds.to_dataframe()

def read_json(file_path):
    with open(file_path) as data_file:
        json_data = json.load(data_file)

        def node_key_func(x):
            return int(x[3:])

        already_read_first_row = False

        for key_wb, value_params in json_data.items():
            key_wb = node_key_func(key_wb)

            if not already_read_first_row:
                df_main = pd.json_normalize(value_params)
                df_main['ID'] = key_wb
                already_read_first_row = True

            else:
                df = pd.json_normalize(value_params)
                df['ID'] = key_wb
                df_main = pd.concat([df_main, df], ignore_index=True)

    return df_main


def read_csv(geo_file_path, header="infer", layer_string=None):
    if geo_file_path.suffix == ".zip":
        if layer_string is None:
            raise ValueError("layer_string is needed if reading from compressed csv")
        with zipfile.ZipFile(geo_file_path, "r") as zcsv:
            with zcsv.open(layer_string) as csv:
                return pd.read_csv(csv, header=header)
    else:
        return pd.read_csv(geo_file_path, header=header)


def read_geopandas(geo_file_path, layer_string=None, driver_string=None):
    return gpd.read_file(geo_file_path, driver=driver_string, layer_string=layer_string)


def read(geo_file_path, layer_string=None, driver_string=None):
    if geo_file_path.suffix == ".nc":
        return read_netcdf(geo_file_path)
    elif geo_file_path.suffix == ".json":
        return read_json(geo_file_path)
    else:
        return read_geopandas(
            geo_file_path, layer_string=layer_string, driver_string=driver_string
        )


def read_mask(path, layer_string=None):
    return read_csv(path, header=None, layer_string=layer_string)


def read_custom_input_new(custom_input_file):
    if custom_input_file[-4:] == "yaml":
        with open(custom_input_file) as custom_file:
            data = yaml.load(custom_file, Loader=yaml.SafeLoader)
    else:
        with open(custom_input_file) as custom_file:
            data = json.load(custom_file)
    log_parameters = data.get("log_parameters", {})
    network_topology_parameters = data.get("network_topology_parameters", None)
    supernetwork_parameters = network_topology_parameters.get(
        "supernetwork_parameters", None
    )
    waterbody_parameters = network_topology_parameters.get("waterbody_parameters", None)
    compute_parameters = data.get("compute_parameters", {})
    forcing_parameters = compute_parameters.get("forcing_parameters", {})
    restart_parameters = compute_parameters.get("restart_parameters", {})
    diffusive_parameters = compute_parameters.get("diffusive_parameters", {})
    data_assimilation_parameters = compute_parameters.get(
        "data_assimilation_parameters", {}
    )
    output_parameters = data.get("output_parameters", {})
    parity_parameters = output_parameters.get("wrf_hydro_parity_check", {})

    # TODO: add error trapping for potentially missing files
    return (
        log_parameters,
        supernetwork_parameters,
        waterbody_parameters,
        compute_parameters,
        forcing_parameters,
        restart_parameters,
        diffusive_parameters,
        output_parameters,
        parity_parameters,
        data_assimilation_parameters,
    )


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
    diffusive_parameters = data.get("diffusive_parameters", {})
    coastal_parameters = data.get("coastal_parameters", {})

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
        diffusive_parameters,
        coastal_parameters,
    )


def read_nexus_to_downstream_flowpath_mapping(nexus_to_downstream_flowpath_mapping_file_path):
    with open(nexus_to_downstream_flowpath_mapping_file_path) as json_file:
        nexus_to_downstream_flowpath_mapping_dict = json.load(json_file)

    return nexus_to_downstream_flowpath_mapping_dict


def read_nexus_file(nexus_file_path):

    #Currently reading data in format:
    #[
       #{
         #"ID": "wb-44",
         #"toID": "nex-45"
         #},
    with open(nexus_file_path) as data_file:
        json_data_list = json.load(data_file)

    nexus_to_downstream_flowpath_dict_str = {}

    for id_dict in json_data_list:
        if "nex" in id_dict['ID']:
            nexus_to_downstream_flowpath_dict_str[id_dict['ID']] = id_dict['toID']

    def node_key_func_nexus(x):
        return int(x[4:])

    def node_key_func_wb(x):
        return int(x[3:])

    # Extract the ID integer values
    nexus_to_downstream_flowpath_dict = {node_key_func_nexus(k): node_key_func_wb(v) for k, v in nexus_to_downstream_flowpath_dict_str.items()}

    return nexus_to_downstream_flowpath_dict


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


def read_reservoir_parameter_file(
    reservoir_parameter_file, lake_index_field="lake_id", lake_id_mask=None
):
    """
    Reads reservoir parameter file, which is separate from the LAKEPARM file.
    This function is only called if Hybrid Persistence or RFC type reservoirs
    are active.
    """
    with xr.open_dataset(reservoir_parameter_file) as ds:
        ds = ds.swap_dims({"feature_id": lake_index_field})

        ds_new = ds["reservoir_type"]

        df1 = ds_new.sel({lake_index_field: list(lake_id_mask)}).to_dataframe()

    return df1


def get_ql_from_csv(qlat_input_file, index_col=0):
    """
    qlat_input_file: comma delimted file with header giving timesteps, rows for each segment
    index_col = 0: column/field in the input file with the segment/link id
    """
    ql = pd.read_csv(qlat_input_file, index_col=index_col)
    ql.index = ql.index.astype(int)
    ql = ql.sort_index(axis="index")
    return ql.astype("float32")


def get_nexus_flows_from_csv(nexus_input_file, index_col=0):
    """
    nexus_input_file: comma delimted file for a single nexus flow with timestep rows
    index_col = 0:
    col1: date-time
    col2: flows
    """
    nexus_flows = pd.read_csv(nexus_input_file, index_col=0, header=None)
    return nexus_flows


def read_qlat(path):
    """
    retained for backwards compatibility with early v02 files
    """
    return get_ql_from_csv(path)


# TODO: Generalize this name -- perhaps `read_wrf_hydro_chrt_mf()`
def get_ql_from_wrf_hydro_mf(
    qlat_files,
    index_col="feature_id",
    value_col="q_lateral",
):
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
        try:
            ql = pd.DataFrame(
                ds[value_col].values.T,
                index=ds[index_col].values[0],
                columns=ds.time.values,
                # dtype=float,
            )
        except:
            ql = pd.DataFrame(
                ds[value_col].values.T,
                index=ds[index_col].values,
                columns=ds.time.values,
                # dtype=float,
            )

    return ql


def drop_all_coords(ds):
    return ds.reset_coords(drop=True)


def write_q_to_wrf_hydro(
    flowveldepth,
    chrtout_files,
    output_folder,
    qts_subdivisions,
    new_extension="TRTE"
):
    """
    Write t-route simulated flows to WRF-Hydro CHRTOUT files.

    Arguments:
        flowveldepth (pandas Data Frame): t-route simulated flow, velocity and depth
        chrtout_files (list): chrtout filepaths
        output_folder (pathlib.Path): folder where updated chrtout files will be written
        qts_subdivisions (int): number of t-route timesteps per WRF-hydro timesteps
    """
    
    # count the number of simulated timesteps
    nsteps = len(flowveldepth.loc[:,::3].columns)
    
    # determine how many files to write results out to
    nfiles_to_write = int(np.floor(nsteps / qts_subdivisions))
    
    if nfiles_to_write > 1:
    
        # open all CHRTOUT files as a single xarray dataset
        with xr.open_mfdataset(chrtout_files[:nfiles_to_write], combine="by_coords") as chrtout:

            # !!NOTE: If break_at_waterbodies == True, segment feature_ids coincident with water bodies do
            # not show up in the flowveldepth dataframe. Re-indexing inserts these missing feature_ids and
            # populates columns with NaN values.
            flowveldepth_reindex = flowveldepth.reindex(chrtout.feature_id.values)

            # unpack, subset, and transpose t-route flow data
            qtrt = flowveldepth_reindex.loc[:, ::3].to_numpy().astype("float32")
            qtrt = qtrt[:, qts_subdivisions-1::qts_subdivisions]
            qtrt = np.transpose(qtrt)

            # construct DataArray for t-route flows, dims, coords, and attrs consistent with CHRTOUT
            qtrt_DataArray = xr.DataArray(
                data=da.from_array(qtrt),
                dims=["time", "feature_id"],
                coords=dict(time=chrtout.time.values, feature_id=chrtout.feature_id.values),
                attrs=dict(description="River Flow, t-route", units="m3 s-1",),
            )

            # add t-route DataArray to CHRTOUT dataset
            chrtout["streamflow_troute"] = qtrt_DataArray

            # group by time
            grp_object = chrtout.groupby("time")

        # build a list of datasets, one for each timestep
        dataset_list = []
        for grp, vals in iter(grp_object):
            dataset_list.append(vals)

        # save a new set of chrtout files to disk that contail t-route simulated flow
        chrtout_files_new = [
            output_folder / (s.name + "." + new_extension) for s in chrtout_files
        ]

        # mfdataset solution - can theoretically be parallelised via dask.distributed
        xr.save_mfdataset(dataset_list, paths=chrtout_files_new)

#     # pure serial solution - saving for timing tests against mfdataset
#     for i, dat in enumerate(dataset_list):
#         dat.to_netcdf(chrtout_files_new[i])


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


def read_netcdfs(paths, dim, transform_func=None):
    def process_one_path(path):
        with xr.open_dataset(path) as ds:
            if transform_func is not None:
                ds = transform_func(ds)
            ds.load()
            return ds

    datasets = [process_one_path(p) for p in paths]
    combined = xr.concat(datasets, dim, combine_attrs = "override")
    return combined


def preprocess_time_station_index(xd):
    stationId_da_mask = list(
        map(compose(bytes.isalnum, bytes.strip), xd.stationId.values)
    )
    stationId = list(map(bytes.strip, xd.stationId[stationId_da_mask].values))
    #stationId_int = xd.stationId[stationId_da_mask].values.astype(int)

    unique_times_str = np.unique(xd.time.values).tolist()

    unique_times = np.array(unique_times_str, dtype="str")

    center_time = xd.sliceCenterTimeUTC

    tmask = []
    for t in unique_times_str:
        tmask.append(xd.time == t)

    data_var_dict = {}
    # TODO: make this input parameters
    data_vars = ("discharge", "discharge_quality")

    for v in data_vars:
        vals = []
        for i, t in enumerate(unique_times_str):
            vals.append(np.where(tmask[i],xd[v].values[stationId_da_mask],np.nan))
        combined = np.vstack(vals).T
        data_var_dict[v] = (["stationId","time"], combined)

    return xr.Dataset(
        data_vars=data_var_dict,
        coords={"stationId": stationId, "time": unique_times},
        attrs={"sliceCenterTimeUTC": center_time},
    )


def get_nc_attributes(nc_list, attribute, file_selection=[0, -1]):
    rv = []
    for fs in file_selection:
        try:
            rv.append(get_attribute(nc_list[fs], attribute))
        except:
            rv.append(-1)
    return rv


def get_attribute(nc_file, attribute):
    # TODO: consider naming as get_nc_attribute
    # -- this is really specific to a netcdf file.

    with xr.open_dataset(nc_file) as xd:
        return xd.attrs[attribute]


def build_filtered_gage_df(segment_gage_df, gage_col="gages"):
    """
    segment_gage_df - dataframe indexed by segment with at least
        one column, gage_col, with the gage ids by segment.
    gage_col - the name of the column containing the gages
        to filter.
    """
    # TODO: use this function to filter the inputs
    # coming into the da read routines below.
    gage_list = list(map(bytes.strip, segment_gage_df[gage_col].values))
    gage_mask = list(map(bytes.isalnum, gage_list))
    segment_gage_df = segment_gage_df.loc[gage_mask, [gage_col]]
    segment_gage_df[gage_col] = segment_gage_df[gage_col].map(bytes.strip)
    return segment_gage_df.to_dict()


def build_lastobs_df(
        lastobsfile,
        routelink,
        wrf_lastobs_flag,
        time_shift = 0,
        gage_id = "gages",
        link_id = "link",
        model_discharge_id = "model_discharge",
        obs_discharge_id = "discharge",
        time_idx_id = "timeInd",
        station_id = "stationId",
        station_idx_id = "stationIdInd",
        time_id = "time",
        discharge_nan = -9999.0,
        ref_t_attr_id = "modelTimeAtOutput",
        route_link_idx = "feature_id",
        # last_nudge_id = "last_nudge",
    ):

    standard_columns = {
        "lastobs_discharge": obs_discharge_id,
        "time_since_lastobs": time_id,
        "gages": gage_id,
        "last_model_discharge": model_discharge_id
    }

    """
    Open lastobs file, import the segment keys from the routelink_file
    and extract discharges.
    """
    # TODO: We should already know the link/gage relationship by this point and can require that as an input
    # TODO: ... so we could get rid of the following handful of lines.
    with xr.open_dataset(routelink) as ds1:
        gage_list = list(map(bytes.strip, ds1.gages.values))
        gage_mask = list(map(bytes.isalnum, gage_list))

        gage_da = list(map(bytes.strip, ds1.gages[gage_mask].values))
        # gage_da = ds1.gages[gage_mask].values.astype(int)

        data_var_dict = {}
        data_vars = ("link", "to", "ascendingIndex")
        for v in data_vars:
            data_var_dict[v] = (["gages"], ds1[v].values[gage_mask])
        ds1 = xr.Dataset(data_vars=data_var_dict, coords={"gages": gage_da})
        station_gage_df = ds1.to_dataframe()

    with xr.open_dataset(lastobsfile) as ds:
        model_discharge_last_ts = ds[model_discharge_id][:,-1].to_dataframe()

        # TODO: Determine if the df_discharges extractions can be performed
        # exclusively on the dataset with no
        # transformation to dataframe.
        # ... like how we are doing for the model_discharge_last_ts
        # NOTE: The purpose of the interpolation is to pull the last non-nan
        # value in the time series forward to the end of the series for easy extraction.
        # Caution should be exercised not to compare modeled values from the last
        # time step with the last valid observed values, as these may not
        # correspond to the same times. (See additional comment below...)

        df_discharges = ds[obs_discharge_id].to_dataframe()
        last_ts = df_discharges.index.get_level_values(time_idx_id)[-1]
        df_discharges[df_discharges[obs_discharge_id] != discharge_nan] = df_discharges[
            df_discharges[obs_discharge_id] != discharge_nan
        ].interpolate(method="linear", axis=1)

        discharge_last_ts = df_discharges[
            df_discharges.index.get_level_values(time_idx_id) == last_ts
        ]
        # ref_time = ds.attrs[ref_t_attr_id]
        ref_time = datetime.strptime(ds.attrs[ref_t_attr_id], "%Y-%m-%d_%H:%M:%S")
        # lastobs_times = ds[time_id][:,-1]
        # lastobs_times = ds[time_id].str.decode("utf-8").to_dataframe()
        lastobs_times = pd.to_datetime(
            ds[time_id][:,-1].to_dataframe()[time_id].str.decode("utf-8"),
            format="%Y-%m-%d_%H:%M:%S",
            errors='coerce',
        )
        lastobs_times = (lastobs_times - ref_time).dt.total_seconds()
        lastobs_times = lastobs_times - time_shift

        lastobs_stations = ds[station_id].to_dataframe()
        lastobs_stations[station_id] = lastobs_stations[station_id].map(bytes.strip)

        ## END OF CONTEXT (Remaining items could be outdented...)
        model_discharge_last_ts = model_discharge_last_ts.join(lastobs_stations)
        model_discharge_last_ts = model_discharge_last_ts.join(lastobs_times)
        model_discharge_last_ts = model_discharge_last_ts.join(discharge_last_ts)
        model_discharge_last_ts = model_discharge_last_ts.loc[
            model_discharge_last_ts[model_discharge_id] != discharge_nan
        ]
        model_discharge_last_ts = model_discharge_last_ts.reset_index().set_index(
            station_id
        )
        model_discharge_last_ts = model_discharge_last_ts.drop(
            [station_idx_id, time_idx_id], axis=1
        )

        model_discharge_last_ts[obs_discharge_id] = model_discharge_last_ts[
            obs_discharge_id
        ].to_frame()

        # TODO: Remove any of the following comments that are no longer needed
        # If predict from lastobs file use last obs file results
        # if lastobs_file == "error-based":
        # elif lastobs_file == "obs-based":  # the wrf-hydro default
        # NOTE:  The following would compare potentially mis-matched
        # obs/model pairs, so it is commented until we can figure out
        # a more robust bias-type persistence.
        # For now, we use only obs-type persistence.
        # It would be possible to preserve a 'last_valid_bias' which would
        # presumably correspond to the last_valid_time.
        # # if wrf_lastobs_flag:
        # #     model_discharge_last_ts[last_nudge_id] = (
        # #         model_discharge_last_ts[obs_discharge_id]
        # #         - model_discharge_last_ts[model_discharge_id]
        # #     )

        # final_df = station_gage_df.join(model_discharge_last_ts[obs_discharge_id])  # Preserve all columns
        final_df = station_gage_df.join(model_discharge_last_ts)
        final_df = final_df.reset_index()
        final_df = final_df.set_index(link_id)
        # final_df = final_df.drop([gage_id], axis=1)  # Not needed -- gage id could be useful
        # TODO: What effect does this dropna have?
        final_df = final_df.dropna()
        # Translate to standard column names
        final_df = final_df.rename(columns=reverse_dict(standard_columns))

        # Else predict from the model outputs from t-route if index doesn't match interrupt computation as the results won't be valid
        # else:
        #     fvd_df = fvd_df
        #     if len(model_discharge_last_ts.index) == len(fvd_df.index):
        #         model_discharge_last_ts["last_nudge"] = (
        #             model_discharge_last_ts["discharge"] - fvd_df[fvd_df.columns[0]]
        #         )
        #     else:
        #         print("THE NUDGING FILE IDS DO NOT MATCH THE FLOWVELDEPTH IDS")
        #         sys.exit()
        # # Predictions created with continuously decreasing deltas until near 0 difference
        # a = 120
        # prediction_df = pd.DataFrame(index=model_discharge_last_ts.index)

        # for time in range(0, 720, 5):
        #     weight = math.exp(time / -a)
        #     delta = pd.DataFrame(
        #         model_discharge_last_ts["last_nudge"] / weight)

        #     if time == 0:
        #         prediction_df[str(time)] = model_discharge_last_ts["last_nudge"]
        #         weight_diff = prediction_df[str(time)] - prediction_df[str(time)]
        #     else:
        #         if weight > 0.1:
        #             prediction_df[str(time)] = (
        #                 delta["last_nudge"] + model_discharge_last_ts["model_discharge"]
        #             )
        #         elif weight < -0.1:
        #             prediction_df[str(time)] = (
        #                 delta["last_nudge"] + model_discharge_last_ts["model_discharge"]
        #             )
        # prediction_df["0"] = model_discharge_last_ts["model_discharge"]
        return final_df


def get_usgs_df_from_csv(usgs_csv, routelink_subset_file, index_col="link"):
    """
    routelink_subset_file - provides the gage-->segment crosswalk. Only gages that are represented in the
    crosswalk will be brought into the evaluation.
    usgs_csv - csv file with SEGMENT IDs in the left-most column labeled with "link",
                        and date-headed values from time-slice files in the format
                        "2018-09-18 00:00:00"

    It is assumed that the segment crosswalk and interpolation have both
    already been performed, so we do not need to comprehend
    the potentially non-numeric byte-strings associated with gage IDs, nor
    do we need to interpolate anything here as when we read from the timeslices.

    If that were necessary, we might use a solution such as proposed here:
    https://stackoverflow.com/a/35058538
    note that explicit typing of the index cannot be done on read and
    requires a two-line solution such as:
    ```
    df2 = pd.read_csv(usgs_csv, dtype={index_col:bytes})
    df2 = df2.set_index(index_col)
    ```
    """

    df2 = pd.read_csv(usgs_csv, index_col=index_col)

    with xr.open_dataset(routelink_subset_file) as ds:
        gage_list = list(map(bytes.strip, ds.gages.values))
        gage_mask = list(map(bytes.isdigit, gage_list))

        gage_da = ds[index_col][gage_mask].values.astype(int)

        data_var_dict = {}
        data_vars = ("gages", "to", "ascendingIndex")
        for v in data_vars:
            data_var_dict[v] = ([index_col], ds[v].values[gage_mask])
        ds = xr.Dataset(data_vars=data_var_dict, coords={index_col: gage_da})
        df = ds.to_dataframe()

    usgs_df = df.join(df2)
    usgs_df = usgs_df.drop(["gages", "ascendingIndex", "to"], axis=1)

    return usgs_df


def get_usgs_from_time_slices_folder(
    routelink_subset_file,
    usgs_files,
    qc_threshold,
    max_fill_1min,
    dt,
    t0 = None,
):

    """
    routelink_subset_file - provides the gage-->segment crosswalk.
        Only gages that are represented in the
        crosswalk will be brought into the evaluation.
    usgs_files - list of "time-slice" files containing observed values
    qc_threshold - sets the lowest acceptable quality value;
        lower values will cause the associated obs value to be discarded
        and replaced with NaN.
    max_fill_1min - sets the maximum interpolation length
    t0 - optional date parameter to trim the front of the files -- if not provided,
        the interpolated values are truncated so that the first value returned
        corresponds to the first center date of the first provided file.
    """
    frequency = str(int(dt/60))+"min"
    with read_netcdfs(usgs_files, "time", preprocess_time_station_index,) as ds2:

        # dataframe containing discharge observations
        df2 = pd.DataFrame(
            ds2["discharge"].values,
            index=ds2["stationId"].values,
            columns=ds2.time.values,
        )

        # dataframe containing discharge quality flags [0,1]
        df_qual = pd.DataFrame(
            ds2["discharge_quality"].values/100,
            index=ds2["stationId"].values,
            columns=ds2.time.values,
        )

    with xr.open_dataset(routelink_subset_file) as ds:
        gage_list = list(map(bytes.strip, ds.gages.values))
        gage_mask = list(map(bytes.isalnum, gage_list))

        gage_da = list(map(bytes.strip, ds.gages[gage_mask].values))
        # gage_da = ds.gages[gage_mask].values.astype(int)

        data_var_dict = {}
        data_vars = ("link", "to", "ascendingIndex")
        for v in data_vars:
            data_var_dict[v] = (["gages"], ds[v].values[gage_mask])
        ds = xr.Dataset(data_vars=data_var_dict, coords={"gages": gage_da})
    df = ds.to_dataframe()

    usgs_df = (df.join(df2).
               reset_index().
               rename(columns={"index": "gages"}).
               set_index("link").
               drop(["gages", "ascendingIndex", "to"], axis=1))

    usgs_qual_df = (df.join(df_qual).
               reset_index().
               rename(columns={"index": "gages"}).
               set_index("link").
               drop(["gages", "ascendingIndex", "to"], axis=1))

    # Start and end times of the data obtained from the timeslice dataset
    date_time_strs = usgs_df.columns.tolist()
    date_time_data_start = datetime.strptime(date_time_strs[0], "%Y-%m-%d_%H:%M:%S")
    date_time_data_end = datetime.strptime(date_time_strs[-1], "%Y-%m-%d_%H:%M:%S")

    # Nominal start and end times of the data in timeslice dataset
    # TODO: Consider the case of missing timeslice files...
    # The current method could be fragile in the event of a
    # missing first or last file.
    (first_center_time, last_center_time) = get_nc_attributes(usgs_files, "sliceCenterTimeUTC", (0,-1))
    date_time_center_start = datetime.strptime(first_center_time, "%Y-%m-%d_%H:%M:%S")
    date_time_center_end = datetime.strptime(last_center_time, "%Y-%m-%d_%H:%M:%S")

    dates = []
    for j in pd.date_range(date_time_center_start, date_time_center_end, freq=frequency):
        dates.append(j)
    """
    # dates_to_drop = ~usgs_df.columns.isin(dates)
    OR
    # dates_to_drop = usgs_df.columns.difference(dates)
    # dates_to_add = pd.Index(dates).difference(usgs_df.columns)
    """

    # TODO: Implement a robust test verifying the intended output of the interpolation
    """
    The idea would be to have a function in Cython which could be called in a testable
    framework with inputs such as the following (and corresponding expected outputs for
    the various combinations) for use with pytest or the like:
    ```
    obs = [ 10, 11, 14, 18, 30, 32, 26, 20, 14, 12, 11, 10, 10, 10, 10]
    obs_gap1 = [ None, None, None, 18, 30, 32, 26, 20, 14, 12, 11, 10, 10, 10, 10]
    obs_gap2 = [ 10, None, None, 18, 30, 32, 26, 20, 14, 12, 11, 10, 10, 10, 10]
    obs_gap3 = [ 10, 11, 14, None, None, None, 26, 20, 14, 12, 11, 10, 10, 10, 10]
    obs_gap4 = [ 10, 11, 14, 18, 30, 32, 26, None, None, None, 11, 10, 10, 10, 10]
    obs_gap5 = [ 10, 11, 14, 18, 30, 32, 26, 20, 14, 12, 11, None, None, None, None]
    modeled_low = [ 8, 9, 12, 16, 28, 30, 24, 18, 12, 10, 9, 8, 8, 8, 8]
    modeled_high = [ 12, 13, 16, 20, 32, 34, 28, 22, 16, 14, 13, 12, 12, 12, 12]
    modeled_shift_late = [ 10, 10, 10, 11, 14, 18, 30, 32, 26, 20, 14, 12, 11, 10, 10]
    modeled_shift_late = [ 11, 14, 18, 30, 32, 26, 20, 14, 12, 11, 10, 10, 10, 10, 10]
    lastobs = {"obs":9.5, "time":0}  # Most recent observation at simulation start
    lastobs_old = {"obs":9.5, "time":-3600}  # Most recent observation 1 hour ago
    lastobs_NaN = {"obs":NaN, "time":NaT}  # No valid recent observation
    ```
    """

    #TODO: separate the interpolation into a function; eventually, the data source
    # could be something other than the time-slice files, but the interpolation
    # might be the same and the function would facilitate taking advantage of that.

    # ---- Laugh testing ------
    # screen-out erroneous qc flags
    usgs_qual_df = usgs_qual_df.mask(usgs_qual_df < 0, np.nan)
    usgs_qual_df = usgs_qual_df.mask(usgs_qual_df > 1, np.nan)

    # screen-out poor quality flow observations
    usgs_df = usgs_df.mask(usgs_qual_df < qc_threshold, np.nan)

    # screen-out erroneous flow observations
    usgs_df = usgs_df.mask(usgs_df <= 0, np.nan)

    # ---- Interpolate USGS observations to time discretization of the simulation ----
    usgs_df_T = usgs_df.transpose()
    usgs_df_T.index = pd.to_datetime(usgs_df_T.index, format = "%Y-%m-%d_%H:%M:%S")

    """
    Note: The max_fill is applied when the series is being considered at a 1 minute interval
    so 14 minutes ensures no over-interpolation with 15-minute gage records, but creates
    square-wave signals at gages reporting hourly...
    therefore, we use a 59 minute gap filling tolerance.
    """
    if t0:
        date_time_center_start = t0
    # TODO: Add reporting interval information to the gage preprocessing (timeslice generation)
    usgs_df_T = (usgs_df_T.resample('min').
                 interpolate(limit = max_fill_1min, limit_direction = 'both').
                 resample(frequency).
                 asfreq().
                 loc[date_time_center_start:,:])

    # usgs_df_T.reindex(dates)
    usgs_df_new = usgs_df_T.transpose()

    return usgs_df_new


def get_param_str(target_file, param):
    # TODO: remove this duplicate function
    return get_attribute(target_file, param)


# TODO: Move channel restart above usgs to keep order with execution script
def get_channel_restart_from_csv(
    channel_initial_states_file,
    index_col=0,
    default_us_flow_column="qu0",
    default_ds_flow_column="qd0",
    default_depth_column="h0",
):
    """
    channel_initial_states_file: CSV standard restart file
    index_col = 0: column/field in the input file with the segment/link id
    NOT USED YET default_us_flow_column: name used in remainder of program to refer to this column of the dataset
    NOT USED YET default_ds_flow_column: name used in remainder of program to refer to this column of the dataset
    NOT USED YET default_depth_column: name used in remainder of program to refer to this column of the dataset
    """
    q0 = pd.read_csv(channel_initial_states_file, index_col=index_col)
    q0.index = q0.index.astype(int)
    q0 = q0.sort_index(axis="index")
    return q0.astype("float32")


def get_channel_restart_from_wrf_hydro(
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


def write_channel_restart_to_wrf_hydro(
    data,
    restart_files,
    output_folder,
    channel_initial_states_file,
    dt_troute,
    nts_troute,
    t0,
    crosswalk_file,
    channel_ID_column,
    new_extension,
    restart_file_dimension_var="links",
    troute_us_flow_var_name="qlink1_troute",
    troute_ds_flow_var_name="qlink2_troute",
    troute_depth_var_name="hlink_troute",
):
    """
    Write t-route flow and depth data to WRF-Hydro restart files. New WRF-Hydro restart
    files are created that contain all of the data in the original files, plus t-route
    flow and depth data.

    Agruments
    ---------
        data (Data Frame): t-route simulated flow, velocity and depth data
        restart_files (list): globbed list of WRF-Hydro restart files
        output_folder (pathlib.Path): folder where updated restart files will be written
        channel_initial_states_file (str): WRF-HYDRO standard restart file used to initiate t-route simulation
        dt_troute (int): timestep of t-route simulation (seconds)
        nts_troute (int): number of t-route simulation timesteps
        crosswalk_file (str): File containing reservoir IDs IN THE ORDER of the Restart File
        channel_ID_column (str): field in the crosswalk file to assign as the index of the restart values
        restart_file_dimension_var (str): name of flow and depth data dimension in the Restart File
        troute_us_flow_var_name (str):
        troute_ds_flow_var_name (str):
        troute_depth_var_name (str):

    Returns
    -------
    """
    # create t-route simulation timestamp array
    # TODO: Replace this with a more general function to identify
    # the model timesteps of a particular run.
    # consider as a candidate something like:
    # pd.date_range(start = datetime.strptime('2017-12-31T06:00:00',"%Y-%m-%dT%H:%M:%S"), periods = nts_troute, freq = '300s')
    # or the equivalent resulting from
    # pd.date_range(start = np.datetime64('2017-12-31T06:00:00'), periods = nts_troute, freq = '300s')
    t_array = np.array(t0, dtype=np.datetime64)
    troute_dt = np.timedelta64(dt_troute, "s")
    troute_timestamps = t_array + np.arange(nts_troute) * troute_dt

    # extract ordered feature_ids from crosswalk file
    # TODO: Find out why we re-index this dataset when it
    # already has a segment index.
    with xr.open_dataset(crosswalk_file) as xds:
        xdf = xds[channel_ID_column].to_dataframe()
    xdf = xdf.reset_index()
    xdf = xdf[[channel_ID_column]]

    # reindex flowveldepth array
    flowveldepth_reindex = data.reindex(xdf.link)

    # TODO: The comment "revise do this one at a time" appears to be done... is it?
    # get restart timestamps - revise do this one at a time
    for f in restart_files:
        with xr.open_dataset(f) as ds:

            # get timestamp from restart file
            t = np.array(ds.Restart_Time.replace("_", " "), dtype=np.datetime64)

            # find index troute_timestamp value that matches restart file timestamp
            a = np.where(troute_timestamps == t)[0].tolist()

            # if the restart timestamp exists in the t-route simulatuion
            if len(a) > 0:

                # TODO: consider Packaging the following into a function
                # e.g., def get_restart_vals_from_fvd(a, flowveldepth, idx=0,):
                # pull flow data from flowveldepth array, package into DataArray
                # !! TO DO - is there a more percise way to slice flowveldepth array?
                qtrt = (
                    flowveldepth_reindex.iloc[:, ::3]
                    .iloc[:, a]
                    .to_numpy()
                    .astype("float32")
                )
                qtrt = qtrt.reshape((len(flowveldepth_reindex,)))
                qtrt_DataArray = xr.DataArray(
                    data=qtrt, dims=[restart_file_dimension_var],
                )

                # pull depth data from flowveldepth array, package into DataArray
                # !! TO DO - is there a more percise way to slice flowveldepth array?
                htrt = (
                    flowveldepth_reindex.iloc[:, 2::3]
                    .iloc[:, a]
                    .to_numpy()
                    .astype("float32")
                )
                htrt = htrt.reshape((len(flowveldepth_reindex,)))
                htrt_DataArray = xr.DataArray(
                    data=htrt, dims=[restart_file_dimension_var],
                )

                # insert troute data into restart dataset
                ds[troute_us_flow_var_name] = qtrt_DataArray
                ds[troute_ds_flow_var_name] = qtrt_DataArray
                ds[troute_depth_var_name] = htrt_DataArray

                # write edited to disk with new filename
                ds.to_netcdf(output_folder / (f.name + "." + new_extension))


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


def build_coastal_dataframe(coastal_boundary_elev):
    coastal_df = pd.read_csv(
        coastal_boundary_elev, sep="  ", header=None, engine="python"
    )
    return coastal_df


def build_coastal_ncdf_dataframe(coastal_ncdf):
    with xr.open_dataset(coastal_ncdf) as ds:
        coastal_ncdf_df = ds[["elev", "depth"]]
        return coastal_ncdf_df.to_dataframe()


def lastobs_df_output(
    lastobs_df,
    dt,
    nts,
    t0,
    gages,
    lastobs_output_folder=False,
):

    # join gageIDs to lastobs_df
    lastobs_df = lastobs_df.join(gages)

    # timestamp of last simulation timestep
    modelTimeAtOutput = t0 + timedelta(seconds = nts * dt)
    modelTimeAtOutput_str = modelTimeAtOutput.strftime('%Y-%m-%d_%H:%M:%S')

    # timestamp of last observation
    var = [timedelta(seconds=d) for d in lastobs_df.time_since_lastobs.fillna(0)]
    lastobs_timestamp = [modelTimeAtOutput - d for d in var]
    lastobs_timestamp_str = [d.strftime('%Y-%m-%d_%H:%M:%S') for d in lastobs_timestamp]

    # create xarray Dataset similarly structured to WRF-generated lastobs netcdf files
    ds = xr.Dataset(
        {
            "stationId": (["stationIdInd"], lastobs_df["gages"].to_numpy(dtype = '|S15')),
            "time": (["stationIdInd"], np.asarray(lastobs_timestamp_str,dtype = '|S19')),
            "discharge": (["stationIdInd"], lastobs_df["lastobs_discharge"].to_numpy()),
        }
    )
    ds.attrs["modelTimeAtOutput"] = "example attribute"

    # write-out LastObs file as netcdf
    output_path = pathlib.Path(lastobs_output_folder + "/nudgingLastObs." + modelTimeAtOutput_str + ".nc").resolve()
    ds.to_netcdf(str(output_path))
