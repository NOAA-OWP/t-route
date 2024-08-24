import yaml
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List
import shutil
import argparse
import geopandas as gpd

import pandas as pd
import numpy as np

from app.core.settings import Settings

def edit_yaml(original_file: Path, params: Dict[str, str], restart_file: Path):
    tmp_yaml = original_file.with_name(original_file.stem + '_tmp_' + params["lid"] + original_file.suffix)
    with open(original_file, 'r') as file:
        data = yaml.safe_load(file)

    output_dir = Path(data["output_parameters"]["stream_output"]["stream_output_directory"].format(params["lid"]))
    output_dir.mkdir(exist_ok=True)
    
    data["network_topology_parameters"]["supernetwork_parameters"]["geo_file_path"] = params["geo_file_path"]
    
    data["compute_parameters"]["restart_parameters"]["start_datetime"] = params["start_datetime"]
    data["compute_parameters"]["restart_parameters"]["lite_channel_restart_file"] = restart_file.__str__()
    data["compute_parameters"]["forcing_parameters"]["nts"] = params["nts"]
    data["compute_parameters"]["forcing_parameters"]["qlat_input_folder"] = params["qlat_input_folder"]
    
    data["output_parameters"]["stream_output"]["stream_output_directory"] = output_dir.__str__()
    
    # Write the edited data to the new file
    with open(tmp_yaml, 'w') as file:
        yaml.dump(data, file)
    
    return tmp_yaml


def create_params(
    lid: str,
    feature_id: str, 
    hy_id: str,
    initial_start: float,
    start_time: str,
    num_forecast_days: int,
    settings: Settings
) -> Dict[str, str]:
    dt = datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%S")
    start_datetime = dt.strftime("%Y-%m-%d_%H:%M")

    nts = 288 * num_forecast_days

    geo_file_path = settings.geofile_path.format(feature_id)
    qlat_input_folder = settings.qlat_input_path.format(lid)
    return {
        "lid": lid,
        "hy_id": hy_id,
        "initial_start": initial_start,
        "start_datetime": start_datetime,
        "geo_file_path": geo_file_path,
        "nts": nts,
        "qlat_input_folder": qlat_input_folder,
    }

def create_initial_start_file(params: Dict[str, str], settings: Settings) -> Path:
    start_datetime = datetime.strptime(params["start_datetime"], "%Y-%m-%d_%H:%M")
    formatted_datetime = start_datetime.strftime("%Y-%m-%d_%H:%M")

    # Pulling the keys out of the gpkg file
    gdf = gpd.read_file(params["geo_file_path"], layer="network")
    mask = gdf["divide_id"].isna()
    keys = [int(val[4:]) for val in set(gdf[~mask]["divide_id"].values.tolist())]

    discharge_upstream = np.zeros([len(keys)])
    discharge_downstream = np.zeros([len(keys)])
    height = np.zeros([len(keys)])
    idx = keys.index(int(params["hy_id"]))
    discharge_upstream[idx] = float(params["initial_start"])

    time_array = np.array([pd.to_datetime(formatted_datetime, format="%Y-%m-%d_%H:%M")] * len(keys))

    df = pd.DataFrame({
        "time": time_array,
        "key": np.array(keys),
        "qu0": discharge_upstream,
        "qd0": discharge_downstream,
        "h0": height, # Todo look into adding stage here

    })
    df.set_index('key', inplace=True)
    restart_path = Path(settings.restart_path.format(params["lid"]))
    restart_path.mkdir(exist_ok=True)
    restart_full_path = restart_path / settings.restart_file.format(formatted_datetime)
    df.to_pickle(restart_full_path)
    return restart_full_path