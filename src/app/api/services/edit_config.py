import yaml
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List
import shutil
import argparse

from app.core.settings import Settings

def edit_yaml(original_file: Path, params: Dict[str, str]):
    tmp_yaml = original_file.with_name(original_file.stem + '_tmp_' + params["lid"] + original_file.suffix)
    with open(original_file, 'r') as file:
        data = yaml.safe_load(file)

    output_dir = Path(data["output_parameters"]["stream_output"]["stream_output_directory"].format(params["lid"]))
    output_dir.mkdir(exist_ok=True)
    
    data["network_topology_parameters"]["supernetwork_parameters"]["geo_file_path"] = params["geo_file_path"]
    data["compute_parameters"]["restart_parameters"]["start_datetime"] = params["start_datetime"]
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
        "start_datetime": start_datetime,
        "geo_file_path": geo_file_path,
        "nts": nts,
        "qlat_input_folder": qlat_input_folder,
    }