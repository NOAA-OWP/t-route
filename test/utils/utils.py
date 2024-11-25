import argparse
import os
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

def get_fid_data(output_dir: Path, fid_list: List[str]) -> Tuple[np.ndarray, np.ndarray]:
    """Reads a folder and gets the flow data for a lid

    Parameters
    ----------
    output_dir: Path
        the path to the output directory
    
    Returns
    -------

    """
    list_ds = list(output_dir.glob("*.nc"))
    sorted_ds = sorted(list_ds, key=lambda x: os.path.getctime(x))

    flow_array = np.zeros([len(sorted_ds), len(fid_list)])
    time_delta = []
    for idx, file in enumerate(sorted_ds):
        ds = xr.open_dataset(file, engine="netcdf4") 
        try:
            flow_array[idx] = ds.sel(feature_id=fid_list).flow.values.squeeze()
        except KeyError as e:
            raise e
        file_time = file.stem.split("_")[-1]
        time_delta.append(file_time)
    
    datetime_series = pd.to_datetime(time_delta, format='%Y%m%d%H%M').to_numpy()

    return flow_array, datetime_series

def get_rfc_data(input_dir: Path, fid_list: str) -> Tuple[np.ndarray, np.ndarray]:
    list_ds = list(input_dir.glob("*.csv"))
    sorted_ds = sorted(list_ds, key=lambda x: os.path.getctime(x))

    flow_array = []
    time_delta = []
    for idx, file in enumerate(sorted_ds):
        df = pd.read_csv(file) 
        columns = df.columns
        try:
            flow_array.append(df[df["feature_id"].isin(fid_list)][columns[1]].values.tolist())
        except KeyError as e:
            raise e
        file_time = file.stem.split(".")[0]
        time_delta.append(file_time)
    
    datetime_series = pd.to_datetime(time_delta, format='%Y%m%d%H%M').to_numpy()
    flow_array = np.array(flow_array)
    return flow_array, datetime_series

def plot_data(args: argparse.Namespace) -> None:
    output_dir = Path(args.troute_output)
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True) 
    fid_list = args.feature_id
    print("Reading t-route output")
    flow, time_delta = get_fid_data(output_dir=output_dir, fid_list=fid_list)
    print("Reading RnR output")
    # input_flow_cfs, input_time_delta = get_rfc_data(input_dir=input_dir, fid_list=fid_list)
    flow_cfs = flow * 35.3147  # converting to cfs
    # input_flow_cfs = input_flow * 35.3147  # converting to cfs

    for idx, fid in enumerate(fid_list):
        _flow = flow_cfs[:, idx]
        # _input_flow = input_flow_cfs[:, idx]

        # sorted_indices = np.argsort(input_time_delta)
        # sorted_time = input_time_delta[sorted_indices]
        # sorted_flow = _input_flow[sorted_indices]
        fig, ax = plt.subplots(figsize=(12, 8), sharex=True)
        plt.plot(time_delta, _flow, c="tab:blue", label="T-Route Flow")
        plt.xlabel("timedelta64[ns]")
        plt.ylabel("discharge [cfs]")
        plt.title(f"T-Route Routed flow for feature: {fid}")
        plt.legend()

        # ax2.plot(sorted_time, sorted_flow, c="k", label="RFC Flow")
        # ax2.set_xlabel("timedelta64[ns]")
        # ax2.set_ylabel("discharge [cfs]")
        # ax2.set_title(f"RFC flow for feature: {fid}")
        # ax2.legend()

        # plt.gcf().autofmt_xdate()  # Rotate and align the tick labels
        
        # Adjust the layout and display the plot
        plt.tight_layout()
        # plt.show()
        plt.savefig(plot_dir / f"RFC_test_output_{fid}.png")
        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A script to create the RFC domain files for T-Route")
    parser.add_argument("troute_output", type=str, help="The directory where t-route output is saved")
    parser.add_argument(
        '--feature_id', 
        nargs='+', 
        help='A list of feature IDs for us to plot'
    )    
    args = parser.parse_args()

    plot_data(args)
