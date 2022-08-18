from abc import ABC, abstractmethod
import troute.nhd_io as nhd_io
import pandas as pd
import pathlib
import time

#FIXME parameterize into construciton
showtiming = True
verbose = True

def build_data_assimilation_csv(data_assimilation_parameters):

    return nhd_io.get_usgs_df_from_csv(
        data_assimilation_parameters["data_assimilation_csv"],
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
    )

def build_data_assimilation_folder(data_assimilation_parameters, run_parameters):

    usgs_timeslices_folder = pathlib.Path(
        data_assimilation_parameters["data_assimilation_timeslices_folder"],
    ).resolve()
    if "data_assimilation_filter" in data_assimilation_parameters:
        da_glob_filter = data_assimilation_parameters["data_assimilation_filter"]
        usgs_files = sorted(usgs_timeslices_folder.glob(da_glob_filter))
    elif "usgs_timeslice_files" in data_assimilation_parameters:
        usgs_files = data_assimilation_parameters.get("usgs_timeslice_files", None)
        usgs_files = [usgs_timeslices_folder.joinpath(f) for f in usgs_files]
    else:
        print("No Files Found for DA")
        # TODO: Handle this with a real exception

    return nhd_io.get_usgs_from_time_slices_folder(
        data_assimilation_parameters["wrf_hydro_da_channel_ID_crosswalk_file"],
        usgs_files,
        data_assimilation_parameters.get("qc_threshold", 1),
        data_assimilation_parameters.get("data_assimilation_interpolation_limit", 59),
        #FIXME/TODO collapse these into da parameters
        run_parameters["dt"],
        run_parameters["t0"],
    )

class DataAssimilation(ABC):
    """
    
    """

    @property
    @abstractmethod
    def usgs_df(self):
        pass

class NudgingDA(DataAssimilation):
    """
    
    """
    slots = ["_usgs_df", "_lastobs_df", "_da_params"]
    def __init__(self, data_assimilation_parameters, run_parameters):
        data_assimilation_csv = data_assimilation_parameters.get(
            "data_assimilation_csv", None
        )
        data_assimilation_folder = data_assimilation_parameters.get(
            "data_assimilation_timeslices_folder", None
        )
        # TODO: Copy comments from nhd_network_utilitys `build_data_assimilation_lastobs`
        lastobs_file = data_assimilation_parameters.get("wrf_hydro_lastobs_file", None)
        lastobs_start = data_assimilation_parameters.get(
            "wrf_hydro_lastobs_lead_time_relative_to_simulation_start_time", 0
        )
        lastobs_type = data_assimilation_parameters.get("wrf_lastobs_type", "error-based")
        lastobs_crosswalk_file = data_assimilation_parameters.get(
            "wrf_hydro_da_channel_ID_crosswalk_file", None
        )

        self._da_params = {}

        if data_assimilation_csv or data_assimilation_folder or lastobs_file:
            self._da_params["da_decay_coefficient"] = data_assimilation_parameters.get("da_decay_coefficient", 120)
            # TODO: Add parameters here for interpolation length (14/59), QC threshold (1.0)

            if showtiming:
                start_time = time.time()
            if verbose:
                print("creating usgs time_slice data array ...")
            self._last_obs_df = nhd_io.build_lastobs_df(
                lastobs_file,
                lastobs_crosswalk_file,
                lastobs_type,  # TODO: Confirm that we are using this; delete it if not.
                lastobs_start,
            )
            if data_assimilation_csv:
                self._usgs_df = build_data_assimilation_csv(data_assimilation_parameters)
            elif data_assimilation_folder:
                self._usgs_df = build_data_assimilation_folder(data_assimilation_parameters, run_parameters)

            if not self._last_obs_df.index.empty:
                if not self._usgs_df.empty and not self._usgs_df.index.equals(self._last_obs_df.index):
                    print("USGS Dataframe Index Does Not Match Last Observations Dataframe Index")
                    self._usgs_df = self._usgs_df.loc[self._last_obs_df.index]
            if verbose:
                print("usgs array complete")
            if showtiming:
                print("... in %s seconds." % (time.time() - start_time))
        else:
            self._last_obs_df = pd.DataFrame()
            self._usgs_df = pd.DataFrame()

    @property
    def asssimilation_parameters(self):
        return self._da_params
    
    @property
    def last_obs(self):
        return self._last_obs_df

    @property
    def usgs_df(self):
        return self._usgs_df