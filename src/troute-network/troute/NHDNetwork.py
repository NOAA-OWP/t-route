 from  .AbstractNetwork import AbstractNetwork
import xarray as xr
import pathlib
from collections import defaultdict
import troute.nhd_io as nhd_io
import pandas as pd
import time

__showtiming__ = True #FIXME pass flag
__verbose__ = True #FIXME pass verbosity

def read_qlats(forcing_parameters, segment_index):
    # STEP 5: Read (or set) QLateral Inputs
    if __showtiming__:
        start_time = time.time()
    if __verbose__:
        print("creating qlateral array ...")
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", 1)
    nts = forcing_parameters.get("nts", 1)
    qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
    qlat_input_file = forcing_parameters.get("qlat_input_file", None)
    if qlat_input_folder:
        qlat_input_folder = pathlib.Path(qlat_input_folder)
        if "qlat_files" in forcing_parameters:
            qlat_files = forcing_parameters.get("qlat_files")
            qlat_files = [qlat_input_folder.joinpath(f) for f in qlat_files]
        elif "qlat_file_pattern_filter" in forcing_parameters:
            qlat_file_pattern_filter = forcing_parameters.get(
                "qlat_file_pattern_filter", "*CHRT_OUT*"
            )
            qlat_files = sorted(qlat_input_folder.glob(qlat_file_pattern_filter))

        qlat_file_index_col = forcing_parameters.get(
            "qlat_file_index_col", "feature_id"
        )
        qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", "q_lateral")

        qlat_df = nhd_io.get_ql_from_wrf_hydro_mf(
            qlat_files=qlat_files,
            #ts_iterator=ts_iterator,
            #file_run_size=file_run_size,
            index_col=qlat_file_index_col,
            value_col=qlat_file_value_col,
        )
        qlat_df = qlat_df[qlat_df.index.isin(segment_index)]
    elif qlat_input_file:
        qlat_df = nhd_io.get_ql_from_csv(qlat_input_file)
    else:
        qlat_const = forcing_parameters.get("qlat_const", 0)
        qlat_df = pd.DataFrame(
            qlat_const,
            index=segment_index,
            columns=range(nts // qts_subdivisions),
            dtype="float32",
        )
    
    max_col = 1 + nts // qts_subdivisions
    
    if len(qlat_df.columns) > max_col:
        qlat_df.drop(qlat_df.columns[max_col:], axis=1, inplace=True)
    
    if not segment_index.empty:
        qlat_df = qlat_df[qlat_df.index.isin(segment_index)]

    if __verbose__:
        print("qlateral array complete")
    if __showtiming__:
        print("... in %s seconds." % (time.time() - start_time))

    return qlat_df

class NHDNetwork(AbstractNetwork):
    """
    
    """
    
    def __init__(self, supernetwork_parameters, waterbody_parameters=None, restart_parameters=None, forcing_parameters=None, verbose=False, showtiming=False, layer_string=None, driver_string=None,):
        """
        
        """
        global __verbose__, __showtiming__
        __verbose__ = verbose
        __showtiming__ = showtiming
        if __verbose__:
            print("creating supernetwork connections set")
        if __showtiming__:
            start_time = time.time()
        geo_file_path = pathlib.Path(supernetwork_parameters["geo_file_path"])
        cols = supernetwork_parameters["columns"]
        terminal_code = supernetwork_parameters.get("terminal_code", 0)
        mask = supernetwork_parameters.get("mask_file_path", None)
        mask_layer = supernetwork_parameters.get("mask_layer_string", None)
        mask_key = supernetwork_parameters.get("mask_key", None)
        break_network_at_waterbodies = supernetwork_parameters.get(
        "break_network_at_waterbodies", False
        )
        break_network_at_gages = supernetwork_parameters.get(
            "break_network_at_gages", False
        )
        break_points = {"break_network_at_waterbodies": break_network_at_waterbodies,
                        "break_network_at_gages": break_network_at_gages}
        with xr.open_dataset(geo_file_path) as ds:
            self._dataframe = ds.to_dataframe()
        
        if mask:
            data_mask = nhd_io.read_mask(
                pathlib.Path(mask),
                layer_string=mask_layer,
            )
            self._dataframe = self._dataframe.filter(
                data_mask.iloc[:, mask_key], axis=0
            )

        self._waterbody_types_df = pd.DataFrame()
        self._waterbody_df = pd.DataFrame()
        #FIXME the base class constructor is finiky
        #as it requires the _dataframe, then sets some 
        #initial default properties...which, at the moment
        #are used by the subclass constructor.
        #So it needs to be called at just the right spot...
        super().__init__(cols, terminal_code, break_points)
        #Load waterbody/reservoir info
        if waterbody_parameters:
            #FIXME later, DO ALL LAKE PARAMS BETTER
            levelpool_params = waterbody_parameters.get('level_pool', None)
            if not levelpool_params:
                raise(RuntimeError("No supplied levelpool parameters in routing config"))
            
            lake_id = levelpool_params.get("level_pool_waterbody_id", "lake_id")
            self._waterbody_df = nhd_io.read_level_pool_waterbody_df(
                levelpool_params["level_pool_waterbody_parameter_file_path"],
                lake_id,
                self.waterbody_connections.values()
            )
            # Remove duplicate lake_ids and rows
            self._waterbody_df = (
                self._waterbody_df.reset_index()
                .drop_duplicates(subset=lake_id)
                .set_index(lake_id)
            )
            self._waterbody_df["qd0"] = 0.0
            self._waterbody_df["h0"] = -1e9
            hybrid_params = waterbody_parameters.get('hybrid_and_rfc', None)
            try:
                self._waterbody_types_df = nhd_io.read_reservoir_parameter_file(
                    hybrid_params["reservoir_parameter_file"],
                    lake_id,
                    self.waterbody_connections.values(),
                )
                # Remove duplicate lake_ids and rows
                self._waterbody_types_df = (
                self._waterbody_types_df.reset_index()
                .drop_duplicates(subset=lake_id)
                .set_index(lake_id)
                )
            except:
                self._waterbody_types_df = pd.DataFrame()
        if __verbose__:
            print("supernetwork connections set complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))
            
        if __verbose__:
            print("setting waterbody initial states ...")
        if __showtiming__:
            start_time = time.time()
        waterbodies_initial_states_df = pd.DataFrame(columns=self._waterbody_df.columns, index=self._waterbody_df.index)
        waterbodies_initial_states_df['qd0'] = 0.0 #Create the qd0 column
        waterbodies_initial_states_df['h0'] = -1e9 #Create  the h0 column
        if restart_parameters.get("wrf_hydro_waterbody_restart_file", None):
            waterbodies_initial_states_df = nhd_io.get_reservoir_restart_from_wrf_hydro(
                restart_parameters["wrf_hydro_waterbody_restart_file"],
                restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file"],
                restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file_field_name"],
                restart_parameters["wrf_hydro_waterbody_crosswalk_filter_file"],
                restart_parameters[
                    "wrf_hydro_waterbody_crosswalk_filter_file_field_name"
                ],
            )

        #Merge the data, keep the initial states by defining `how="right"`
        #self._waterbody_df = pd.merge(
        #    self._waterbody_df, waterbodies_initial_states_df, on="lake_id", suffixes=("_default", None)
        #)
        #NJF easier just to copy only what is needed?
        self._waterbody_df['qd0'] = waterbodies_initial_states_df['qd0'].astype('float32')
        self._waterbody_df['h0'] = waterbodies_initial_states_df['h0'].astype('float32')
        print(self._waterbody_df)
        if __verbose__:
            print("waterbody initial states complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))
            start_time = time.time()
        # STEP 4: Handle Channel Initial States
        if __showtiming__:
            start_time = time.time()
        if __verbose__:
            print("setting channel initial states ...")
        #Get channel restarts
        channel_restart_file = restart_parameters.get("channel_restart_file", None)

        wrf_hydro_channel_restart_file = restart_parameters.get(
            "wrf_hydro_channel_restart_file", None
        )

        if channel_restart_file:
            self._q0 = nhd_io.get_channel_restart_from_csv(channel_restart_file)
            self._q0 = self._q0[self._q0.index.isin(self._dataframe.index)]
            # TODO is this the same???
            #self._q0 = self._q0.loc[self._dataframe.index]
        elif wrf_hydro_channel_restart_file:

            self._q0 = nhd_io.get_channel_restart_from_wrf_hydro(
                restart_parameters["wrf_hydro_channel_restart_file"],
                restart_parameters["wrf_hydro_channel_ID_crosswalk_file"],
                restart_parameters["wrf_hydro_channel_ID_crosswalk_file_field_name"],
                restart_parameters["wrf_hydro_channel_restart_upstream_flow_field_name"],
                restart_parameters["wrf_hydro_channel_restart_downstream_flow_field_name"],
                restart_parameters["wrf_hydro_channel_restart_depth_flow_field_name"],
            )
            self._q0 = self._q0[self._q0.index.isin(self._dataframe.index)]
            # TODO is this the same???
            #self._q0 = self._q0.loc[self._dataframe.index]
            self.t0 = nhd_io.get_param_str(wrf_hydro_channel_restart_file, "Restart_Time")
        if __verbose__:
          print("channel initial states complete")
        if __showtiming__:
          print("... in %s seconds." % (time.time() - start_time))
        #Make sure waterbody parameter data is the correct type
        #self._waterbody_df = self._waterbody_df[['']]
        self._qlateral = read_qlats(forcing_parameters, self._dataframe.index)

    def extract_waterbody_connections(rows, target_col, waterbody_null=-9999):
        """Extract waterbody mapping from dataframe.
        TODO deprecate in favor of below property???"""
        return (
            rows.loc[rows[target_col] != waterbody_null, target_col].astype("int").to_dict()
        )

    @property
    def waterbody_connections(self):
        if self._waterbody_connections is None :
            self._waterbody_connections = self._dataframe.loc[
                self._dataframe["waterbody"] != self.waterbody_null, "waterbody"
            ].astype("int").to_dict()
        return self._waterbody_connections

    @property
    def gages(self):
        """
        
        """
        if self._gages is None and "gages" in self._dataframe.columns:
            self._gages = nhd_io.build_filtered_gage_df(self._dataframe[["gages"]])
        else:
            self._gages = {}
        return self._gages

    @property
    def waterbody_null(self):
        return -9999