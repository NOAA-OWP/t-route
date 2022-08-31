from  .AbstractNetwork import AbstractNetwork
import xarray as xr
import pathlib
from collections import defaultdict
import troute.nhd_io as nhd_io
import pandas as pd
import time

from nwm_routing.preprocess import (
    nwm_network_preprocess,
    nwm_initial_warmstate_preprocess,
    nwm_forcing_preprocess,
    unpack_nwm_preprocess_data,
    )
import troute.nhd_network_utilities_v02 as nnu
import logging
LOG = logging.getLogger('')

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
    __slots__ = ["break_network_at_waterbodies", "waterbody_type_specified", "link_lake_crosswalk", "link_gage_df", 
                 "usgs_lake_gage_crosswalk", "usace_lake_gage_crosswalk", "diffusive_network_data", "topobathy", 
                 "refactored_diffusive_domain", "refactored_reaches", "unrefactored_topobathy", "segment_index", 
                 "t0", "lastobs_df", "da_parameter_dict", "run_sets", "da_sets", "usgs_df", "reservoir_usgs_df", 
                 "reservoir_usgs_param_df", "reservoir_usace_df", "reservoir_usace_param_df"]
    #def __init__(self, supernetwork_parameters, waterbody_parameters=None, restart_parameters=None, forcing_parameters=None, verbose=False, showtiming=False, layer_string=None, driver_string=None,):
    def __init__(
                 self, 
                 supernetwork_parameters, 
                 waterbody_parameters, 
                 restart_parameters, 
                 preprocessing_parameters, 
                 forcing_parameters,
                 compute_parameters, 
                 hybrid_parameters,
                 data_assimilation_parameters,
                 output_parameters,
                 cpu_pool,
                 task_times,
                 verbose=False,
                 showtiming=False,              
                ):
        """
        
        """
        if showtiming:
            network_start_time = time.time() 

        if preprocessing_parameters.get('use_preprocessed_data', False): 
        
            # get data from pre-processed file
            (
                self._connections, #connections,
                self._dataframe, #param_df,
                self._waterbody_connections, #wbody_conn,
                self._waterbody_df, #waterbodies_df,
                self._waterbody_types_df, #waterbody_types_df,
                self.break_network_at_waterbodies,
                self.waterbody_type_specified,
                self.link_lake_crosswalk,
                self._independent_networks, #independent_networks,
                self._reaches_by_tw, #reaches_bytw,
                self._reverse_network, #rconn,
                self.link_gage_df,
                self.usgs_lake_gage_crosswalk, 
                self.usace_lake_gage_crosswalk,
                self.diffusive_network_data,
                self.topobathy,
                self.refactored_diffusive_domain,
                self.refactored_reaches,
                self.unrefactored_topobathy,
            ) = unpack_nwm_preprocess_data(
                preprocessing_parameters
            )
        else:
        
        # build data objects from scratch
            (
                self._connections, #connections,
                self._dataframe, #param_df,
                self._waterbody_connections, #wbody_conn,
                self._waterbody_df, #waterbodies_df,
                self._waterbody_types_df, #waterbody_types_df,
                self.break_network_at_waterbodies,
                self.waterbody_type_specified,
                self.link_lake_crosswalk,                
                self._independent_networks, #independent_networks,
                self._reaches_by_tw, #reaches_bytw,
                self._reverse_network, #rconn,                
                self.link_gage_df,
                self.usgs_lake_gage_crosswalk, 
                self.usace_lake_gage_crosswalk,
                self.diffusive_network_data,
                self.topobathy,
                self.refactored_diffusive_domain,
                self.refactored_reaches,
                self.unrefactored_topobathy,
            ) = nwm_network_preprocess(
                supernetwork_parameters,
                waterbody_parameters,
                preprocessing_parameters,
                compute_parameters,
                data_assimilation_parameters,
            )

        # list of all segments in the domain (MC + diffusive)
        self.segment_index = self._dataframe.index  #param_df.index
        if self.diffusive_network_data:
            for tw in self.diffusive_network_data:
                self.segment_index = self.segment_index.append(
                    pd.Index(self.diffusive_network_data[tw]['mainstem_segs'])
                ) 

        # TODO: This function modifies one of its arguments (waterbodies_df), which is somewhat poor practice given its otherwise functional nature. Consider refactoring
        (
            self._waterbody_df, #waterbodies_df, 
            self._q0, #q0, 
            self.t0, 
            self.lastobs_df, 
            self.da_parameter_dict 
        ) = nwm_initial_warmstate_preprocess(
            self.break_network_at_waterbodies,
            restart_parameters,
            data_assimilation_parameters,
            self.segment_index,
            self._waterbody_df, #waterbodies_df,
            self.link_lake_crosswalk,
        )

        if showtiming:
            network_end_time = time.time()
            task_times['network_time'] = network_end_time - network_start_time
        
        if showtiming:
            ic_end_time = time.time()
            task_times['initial_condition_time'] += ic_end_time - network_end_time

        # Create run_sets: sets of forcing files for each loop
        self.run_sets = nnu.build_forcing_sets(forcing_parameters, self.t0)

        # Create da_sets: sets of TimeSlice files for each loop
        if "data_assimilation_parameters" in compute_parameters:
            self.da_sets = nnu.build_da_sets(data_assimilation_parameters, self.run_sets, self.t0)
        
        (
            _, #self._qlateral, #qlats, 
            self.usgs_df, 
            self.reservoir_usgs_df, 
            self.reservoir_usgs_param_df,
            self.reservoir_usace_df,
            self.reservoir_usace_param_df,
            self._coastal_boundary_depth
        ) = nwm_forcing_preprocess(
            self.run_sets[0],
            forcing_parameters,
            hybrid_parameters,
            self.da_sets[0] if data_assimilation_parameters else {},
            data_assimilation_parameters,
            self.break_network_at_waterbodies,
            self.segment_index,
            self.link_gage_df,
            self.usgs_lake_gage_crosswalk, 
            self.usace_lake_gage_crosswalk,
            self.link_lake_crosswalk,
            self.lastobs_df.index,
            cpu_pool,
            self.t0,
        )

        if showtiming:
            forcing_end_time = time.time()
            task_times['forcing_time'] += forcing_end_time - ic_end_time

  
    def _handle_args_v03(argv):
        '''
        Handle command line input argument - filepath of configuration file
        '''
        parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            "-f",
            "--custom-input-file",
            dest="custom_input_file",
            help="Path of a .yaml or .json file containing model configuration parameters. See doc/v3_doc.yaml",
        )
        return parser.parse_args(argv)
        
        
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
   

    def qlateral(self, t0, forcing_parameters, run, cpu_pool, segment_index):
        start_time = time.time()
        LOG.info("Creating a DataFrame of lateral inflow forcings ...")

        # TODO: find a better way to deal with these defaults and overrides.
        dt                           = forcing_parameters.get("dt", None)
        qts_subdivisions             = forcing_parameters.get("qts_subdivisions", None)
        qlat_input_folder            = forcing_parameters.get("qlat_input_folder", None)
        qlat_file_index_col          = forcing_parameters.get("qlat_file_index_col", "feature_id")
        qlat_file_value_col          = forcing_parameters.get("qlat_file_value_col", "q_lateral")
        qlat_file_gw_bucket_flux_col = forcing_parameters.get("qlat_file_gw_bucket_flux_col", "qBucket")
        qlat_file_terrain_runoff_col = forcing_parameters.get("qlat_file_terrain_runoff_col", "qSfcLatRunoff")

        # TODO: find a better way to deal with these defaults and overrides.
        run["t0"]                           = t0
        run["dt"]                           = run.get("dt", dt)
        run["qts_subdivisions"]             = run.get("qts_subdivisions", qts_subdivisions)
        run["qlat_input_folder"]            = run.get("qlat_input_folder", qlat_input_folder)
        run["qlat_file_index_col"]          = run.get("qlat_file_index_col", qlat_file_index_col)
        run["qlat_file_value_col"]          = run.get("qlat_file_value_col", qlat_file_value_col)
        run["qlat_file_gw_bucket_flux_col"] = run.get("qlat_file_gw_bucket_flux_col", qlat_file_gw_bucket_flux_col)
        run["qlat_file_terrain_runoff_col"] = run.get("qlat_file_terrain_runoff_col", qlat_file_terrain_runoff_col)
        
        self._qlateral = nnu.build_qlateral_array(
            run,
            cpu_pool,
            segment_index,
            )

        LOG.debug("lateral inflow DataFrame creation complete in %s seconds." \
            % (time.time() - start_time))
        return self._qlateral
