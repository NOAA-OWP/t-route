from .AbstractNetwork import AbstractNetwork
import xarray as xr
import pathlib
from collections import defaultdict
import troute.nhd_io as nhd_io
import troute.nhd_preprocess as nhd_prep
import pandas as pd
import time
from datetime import datetime, timedelta

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
    __slots__ = ["_waterbody_type_specified", "_link_lake_crosswalk", "_link_gage_df",
                 "_usgs_lake_gage_crosswalk", "_usace_lake_gage_crosswalk",
                 "_diffusive_network_data", "_topobathy", "_refactored_diffusive_domain",
                 "_refactored_reaches", "_nonrefactored_topobathy", "_segment_index",
                 "_coastal_boundary_depth_df"]
    
    def __init__(self, supernetwork_parameters, waterbody_parameters=None, restart_parameters=None, forcing_parameters=None, compute_parameters=None, data_assimilation_parameters=None, preprocessing_parameters=None, verbose=False, showtiming=False, layer_string=None, driver_string=None,):
        """
        
        """
        global __verbose__, __showtiming__
        __verbose__ = verbose
        __showtiming__ = showtiming
        if __verbose__:
            print("creating supernetwork connections set")
        if __showtiming__:
            start_time = time.time()
        
        
        #------------------------------------------------
        # Preprocess network attributes
        #------------------------------------------------
        
        (self._connections,
         self._dataframe,
         self._waterbody_connections,
         self._waterbody_df,
         self._waterbody_types_df,
         break_network_at_waterbodies,
         self._waterbody_type_specified,
         self._link_lake_crosswalk,
         self._independent_networks,
         self._reaches_by_tw,
         self._reverse_network,
         self._link_gage_df,
         self._usgs_lake_gage_crosswalk, 
         self._usace_lake_gage_crosswalk,
         self._diffusive_network_data,
         self._topobathy,
         self._refactored_diffusive_domain,
         self._refactored_reaches,
         self._nonrefactored_topobathy
        ) = nhd_prep.build_nhd_network(
            supernetwork_parameters,
            waterbody_parameters,
            preprocessing_parameters,
            compute_parameters,
            data_assimilation_parameters
        )
        
        # list of all segments in the domain (MC + diffusive)
        self._segment_index = self._dataframe.index
        if self._diffusive_network_data:
            for tw in self._diffusive_network_data:
                self._segment_index = self._segment_index.append(
                    pd.Index(self._diffusive_network_data[tw]['mainstem_segs'])
                ) 

        '''
        #FIXME the base class constructor is finiky
        #as it requires the _dataframe, then sets some 
        #initial default properties...which, at the moment
        #are used by the subclass constructor.
        #So it needs to be called at just the right spot...
        cols = supernetwork_parameters.get(
            'columns', 
            {
                'key'       : 'link',
                'downstream': 'to',
                'dx'        : 'Length',
                'n'         : 'n',
                'ncc'       : 'nCC',
                's0'        : 'So',
                'bw'        : 'BtmWdth',
                'waterbody' : 'NHDWaterbodyComID',
                'gages'     : 'gages',
                'tw'        : 'TopWdth',
                'twcc'      : 'TopWdthCC',
                'alt'       : 'alt',
                'musk'      : 'MusK',
                'musx'      : 'MusX',
                'cs'        : 'ChSlp',
            }
        )
        terminal_code = supernetwork_parameters.get("terminal_code", 0)
        break_network_at_waterbodies = supernetwork_parameters.get(
        "break_network_at_waterbodies", False
        )
        break_network_at_gages = supernetwork_parameters.get(
            "break_network_at_gages", False
        )
        
        break_points = {"break_network_at_waterbodies": break_network_at_waterbodies,
                        "break_network_at_gages": break_network_at_gages}
        super().__init__(cols, terminal_code, break_points)
        '''
        
        if __verbose__:
            print("supernetwork connections set complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))
            
        #-----------------------------------------------------
        # Set initial waterbody and channel states
        #-----------------------------------------------------
        
        if __verbose__:
            print("setting waterbody and channel initial states ...")
        if __showtiming__:
            start_time = time.time()
        
        (self._waterbody_df,
         self._q0,
         self._t0,) = nhd_prep.nhd_initial_warmstate_preprocess(
            break_network_at_waterbodies,
            restart_parameters,
            data_assimilation_parameters,
            self._segment_index,
            self._waterbody_df,
            self._link_lake_crosswalk,
        )
        
        if __verbose__:
            print("waterbody and channel initial states complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))
            start_time = time.time()
        
        # Create empty dataframe for coastal_boundary_depth_df. This way we can check if
        # it exists, and only read in SCHISM data during 'assemble_forcings' if it doesn't
        self._coastal_boundary_depth_df = pd.DataFrame()
        

    def assemble_forcings(self, run, forcing_parameters, hybrid_parameters, cpu_pool):
        """
        Assembles model forcings for hydrological lateral inflows and coastal boundary 
        depths (hybrid simulations). Run this function after network initialization
        and after any iteration loop in main.
        """
        (self._qlateral, 
         self._coastal_boundary_depth_df
        ) = nhd_prep.nhd_forcing(
            run, 
            forcing_parameters, 
            hybrid_parameters, 
            self._segment_index, 
            cpu_pool,
            self._t0,
            self._coastal_boundary_depth_df,
        )
    
    def new_nhd_q0(self, run_results):
        """
        Prepare a new q0 dataframe with initial flow and depth to act as
        a warmstate for the next simulation chunk.
        """
        self._q0 = pd.concat(
            [
                pd.DataFrame(
                    r[1][:, [-3, -3, -1]], index=r[0], columns=["qu0", "qd0", "h0"]
                )
                for r in run_results
            ],
            copy=False,
        )
    
    def update_waterbody_water_elevation(self):
        """
        Update the starting water_elevation of each lake/reservoir
        with depth values from q0
        """
        self._waterbody_df.update(self._q0)
    
    def new_t0(self, dt, nts):
        """
        Update t0 value for next loop iteration
        """
        self._t0 += timedelta(seconds = dt * nts)

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

    @property
    def wbody_conn(self):
        return self._waterbody_connections
    
    @property
    def waterbody_type_specified(self):
        return self._waterbody_type_specified
    
    @property
    def link_lake_crosswalk(self):
        return self._link_lake_crosswalk
    
    @property
    def link_gage_df(self):
        return self._link_gage_df
    
    @property
    def usgs_lake_gage_crosswalk(self):
        return self._usgs_lake_gage_crosswalk
    
    @property
    def usace_lake_gage_crosswalk(self):
        return self._usace_lake_gage_crosswalk
    
    @property
    def diffusive_network_data(self):
        return self._diffusive_network_data
    
    @property
    def topobathy(self):
        return self._topobathy
    
    @property
    def refactored_diffusive_domain(self):
        return self._refactored_diffusive_domain
    
    @property
    def refactored_reaches(self):
        return self._refactored_reaches
    
    @property
    def nonrefactored_topobathy(self):
        return self._nonrefactored_topobathy
    
    @property
    def segment_index(self):
        return self._segment_index
    
    @property
    def coastal_boundary_depth_df(self):
        return self._coastal_boundary_depth_df