from abc import ABC, abstractmethod
from functools import partial
import pandas as pd
from datetime import datetime
import time

from troute.nhd_network import reverse_dict, extract_connections, replace_waterbodies_connections, reverse_network, reachable_network, split_at_waterbodies_and_junctions, split_at_junction, dfs_decomposition
import troute.nhd_io as nhd_io
import troute.abstractnetwork_preprocess as abs_prep 
__verbose__ = False
__showtiming__ = False

class AbstractNetwork(ABC):
    """
    
    """
    __slots__ = ["_dataframe", "_waterbody_connections", "_gages",  
                "_terminal_codes", "_connections", "_waterbody_df", 
                "_waterbody_types_df", "_waterbody_type_specified",
                "_independent_networks", "_reaches_by_tw", 
                "_reverse_network", "_q0", "_t0", 
                "_qlateral", "_break_segments", "_coastal_boundary_depth_df"]
    
    def __init__(
        self, 
        compute_parameters, 
        waterbody_parameters,
        restart_parameters,
        cols=None, 
        terminal_code=None, 
        break_points=None, 
        verbose=False, 
        showtiming=False
        ):
        
        global __verbose__, __showtiming__
        __verbose__ = verbose
        __showtiming__ = showtiming
        if cols:
            self._dataframe = self._dataframe[list(cols.values())]
            # Rename parameter columns to standard names: from route-link names
            #        key: "link"
            #        downstream: "to"
            #        dx: "Length"
            #        n: "n"  # TODO: rename to `manningn`
            #        ncc: "nCC"  # TODO: rename to `mannningncc`
            #        s0: "So"  # TODO: rename to `bedslope`
            #        bw: "BtmWdth"  # TODO: rename to `bottomwidth`
            #        waterbody: "NHDWaterbodyComID"
            #        gages: "gages"
            #        tw: "TopWdth"  # TODO: rename to `topwidth`
            #        twcc: "TopWdthCC"  # TODO: rename to `topwidthcc`
            #        alt: "alt"
            #        musk: "MusK"
            #        musx: "MusX"
            #        cs: "ChSlp"  # TODO: rename to `sideslope`
            self._dataframe = self._dataframe.rename(columns=reverse_dict(cols))
            self.set_index("key")
            self.sort_index()
        self._waterbody_connections = {}
        self._gages = None
        self._connections = None
        self._independent_networks = None
        self._reverse_network = None
        self._reaches_by_tw = None
        self._q0 = None
        self._t0 = None
        self._qlateral = None
        self._waterbody_type_specified = None
        #qlat_const = forcing_parameters.get("qlat_const", 0)
        #FIXME qlat_const
        """ Figure out a good way to default initialize to qlat_const/c
        qlat_const = 1.0
        self._qlateral = pd.DataFrame(
            qlat_const,
            index=self._dataframe.index,
            columns=range(nts // qts_subdivisions),
            dtype="float32",
        )
        """
        # there may be off-domain nodes that are not explicitly identified
        # but which are terminal (i.e., off-domain) as a result of a mask or some other
        # an interior domain truncation that results in a
        # otherwise valid node value being pointed to, but which is masked out or
        # being intentionally separated into another domain.
        self._terminal_codes = set(
            self._dataframe[
                ~self._dataframe["downstream"].isin(self._dataframe.index)
            ]["downstream"].values
        )
        
         # There can be an externally determined terminal code -- that's this value
        self._terminal_codes.add(terminal_code)
        
        self._break_segments = set()
        if break_points:
            if break_points["break_network_at_waterbodies"]:
                self._break_segments = self._break_segments | set(self.waterbody_connections.values())
            if break_points["break_network_at_gages"]:
                self._break_segments = self._break_segments | set(self.gages.values())
        
        self._connections = extract_connections(self._dataframe, 'downstream', self._terminal_codes)
        
        (
            self._dataframe,
            self._connections,
            self.diffusive_network_data,
            self.topobathy_df,
            self.refactored_diffusive_domain,
            self.refactored_reaches,
            self.unrefactored_topobathy_df
        ) = abs_prep.build_diffusive_domain(
            compute_parameters,
            self._dataframe,
            self._connections,
            )

        (
            self._independent_networks,
            self._reaches_by_tw, 
            self._reverse_network
        ) = abs_prep.create_independent_networks(
            waterbody_parameters, 
            self._connections, 
            self._waterbody_connections, 
            #gages, #TODO update how gages are provided when we figure out DA
            )
        
        (
            self._waterbody_df,
            self._q0, 
            self._t0
        ) = abs_prep.initial_warmstate_preprocess(
            break_points["break_network_at_waterbodies"],
            restart_parameters,
            self._dataframe.index,
            self._waterbody_df,
            )
        
    def assemble_forcings(self, run, forcing_parameters, hybrid_parameters, supernetwork_parameters, cpu_pool):
        """
        Assemble model forcings. Forcings include hydrological lateral inflows (qlats)
        and coastal boundary depths for hybrid runs
        
        Aguments
        --------
        - run                      (dict): List of forcing files pertaining to a 
                                    single run-set
        - forcing_parameters       (dict): User-input simulation forcing parameters
        - hybrid_parameters        (dict): User-input simulation hybrid parameters
        - supernetwork_parameters  (dict): User-input simulation supernetwork parameters
        - segment_index           (Int64): Reach segment ids
        - cpu_pool                  (int): Number of CPUs in the process-parallel pool

        Returns
        -------
        - qlats_df                 (Pandas DataFrame): Lateral inflow data, indexed by 
                                                    segment ID
        - coastal_bounary_depth_df (Pandas DataFrame): Coastal boundary water depths,
                                                    indexed by segment ID
        
        Notes
        -----
        
        """
    
        # Unpack user-specified forcing parameters
        dt                           = forcing_parameters.get("dt", None)
        qts_subdivisions             = forcing_parameters.get("qts_subdivisions", None)
        qlat_input_folder            = forcing_parameters.get("qlat_input_folder", None)
        qlat_file_index_col          = forcing_parameters.get("qlat_file_index_col", "feature_id")
        qlat_file_value_col          = forcing_parameters.get("qlat_file_value_col", "q_lateral")
        qlat_file_gw_bucket_flux_col = forcing_parameters.get("qlat_file_gw_bucket_flux_col", "qBucket")
        qlat_file_terrain_runoff_col = forcing_parameters.get("qlat_file_terrain_runoff_col", "qSfcLatRunoff")

    
        # TODO: find a better way to deal with these defaults and overrides.
        run["t0"]                           = run.get("t0", self.t0)
        run["nts"]                          = run.get("nts")
        run["dt"]                           = run.get("dt", dt)
        run["qts_subdivisions"]             = run.get("qts_subdivisions", qts_subdivisions)
        run["qlat_input_folder"]            = run.get("qlat_input_folder", qlat_input_folder)
        run["qlat_file_index_col"]          = run.get("qlat_file_index_col", qlat_file_index_col)
        run["qlat_file_value_col"]          = run.get("qlat_file_value_col", qlat_file_value_col)
        run["qlat_file_gw_bucket_flux_col"] = run.get("qlat_file_gw_bucket_flux_col", qlat_file_gw_bucket_flux_col)
        run["qlat_file_terrain_runoff_col"] = run.get("qlat_file_terrain_runoff_col", qlat_file_terrain_runoff_col)
        
        #---------------------------------------------------------------------------
        # Assemble lateral inflow data
        #---------------------------------------------------------------------------

        # Place holder, if reading qlats from a file use this.
        # TODO: add an option for reading qlat data from BMI/model engine
        from_file = True
        if from_file:
            self._qlateral = abs_prep.build_qlateral_array(
                run,
                cpu_pool,
                self._flowpath_dict,
                supernetwork_parameters, 
                self._dataframe.index,
            )

        #---------------------------------------------------------------------
        # Assemble coastal coupling data [WIP]
        #---------------------------------------------------------------------
        # Run if coastal_boundary_depth_df has not already been created:
        if self._coastal_boundary_depth_df.empty:
            coastal_boundary_elev_files = forcing_parameters.get('coastal_boundary_input_file', None) 
            coastal_boundary_domain_files = hybrid_parameters.get('coastal_boundary_domain', None)    
            
            if coastal_boundary_elev_files:
                #start_time = time.time()    
                #LOG.info("creating coastal dataframe ...")
                
                coastal_boundary_domain   = nhd_io.read_coastal_boundary_domain(coastal_boundary_domain_files)          
                self._coastal_boundary_depth_df = nhd_io.build_coastal_ncdf_dataframe(
                    coastal_boundary_elev_files,
                    coastal_boundary_domain,
                )
                    
                #LOG.debug(
                #    "coastal boundary elevation observation DataFrame creation complete in %s seconds." \
                #    % (time.time() - start_time)
                #)            

    def new_q0(self, run_results):
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
        return self._q0
    
    def update_waterbody_water_elevation(self):           
        """
        Update the starting water_elevation of each lake/reservoir
        with flow and depth values from q0
        """
        self._waterbody_df.update(self._q0)
        
    def new_t0(self, dt, nts):
        """
        Update t0 value for next loop iteration
        """
        self._t0 += timedelta(seconds = dt * nts)

    @property
    def network_break_segments(self):
        """
        """
        return self._break_segments

    @property
    def reverse_network(self):
        """
        
        """
        if self._reverse_network is None:
            self._reverse_network = reverse_network(self.connections)
        return self._reverse_network

    @property
    def independent_networks(self):
        """
        
        """
        if self._independent_networks is None:
            # STEP 2: Identify Independent Networks and Reaches by Network
            if __showtiming__:
                start_time = time.time()
            if __verbose__:
                print("organizing connections into reaches ...")

            self._independent_networks = reachable_network(self.reverse_network)
            
            if __verbose__:
                print("reach organization complete")
            if __showtiming__:
                print("... in %s seconds." % (time.time() - start_time))
        return self._independent_networks
    
    @property
    def reaches_by_tailwater(self):
        """
        
        """
        if self._reaches_by_tw is None:
            self._reaches_by_tw = {}
            for tw, net in self.independent_networks.items():
                if self.network_break_segments:
                    path_func = partial(
                        split_at_waterbodies_and_junctions, self.network_break_segments, net
                    )
                else:
                    path_func = partial(split_at_junction, net)

                self._reaches_by_tw[tw] = dfs_decomposition(net, path_func)
        return self._reaches_by_tw

    @property
    def waterbody_dataframe(self):
        return self._waterbody_df
    
    @property
    def waterbody_types_dataframe(self):
        return self._waterbody_types_df
    
    @property
    def waterbody_type_specified(self):
        if self._waterbody_type_specified is None:
            self._waterbody_type_specified = False
        return self._waterbody_type_specified

    @property
    def connections(self):
        if self._connections is None:
            self._connections = extract_connections(
                self._dataframe, "downstream", terminal_codes=self._terminal_codes
            )
        return self._connections

    @property
    def qlateral(self):
        """
        
        """
        return self._qlateral

    @property
    def q0(self):
        """
            Initial channel segment flow values
            If not set elsewhere, they are 0
        """
        if self._q0 is None:
            self._q0 =  pd.DataFrame(
            0, index=self._dataframe.index, columns=["qu0", "qd0", "h0"], dtype="float32",
            )
        return self._q0

    @property
    def t0(self):
        """
            Time 0 as a datetime object
            If not set elsewhere, it defaults to "2015-08-16_00:00:00"
        """
        if self._t0 is None:
            self._t0 = datetime.strptime("2015-08-16_00:00:00", "%Y-%m-%d_%H:%M:%S")
        return self._t0

    @t0.setter
    def t0(self, value):
        """
        
        """
        if isinstance(value, datetime):
            self._t0 = value
        else:
            self._t0 = datetime.strptime(value, "%Y-%m-%d_%H:%M:%S")

    @property
    @abstractmethod
    def waterbody_connections(self):
        pass

    @property
    @abstractmethod
    def waterbody_null(self):
        pass

    @property
    @abstractmethod
    def gages(self):
        pass

    @property
    def dataframe(self):
        return self._dataframe

    @property
    def terminal_codes(self):
        return self._terminal_codes
    
    @property
    def coastal_boundary_depth_df(self):
        """
        
        """
        return self._coastal_boundary_depth_df

    def set_synthetic_wb_segments(self, synthetic_wb_segments, synthetic_wb_id_offset):
        """
        
        """
        self._dataframe.reset_index(inplace=True) #reset index so key is now column
        # rename the current key column to key32
        key32_d = {"key":"key32"}
        self._dataframe = self._dataframe.rename(columns=key32_d)
        # create a key index that is int64
        # copy the links into the new column
        self._dataframe["key"] = self._dataframe.key32.astype("int64")
        # update the values of the synthetic reservoir segments
        fix_idx = self._dataframe.key.isin(set(synthetic_wb_segments))
        self._dataframe.loc[fix_idx,"key"] = (self._dataframe[fix_idx].key + synthetic_wb_id_offset).astype("int64")
        #Reset key to index
        self.set_index("key")
        self.sort_index()

    def replace_waterbodies(self):
        """

        """
        #Make sure to update held state self._connections
        #but pass the property self.connectionsto ensure it gets properly instansiated
        #in case it hasn't already been
        self._connections = replace_waterbodies_connections(
            self.connections, self._waterbody_connections
        )

    def set_index(self, key):
        """
        
        """
        #If the index name is already `key`, don't bother
        if self._dataframe.index.name != key:
            self._dataframe.set_index(key, inplace=True)
    
    def sort_index(self):
        """
        
        """
        self._dataframe = self._dataframe.sort_index()
    
    def drop(self, key, axis=1):
        """
            FIXME can be problematic to drop keys
            before certain properties are intialized...
        """
        self._dataframe.drop(key, axis=axis, inplace=True)

    def astype(self, type, columns=None):
        """
        
        """
        if columns:
            self._dataframe[columns] = self._dataframe[columns].astype(type)
        else:
            self._dataframe = self._dataframe.astype(type)
            
