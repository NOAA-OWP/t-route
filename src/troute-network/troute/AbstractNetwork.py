from abc import ABC, abstractmethod
from functools import partial
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

import os
import pathlib
import time
import logging
import pyarrow as pa
import pyarrow.parquet as pq
import xarray as xr

from troute.nhd_network import extract_connections, replace_waterbodies_connections, reverse_network, reachable_network, split_at_waterbodies_and_junctions, split_at_junction, dfs_decomposition
from troute.nhd_network_utilities_v02 import organize_independent_networks
import troute.nhd_io as nhd_io 
from .AbstractRouting import MCOnly, MCwithDiffusive, MCwithDiffusiveNatlXSectionNonRefactored, MCwithDiffusiveNatlXSectionRefactored

LOG = logging.getLogger('')

class AbstractNetwork(ABC):
    """
    
    """
    __slots__ = ["_dataframe", "_waterbody_connections", "_gages",  
                "_terminal_codes", "_connections", "_waterbody_df", 
                "_waterbody_types_df", "_waterbody_type_specified", "_link_gage_df",
                "_independent_networks", "_reaches_by_tw", "_flowpath_dict",
                "_reverse_network", "_q0", "_t0", "_link_lake_crosswalk",
                "_usgs_lake_gage_crosswalk", "_usace_lake_gage_crosswalk", "_rfc_lake_gage_crosswalk",
                "_qlateral", "_break_segments", "_segment_index", "_coastal_boundary_depth_df",
                "supernetwork_parameters", "waterbody_parameters","data_assimilation_parameters",
                "restart_parameters", "compute_parameters", "forcing_parameters",
                "hybrid_parameters", "preprocessing_parameters",
                "verbose", "showtiming", "break_points", "_routing"]
    
    def __init__(self,):

        self._independent_networks = None
        self._reverse_network = None
        self._reaches_by_tw = None
        self._q0 = None
        self._t0 = None
        self._qlateral = None
        self._link_gage_df = None
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
        break_network_at_waterbodies = self.waterbody_parameters.get("break_network_at_waterbodies", False)        
        streamflow_da = self.data_assimilation_parameters.get('streamflow_da', False)
        break_network_at_gages       = False       
        if streamflow_da:
            break_network_at_gages   = streamflow_da.get('streamflow_nudging', False)
        self.break_points            = {"break_network_at_waterbodies": break_network_at_waterbodies,
                                        "break_network_at_gages": break_network_at_gages}

        self._break_segments = set()

        if self.break_points["break_network_at_waterbodies"]:
            self._break_segments = self._break_segments | set(self.waterbody_connections.values())
        if self.break_points["break_network_at_gages"]:
            self._break_segments = self._break_segments | set(self.gages.get('gages',{}).keys())
        
        self.initialize_routing_scheme()

        self.create_independent_networks()

        self.initial_warmstate_preprocess()


    def assemble_forcings(self, run,):
        """
        Assemble model forcings. Forcings include hydrological lateral inflows (qlats)
        and coastal boundary depths for hybrid runs
        
        Aguments
        --------
        - run                          (dict): List of forcing files pertaining to a 
                                               single run-set

        Returns
        -------
        
        Notes
        -----
        
        """
    
        # Unpack user-specified forcing parameters
        dt                           = self.forcing_parameters.get("dt", None)
        qts_subdivisions             = self.forcing_parameters.get("qts_subdivisions", None)
        qlat_input_folder            = self.forcing_parameters.get("qlat_input_folder", None)
        qlat_file_index_col          = self.forcing_parameters.get("qlat_file_index_col", "feature_id")
        qlat_file_value_col          = self.forcing_parameters.get("qlat_file_value_col", "q_lateral")
        qlat_file_gw_bucket_flux_col = self.forcing_parameters.get("qlat_file_gw_bucket_flux_col", "qBucket")
        qlat_file_terrain_runoff_col = self.forcing_parameters.get("qlat_file_terrain_runoff_col", "qSfcLatRunoff")

    
        # TODO: find a better way to deal with these defaults and overrides.
        run["t0"]                           = run.get("t0", self.t0)
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
        start_time = time.time()
        LOG.info("Creating a DataFrame of lateral inflow forcings ...")

        self.build_qlateral_array(
            run,
        )
        
        LOG.debug(
            "lateral inflow DataFrame creation complete in %s seconds." \
                % (time.time() - start_time)
                )

        #---------------------------------------------------------------------
        # Assemble coastal coupling data [WIP]
        #---------------------------------------------------------------------
        # Run if coastal_boundary_depth_df has not already been created:
        if self._coastal_boundary_depth_df.empty:
            coastal_boundary_elev_files = self.forcing_parameters.get('coastal_boundary_input_file', None) 
            coastal_boundary_domain_files = self.hybrid_parameters.get('coastal_boundary_domain', None)    
            
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
            else:
                self._coastal_boundary_depth_df = pd.DataFrame()            

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
            if self.showtiming:
                start_time = time.time()
            if self.verbose:
                print("organizing connections into reaches ...")

            self._independent_networks = reachable_network(self.reverse_network)
            
            if self.verbose:
                print("reach organization complete")
            if self.showtiming:
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
    def link_lake_crosswalk(self):
        return self._link_lake_crosswalk

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
    def segment_index(self):
        """
            Segment IDs of all reaches in parameter dataframe
            and diffusive domain.
        """
        # list of all segments in the domain (MC + diffusive)
        self._segment_index = self.dataframe.index
        if self._routing.diffusive_network_data:
            for tw in self._routing.diffusive_network_data:
                self._segment_index = self._segment_index.append(
                    pd.Index(self._routing.diffusive_network_data[tw]['mainstem_segs'])
                )
        return self._segment_index
    
    @property
    def link_gage_df(self):
        if self._link_gage_df is None:
            self._link_gage_df = pd.DataFrame.from_dict(self._gages)
            self._link_gage_df.index.name = 'link'
        return self._link_gage_df.copy()

    @property
    def usgs_lake_gage_crosswalk(self):
        return self._usgs_lake_gage_crosswalk
    
    @property
    def usace_lake_gage_crosswalk(self):
        return self._usace_lake_gage_crosswalk
    
    @property
    def rfc_lake_gage_crosswalk(self):
        return self._rfc_lake_gage_crosswalk

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

    @property
    def diffusive_network_data(self):
        return self._routing.diffusive_network_data
    
    @property
    def topobathy_df(self):
        return self._routing.topobathy_df
    
    @property
    def refactored_diffusive_domain(self):
        return self._routing.refactored_diffusive_domain
    
    @property
    def refactored_reaches(self):
        return self._routing.refactored_reaches
    
    @property
    def unrefactored_topobathy_df(self):
        return self._routing.unrefactored_topobathy_df

        
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
            
    def initialize_routing_scheme(self,):
        '''
        
        '''
        # Get user inputs from configuration file
        run_hybrid = self.hybrid_parameters.get("run_hybrid_routing", False)
        use_topobathy = self.hybrid_parameters.get('use_natl_xsections', False)
        run_refactored = self.hybrid_parameters.get('run_refactored_network', False)

        routing_type = [run_hybrid, use_topobathy, run_refactored]

        _routing_scheme_map = {
            MCOnly: [False, False, False],
            MCwithDiffusive: [True, False, False],
            MCwithDiffusiveNatlXSectionNonRefactored: [True, True, False],
            MCwithDiffusiveNatlXSectionRefactored: [True, True, True],
            }
        
        # Default to MCOnly routing
        routing_scheme = MCOnly

        # Check user input to determine the routing scheme
        for key, value in _routing_scheme_map.items():
            if value==routing_type:
                routing_scheme = key
        
        routing = routing_scheme(self.hybrid_parameters)

        (
            self._dataframe,
            self._connections
        ) = routing.update_routing_domain(self.dataframe, self.connections)

        self._routing = routing
    
    def create_independent_networks(self,):

        LOG.info("organizing connections into reaches ...")
        start_time = time.time() 
        gage_break_segments = set()
        wbody_break_segments = set()
        
        break_network_at_waterbodies = self.waterbody_parameters.get(
            "break_network_at_waterbodies", False
        )
        
        # if streamflow DA, then break network at gages
        #TODO update to work with HYFeatures, need to determine how we'll do DA...
        break_network_at_gages = False
        
        if break_network_at_waterbodies:
            wbody_break_segments = wbody_break_segments.union(self._waterbody_connections.values())
            
        if break_network_at_gages:
            gage_break_segments = gage_break_segments.union(self.gages['gages'].keys())
    
        (
            self._independent_networks, 
            self._reaches_by_tw, 
            self._reverse_network
            ) = organize_independent_networks(
                self.connections,
                wbody_break_segments,
                gage_break_segments,
                )
        
        LOG.debug("reach organization complete in %s seconds." % (time.time() - start_time))

    def initial_warmstate_preprocess(self,):

        '''
        Assemble model initial condition data:
            - waterbody inital states (outflow and pool elevation)
            - channel initial states (flow and depth)
            - initial time
        
        Arguments
        ---------
        
        Returns
        -------
        
        Notes
        -----
        '''

        restart_parameters = self.restart_parameters

        # generalize waterbody ID's to be used with any network
        index_id = self.waterbody_dataframe.index.names[0]

        #----------------------------------------------------------------------------
        # Assemble waterbody initial states (outflow and pool elevation
        #----------------------------------------------------------------------------
        if self.break_points['break_network_at_waterbodies']:

            start_time = time.time()
            LOG.info("setting waterbody initial states ...")

            # if a lite restart file is provided, read initial states from it.
            if restart_parameters.get("lite_waterbody_restart_file", None):
                
                waterbodies_initial_states_df, _ = nhd_io.read_lite_restart(
                    restart_parameters['lite_waterbody_restart_file']
                )
                
            # read waterbody initial states from WRF-Hydro type restart file
            elif restart_parameters.get("wrf_hydro_waterbody_restart_file", None):
                waterbodies_initial_states_df = nhd_io.get_reservoir_restart_from_wrf_hydro(
                    restart_parameters["wrf_hydro_waterbody_restart_file"],
                    restart_parameters["wrf_hydro_waterbody_ID_crosswalk_file"],
                    restart_parameters.get("wrf_hydro_waterbody_ID_crosswalk_file_field_name", index_id),
                    restart_parameters["wrf_hydro_waterbody_crosswalk_filter_file"],
                    restart_parameters.get(
                        "wrf_hydro_waterbody_crosswalk_filter_file_field_name",
                        'NHDWaterbodyComID'
                    ),
                )
            
            # if no restart file is provided, default initial states
            else:
                # TODO: Consider adding option to read cold state from route-link file
                waterbodies_initial_ds_flow_const = 0.0
                waterbodies_initial_depth_const = -1e9
                # Set initial states from cold-state
                waterbodies_initial_states_df = pd.DataFrame(
                    0,
                    index=self.waterbody_dataframe.index,
                    columns=[
                        "qd0",
                        "h0",
                    ],
                    dtype="float32",
                )
                # TODO: This assignment could probably by done in the above call
                waterbodies_initial_states_df["qd0"] = waterbodies_initial_ds_flow_const
                waterbodies_initial_states_df["h0"] = waterbodies_initial_depth_const
                waterbodies_initial_states_df["index"] = range(
                    len(waterbodies_initial_states_df)
                )

            self._waterbody_df = pd.merge(
                self.waterbody_dataframe, waterbodies_initial_states_df, on=index_id
            )

            LOG.debug(
                "waterbody initial states complete in %s seconds."\
                % (time.time() - start_time))
            start_time = time.time()

        #----------------------------------------------------------------------------
        # Assemble channel initial states (flow and depth)
        # also establish simulation initialization timestamp
        # 3 Restart Options:
        #   1. From t-route generated lite restart file (network agnostic)
        #   2. From wrf_hydro_restart file (valid for NHDNetwork only)
        #   3. Cold start, requires user specified start datetime
        #----------------------------------------------------------------------------    
        start_time = time.time()
        LOG.info("setting channel initial states ...")

        # if lite restart file is provided, the read channel initial states from it
        if restart_parameters.get("lite_channel_restart_file", None):
            self._q0, self._t0 = nhd_io.read_lite_restart(
                restart_parameters['lite_channel_restart_file']
            )
        
        elif restart_parameters.get("wrf_hydro_channel_restart_file", None):
            self._q0 = nhd_io.get_channel_restart_from_wrf_hydro(
                restart_parameters["wrf_hydro_channel_restart_file"],
                restart_parameters["wrf_hydro_channel_ID_crosswalk_file"],
                restart_parameters.get("wrf_hydro_channel_ID_crosswalk_file_field_name", 'link'),
                restart_parameters.get("wrf_hydro_channel_restart_upstream_flow_field_name", 'qlink1'),
                restart_parameters.get("wrf_hydro_channel_restart_downstream_flow_field_name", 'qlink2'),
                restart_parameters.get("wrf_hydro_channel_restart_depth_flow_field_name", 'hlink'),
                )

            t0_str = nhd_io.get_param_str(
                    restart_parameters["wrf_hydro_channel_restart_file"], 
                    "Restart_Time"
                )
            
            # convert timestamp from string to datetime
            self._t0 = datetime.strptime(t0_str, "%Y-%m-%d_%H:%M:%S")
        
        else:
            # Set cold initial state
            # assume to be zero
            # 0, index=connections.keys(), columns=["qu0", "qd0", "h0",], dtype="float32"
            self._q0 = pd.DataFrame(
                0, index=self.segment_index, columns=["qu0", "qd0", "h0"], dtype="float32",
                )
            
            # get initial time from user inputs
            self._t0 = restart_parameters.get("start_datetime")

        LOG.debug(
            "channel initial states complete in %s seconds."\
            % (time.time() - start_time)
        )

    def build_forcing_sets(self,):

        forcing_parameters = self.forcing_parameters
        supernetwork_parameters = self.supernetwork_parameters

        run_sets           = forcing_parameters.get("qlat_forcing_sets", None)
        qlat_input_folder  = forcing_parameters.get("qlat_input_folder", None)
        nts                = forcing_parameters.get("nts", None)
        max_loop_size      = forcing_parameters.get("max_loop_size", 12)
        dt                 = forcing_parameters.get("dt", None)

        try:
            qlat_input_folder = pathlib.Path(qlat_input_folder)
            assert qlat_input_folder.is_dir() == True
        except TypeError:
            raise TypeError("Aborting simulation because no qlat_input_folder is specified in the forcing_parameters section of the .yaml control file.") from None
        except AssertionError:
            raise AssertionError("Aborting simulation because the qlat_input_folder:", qlat_input_folder,"does not exist. Please check the the nexus_input_folder variable is correctly entered in the .yaml control file") from None

        forcing_glob_filter = forcing_parameters["qlat_file_pattern_filter"]

        if forcing_glob_filter=="nex-*":
            print("Reformating qlat nexus files as hourly binary files...")
            binary_folder = forcing_parameters.get('binary_nexus_file_folder', None)
            qlat_files = qlat_input_folder.glob(forcing_glob_filter)

            #Check that directory/files specified will work
            if not binary_folder:
                raise(RuntimeError("No output binary qlat folder supplied in config"))
            elif not os.path.exists(binary_folder):
                raise(RuntimeError("Output binary qlat folder supplied in config does not exist"))
            elif len(list(pathlib.Path(binary_folder).glob('*.parquet'))) != 0:
                raise(RuntimeError("Output binary qlat folder supplied in config is not empty (already contains '.parquet' files)"))

            #Add tnx for backwards compatability
            qlat_files_list = list(qlat_files) + list(qlat_input_folder.glob('tnx*.csv'))
            #Convert files to binary hourly files, reset nexus input information
            qlat_input_folder, forcing_glob_filter = nex_files_to_binary(qlat_files_list, binary_folder)
            forcing_parameters["qlat_input_folder"] = qlat_input_folder
            forcing_parameters["qlat_file_pattern_filter"] = forcing_glob_filter
            
        # TODO: Throw errors if insufficient input data are available
        if run_sets:        
            #FIXME: Change it for hyfeature
            '''
            # append final_timestamp variable to each set_list
            qlat_input_folder = pathlib.Path(qlat_input_folder)
            for (s, _) in enumerate(run_sets):
                final_chrtout = qlat_input_folder.joinpath(run_sets[s]['qlat_files'
                        ][-1])
                final_timestamp_str = nhd_io.get_param_str(final_chrtout,
                        'model_output_valid_time')
                run_sets[s]['final_timestamp'] = \
                    datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')
            '''  
        elif qlat_input_folder:        
            # Construct run_set dictionary from user-specified parameters

            # get the first and seconded files from an ordered list of all forcing files
            qlat_input_folder = pathlib.Path(qlat_input_folder)
            all_files          = sorted(qlat_input_folder.glob(forcing_glob_filter))
            first_file         = all_files[0]
            second_file        = all_files[1]

            # Deduce the timeinterval of the forcing data from the output timestamps of the first
            # two ordered CHRTOUT files
            if forcing_glob_filter=="*.CHRTOUT_DOMAIN1":
                t1 = nhd_io.get_param_str(first_file, "model_output_valid_time")
                t1 = datetime.strptime(t1, "%Y-%m-%d_%H:%M:%S")
                t2 = nhd_io.get_param_str(second_file, "model_output_valid_time")
                t2 = datetime.strptime(t2, "%Y-%m-%d_%H:%M:%S")
            elif forcing_glob_filter.startswith('*NEXOUT'):
                t1_str = first_file.name.split('NEXOUT', 1)[0]
                t1 = datetime.strptime(t1_str, '%Y%m%d%H%M')
                t2_str = second_file.name.split('NEXOUT', 1)[0]
                t2 = datetime.strptime(t2_str, '%Y%m%d%H%M')
            else:
                df     = read_file(first_file)
                t1_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
                t1     = datetime.strptime(t1_str,"%Y-%m-%d_%H:%M:%S")
                df     = read_file(second_file)
                t2_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
                t2     = datetime.strptime(t2_str,"%Y-%m-%d_%H:%M:%S")
            
            dt_qlat_timedelta = t2 - t1
            dt_qlat = dt_qlat_timedelta.seconds

            # determine qts_subdivisions
            qts_subdivisions = dt_qlat / dt
            if dt_qlat % dt == 0:
                qts_subdivisions = int(dt_qlat / dt)
            # make sure that qts_subdivisions = dt_qlat / dt
            forcing_parameters['qts_subdivisions']= qts_subdivisions

            # the number of files required for the simulation
            nfiles = int(np.ceil(nts / qts_subdivisions))
            
            # list of forcing file datetimes
            #datetime_list = [t0 + dt_qlat_timedelta * (n + 1) for n in
            #                 range(nfiles)]
            # ** Correction ** Because qlat file at time t is constantly applied throughout [t, t+1],
            #               ** n + 1 should be replaced by n
            datetime_list = [self.t0 + dt_qlat_timedelta * (n) for n in
                            range(nfiles)]        
            datetime_list_str = [datetime.strftime(d, '%Y%m%d%H%M') for d in
                                datetime_list]

            # list of forcing files
            forcing_filename_list = [d_str + forcing_glob_filter[1:] for d_str in
                                    datetime_list_str]
            
            # check that all forcing files exist
            for f in forcing_filename_list:
                try:
                    J = pathlib.Path(qlat_input_folder.joinpath(f))     
                    assert J.is_file() == True
                except AssertionError:
                    raise AssertionError("Aborting simulation because forcing file", J, "cannot be not found.") from None
                    
            # build run sets list
            run_sets = []
            k = 0
            j = 0
            nts_accum = 0
            nts_last = 0
            while k < len(forcing_filename_list):
                run_sets.append({})

                if k + max_loop_size < len(forcing_filename_list):
                    run_sets[j]['qlat_files'] = forcing_filename_list[k:k
                        + max_loop_size]
                else:
                    run_sets[j]['qlat_files'] = forcing_filename_list[k:]

                nts_accum += len(run_sets[j]['qlat_files']) * qts_subdivisions
                if nts_accum <= nts:
                    run_sets[j]['nts'] = int(len(run_sets[j]['qlat_files'])
                                            * qts_subdivisions)
                else:
                    run_sets[j]['nts'] = int(nts - nts_last)

                final_qlat = qlat_input_folder.joinpath(run_sets[j]['qlat_files'][-1]) 
                if forcing_glob_filter=="*.CHRTOUT_DOMAIN1":           
                    final_timestamp_str = nhd_io.get_param_str(final_qlat,'model_output_valid_time')
                elif forcing_glob_filter.startswith('*NEXOUT'):
                    
                    final_timestamp_str = datetime.strptime(
                        final_qlat.name.split('NEXOUT', 1)[0],
                        "%Y%m%d%H%M"
                    ).strftime("%Y-%m-%d_%H:%M:%S")
                else:
                    df = read_file(final_qlat)
                    final_timestamp_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")           
                
                run_sets[j]['final_timestamp'] = \
                    datetime.strptime(final_timestamp_str, '%Y-%m-%d_%H:%M:%S')

                nts_last = nts_accum
                k += max_loop_size
                j += 1

        return run_sets

def nex_files_to_binary(nexus_files, binary_folder):
    for f in nexus_files:
        # read the csv file
        df = pd.read_csv(f, usecols=[1,2], names=['Datetime','qlat'])
        
        # convert and reformat datetime column
        df['Datetime']= pd.to_datetime(df['Datetime']).dt.strftime("%Y%m%d%H%M")

        # reformat the dataframe
        df['feature_id'] = get_id_from_filename(f)
        df = df.pivot(index="feature_id", columns="Datetime", values="qlat")
        df.columns.name = None

        for col in df.columns:
            table_new = pa.Table.from_pandas(df.loc[:, [col]])
            
            if not os.path.exists(f'{binary_folder}/{col}NEXOUT.parquet'):
                pq.write_table(table_new, f'{binary_folder}/{col}NEXOUT.parquet')
            
            else:
                table_old = pq.read_table(f'{binary_folder}/{col}NEXOUT.parquet')
                table = pa.concat_tables([table_old,table_new])
                pq.write_table(table, f'{binary_folder}/{col}NEXOUT.parquet')
    
    nexus_input_folder = binary_folder
    forcing_glob_filter = '*NEXOUT.parquet'

    return nexus_input_folder, forcing_glob_filter

def get_id_from_filename(file_name):
    id = os.path.splitext(file_name)[0].split('-')[1].split('_')[0]
    return int(id)

def read_file(file_name):
    extension = file_name.suffix
    if extension=='.csv':
        df = pd.read_csv(file_name)
    elif extension=='.parquet':
        df = pq.read_table(file_name).to_pandas().reset_index()
        df.index.name = None
    elif extension=='.nc':
        df = xr.open_dataset(file_name).to_pandas().reset_index()
        df.index.name = None
    return df
