from abc import ABC, abstractmethod
from functools import partial
import pandas as pd
import numpy as np
import pyarrow.parquet as pq
from datetime import datetime, timedelta
from joblib import delayed, Parallel
import netCDF4
import time
import logging
import pathlib

from troute.nhd_network import extract_connections, replace_waterbodies_connections, reverse_network, reachable_network, split_at_waterbodies_and_junctions, split_at_junction, dfs_decomposition
from troute.nhd_network_utilities_v02 import organize_independent_networks, build_channel_initial_state, build_refac_connections
import troute.nhd_io as nhd_io 

LOG = logging.getLogger('')

class AbstractNetwork(ABC):
    """
    
    """
    __slots__ = ["_dataframe", "_waterbody_connections", "_gages",  
                "_terminal_codes", "_connections", "_waterbody_df", 
                "_waterbody_types_df", "_waterbody_type_specified",
                "_independent_networks", "_reaches_by_tw", "_flowpath_dict",
                "_reverse_network", "_q0", "_t0", "_link_lake_crosswalk",
                "_qlateral", "_break_segments", "_coastal_boundary_depth_df",
                "_diffusive_network_data", "_topobathy_df", "_refactored_diffusive_domain",
                "_refactored_reaches", "_unrefactored_topobathy_df", "_segment_index",
                "supernetwork_parameters", "waterbody_parameters","data_assimilation_parameters",
                "restart_parameters", "compute_parameters", "verbose", "showtiming", "break_points"]
    
    def __init__(self,):

        self._independent_networks = None
        self._reverse_network = None
        self._reaches_by_tw = None
        self._q0 = None
        self._t0 = None
        self._qlateral = None
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
        self._break_segments = set()
        
        if self.break_points:
            if self.break_points["break_network_at_waterbodies"]:
                self._break_segments = self._break_segments | set(self.waterbody_connections.values())
            if self.break_points["break_network_at_gages"]:
                self._break_segments = self._break_segments | set(self.gages.get('gages').keys())
        
        #self.build_diffusive_domain(compute_parameters)

        self.create_independent_networks()

        self.initial_warmstate_preprocess())


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
            self.build_qlateral_array(
                run,
                cpu_pool,
                supernetwork_parameters,
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
        return self._waterbody_df.sort_index()
    
    @property
    def waterbody_types_dataframe(self):
        return self._waterbody_types_df.sort_index()
    
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
        if self.diffusive_network_data:
            for tw in self.diffusive_network_data:
                self._segment_index = self._segment_index.append(
                    pd.Index(self.diffusive_network_data[tw]['mainstem_segs'])
                )
        return self._segment_index
    
    @property
    def link_gage_df(self):
        link_gage_df = pd.DataFrame.from_dict(self._gages)
        link_gage_df.index.name = 'link'
        return link_gage_df

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
            
    
    def build_diffusive_domain(self, compute_parameters):
        """
        
        """
        hybrid_params = compute_parameters.get("hybrid_parameters", False)
        if hybrid_params:
            # switch parameters
            # if run_hybrid = False, run MC only
            # if run_hybrid = True, if use_topobathy = False, run MC+diffusive on RouteLink.nc
            #    "      "      "  , if use_topobathy = True,  if run_refactored_network = False, run MC+diffusive on original hydrofabric
            #    "      "      "  , if use_topobathy = True,  if run_refactored_network = True,  run MC+diffusive on refactored hydrofabric
            run_hybrid             = hybrid_params.get('run_hybrid_routing', False)
            use_topobathy          = hybrid_params.get('use_natl_xsections', False) 
            run_refactored         = hybrid_params.get('run_refactored_network', False)
            
            # file path parameters of non-refactored hydrofabric defined by RouteLink.nc
            domain_file            = hybrid_params.get("diffusive_domain",   None)
            topobathy_file         = hybrid_params.get("topobathy_domain",   None)
            
            # file path parameters of refactored hydrofabric for diffusive wave channel routing 
            refactored_domain_file    = hybrid_params.get("refactored_domain", None)
            refactored_topobathy_file = hybrid_params.get("refactored_topobathy_domain", None)
            #-------------------------------------------------------------------------
            # for non-refactored hydofabric defined by RouteLink.nc
            # TODO: By default, make diffusive available for both non-refactored and refactored hydrofabric for now. Place a switch in the future. 
            if run_hybrid and domain_file:
                
                LOG.info('reading diffusive domain extent for MC/Diffusive hybrid simulation')
                
                # read diffusive domain dictionary from yaml or json
                diffusive_domain = nhd_io.read_diffusive_domain(domain_file)           
                
                if use_topobathy and topobathy_file:
                    
                    LOG.debug('Natural cross section data on original hydrofabric are provided.')
                    
                    # read topobathy domain netcdf file, set index to 'comid'
                    # TODO: replace 'link' with a user-specified indexing variable name.
                    # ... if for whatever reason there is not a `link` variable in the 
                    # ... dataframe returned from read_netcdf, then the code would break here.
                    self._topobathy_df = (nhd_io.read_netcdf(topobathy_file).set_index('link'))
                    
                    # TODO: Request GID make comID variable an integer in their product, so
                    # we do not need to change variable types, here.
                    self._topobathy_df.index = self._topobathy_df.index.astype(int)
                    
                else:
                    self._topobathy_df = pd.DataFrame()
                    LOG.debug('No natural cross section topobathy data provided. Hybrid simualtion will run on compound trapezoidal geometry.')
                
                # initialize a dictionary to hold network data for each of the diffusive domains
                self._diffusive_network_data = {}
            
            else:
                diffusive_domain       = None
                self._diffusive_network_data = None
                self._topobathy_df     = pd.DataFrame()
                LOG.info('No diffusive domain file specified in configuration file. This is an MC-only simulation')
            self._unrefactored_topobathy_df  = pd.DataFrame()    
            #-------------------------------------------------------------------------
            # for refactored hydofabric 
            if run_hybrid and run_refactored and refactored_domain_file:
                
                LOG.info('reading refactored diffusive domain extent for MC/Diffusive hybrid simulation')
                
                # read diffusive domain dictionary from yaml or json
                self._refactored_diffusive_domain = nhd_io.read_diffusive_domain(refactored_domain_file)           
                
                if use_topobathy and refactored_topobathy_file:
                    
                    LOG.debug('Natural cross section data of refactored hydrofabric are provided.')
                    
                    # read topobathy domain netcdf file, set index to 'comid'
                    # TODO: replace 'link' with a user-specified indexing variable name.
                    # ... if for whatever reason there is not a `link` variable in the 
                    # ... dataframe returned from read_netcdf, then the code would break here.
                    self._topobathy_df = (nhd_io.read_netcdf(refactored_topobathy_file).set_index('link'))

                    # unrefactored_topobaty_data is passed to diffusive kernel to provide thalweg elevation of unrefactored topobathy 
                    # for crosswalking water elevations between non-refactored and refactored hydrofabrics. 
                    self._unrefactored_topobathy_df       = (nhd_io.read_netcdf(topobathy_file).set_index('link'))
                    self._unrefactored_topobathy_df.index = self._unrefactored_topobathy_df.index.astype(int)
                    
                else:
                    self._topobathy_df = pd.DataFrame()
                    LOG.debug('No natural cross section topobathy data of refactored hydrofabric provided. Hybrid simualtion will run on compound trapezoidal geometry.')
                
                # initialize a dictionary to hold network data for each of the diffusive domains
                refactored_diffusive_network_data = {}
            
            else:
                self._refactored_diffusive_domain       = None
                refactored_diffusive_network_data = None
                self._refactored_reaches                = {}
                LOG.info('No refactored diffusive domain file specified in configuration file. This is an MC-only simulation')     
    
        else:
            diffusive_domain                  = None
            self._diffusive_network_data            = None
            self._topobathy_df                = pd.DataFrame()
            self._unrefactored_topobathy_df         = pd.DataFrame() 
            self._refactored_diffusive_domain       = None
            refactored_diffusive_network_data = None   
            self._refactored_reaches                = {}
            LOG.info('No hybrid parameters specified in configuration file. This is an MC-only simulation')
        
        #============================================================================
        # build diffusive domain data and edit MC domain data for hybrid simulation
        
        #
        if diffusive_domain:
            rconn_diff0 = reverse_network(self._connections)
            self._refactored_reaches = {}
            
            for tw in diffusive_domain:
                mainstem_segs = diffusive_domain[tw]['links']
                # we want mainstem_segs start at a mainstem link right after the upstream boundary mainstem link, which is
                # in turn not under any waterbody. This boundary mainstem link should be turned into a tributary segment.
                upstream_boundary_mainstem_link = diffusive_domain[tw]['upstream_boundary_link_mainstem']         
                if upstream_boundary_mainstem_link[0] in mainstem_segs:
                    mainstem_segs.remove(upstream_boundary_mainstem_link[0])
                
                # ===== build diffusive network data objects ==== 
                self._diffusive_network_data[tw] = {}

                # add diffusive domain segments
                self._diffusive_network_data[tw]['mainstem_segs'] =  mainstem_segs

                # diffusive domain tributary segments
                trib_segs = []
                
                for seg in mainstem_segs:
                    us_list = rconn_diff0[seg]
                    for u in us_list:
                        if u not in mainstem_segs:
                            trib_segs.append(u) 

                self._diffusive_network_data[tw]['tributary_segments'] = trib_segs
                # diffusive domain connections object
                self._diffusive_network_data[tw]['connections'] = {k: self._connections[k] for k in (mainstem_segs + trib_segs)}       

                # diffusive domain reaches and upstream connections. 
                # break network at tributary segments
                _, reaches, rconn_diff = organize_independent_networks(
                    self._diffusive_network_data[tw]['connections'],
                    set(trib_segs),
                    set(),
                )
                
                self._diffusive_network_data[tw]['rconn'] = rconn_diff
                self._diffusive_network_data[tw]['reaches'] = reaches[tw]

                # RouteLink parameters
                self._diffusive_network_data[tw]['param_df'] = self._dataframe.filter(
                    (mainstem_segs + trib_segs),
                    axis = 0,
                )
                self._diffusive_network_data[tw]['upstream_boundary_link'] = upstream_boundary_mainstem_link

                if self._refactored_diffusive_domain: 
                    diffusive_parameters = {'geo_file_path': refactored_topobathy_file}
                    refactored_connections = build_refac_connections(diffusive_parameters)

                    # list of stream segments of a single refactored diffusive domain 
                    refac_tw = self._refactored_diffusive_domain[tw]['refac_tw']
                    rlinks_tw = self._refactored_diffusive_domain[tw]['rlinks']
                    refactored_connections_tw = {}   

                    # Subset a connection dictionary (upstream segment as key : downstream segments as values) from refactored_connections
                    # for a single refactored diffusive domain defined by a current tw. 
                    for k in rlinks_tw:
                        if k in refactored_connections.keys() and k != refac_tw:
                            refactored_connections_tw[k] = refactored_connections[k]

                    refactored_diffusive_network_data[refac_tw] = {}                
                    refactored_diffusive_network_data[refac_tw]['tributary_segments'] = trib_segs
                    refactored_diffusive_network_data[refac_tw]['connections'] = refactored_connections_tw                 

                    for k in trib_segs:
                        refactored_diffusive_network_data[refac_tw]['connections'][k]= [self._refactored_diffusive_domain[tw]['incoming_tribs'][k]]

                    # diffusive domain reaches and upstream connections. 
                    # break network at tributary segments
                    _, refactored_reaches_batch, refactored_conn_diff = organize_independent_networks(
                                                                refactored_diffusive_network_data[refac_tw]['connections'],
                                                                set(trib_segs),
                                                                set(),
                                                                )

                    self._refactored_reaches[refac_tw] = refactored_reaches_batch[refac_tw]
                    refactored_diffusive_network_data[refac_tw]['mainstem_segs'] = self._refactored_diffusive_domain[tw]['rlinks']
                    refactored_diffusive_network_data[refac_tw]['upstream_boundary_link'] = self._diffusive_network_data[tw]['upstream_boundary_link'] 
                else:
                    self._refactored_reaches={}

                # ==== remove diffusive domain segs from MC domain ====        
                # drop indices from param_df
                self._dataframe = self._dataframe.drop(mainstem_segs)

                # remove keys from connections dictionary
                for s in mainstem_segs:
                    self._connections.pop(s)

                # update downstream connections of trib segs
                for us in trib_segs:
                    self._connections[us] = []

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
        - break_network_at_waterbodies (bool): If True, waterbody initial states will
                                            be appended to the waterbody parameter
                                            dataframe. If False, waterbodies will
                                            not be simulated and the waterbody
                                            parameter datataframe wil not be changed
        - restart_parameters           (dict): User-input simulation restart 
                                            parameters
        - segment_index        (Pandas Index): All segment IDs in the simulation 
                                            doamin
        - waterbodies_df   (Pandas DataFrame): Waterbody parameters
        
        Returns
        -------
        - waterbodies_df (Pandas DataFrame): Waterbody parameters with initial
                                            states (outflow and pool elevation)
        - q0             (Pandas DataFrame): Initial flow and depth states for each
                                            segment in the model domain
        - t0                     (datetime): Datetime of the model initialization
        
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
        #----------------------------------------------------------------------------    
        start_time = time.time()
        LOG.info("setting channel initial states ...")

        # if lite restart file is provided, the read channel initial states from it
        if restart_parameters.get("lite_channel_restart_file", None):
            # FIXME: Change it for hyfeature!
            self._q0, self._t0 = nhd_io.read_lite_restart(
                restart_parameters['lite_channel_restart_file']
            )
            t0_str = None

        # when a restart file for hyfeature is provied, then read initial states from it.
        elif restart_parameters.get("hyfeature_channel_restart_file", None):        
            self._q0 = build_channel_initial_state(restart_parameters, self.segment_index)        
            channel_initial_states_file = restart_parameters["hyfeature_channel_restart_file"]
            df     = pd.read_csv(channel_initial_states_file)
            t0_str = pd.to_datetime(df.columns[1]).strftime("%Y-%m-%d_%H:%M:%S")
            self._t0     = datetime.strptime(t0_str,"%Y-%m-%d_%H:%M:%S")

        # build initial states from user-provided restart parameters
        else:
            # FIXME: Change it for hyfeature!
            self._q0 = build_channel_initial_state(restart_parameters, self.segment_index)       

            # get initialization time from restart file
            if restart_parameters.get("wrf_hydro_channel_restart_file", None):
                channel_initial_states_file = restart_parameters[
                    "wrf_hydro_channel_restart_file"
                ]
                t0_str = nhd_io.get_param_str(
                    channel_initial_states_file, 
                    "Restart_Time"
                )
            else:
                t0_str = "2015-08-16_00:00:00"

            # convert timestamp from string to datetime
            self._t0 = datetime.strptime(t0_str, "%Y-%m-%d_%H:%M:%S")

        # get initial time from user inputs
        if restart_parameters.get("start_datetime", None):
            t0_str = restart_parameters.get("start_datetime")
            
            def _try_parsing_date(text):
                for fmt in (
                    "%Y-%m-%d_%H:%M", 
                    "%Y-%m-%d_%H:%M:%S", 
                    "%Y-%m-%d %H:%M", 
                    "%Y-%m-%d %H:%M:%S", 
                    "%Y/%m/%d %H:%M", 
                    "%Y/%m/%d %H:%M:%S"
                ):
                    try:
                        return datetime.strptime(text, fmt)
                    except ValueError:
                        pass
                LOG.error('No valid date format found for start_datetime input. Please use format YYYY-MM-DD_HH:MM')
                quit()
                
            self._t0 = _try_parsing_date(t0_str)
        else:
            if t0_str == "2015-08-16_00:00:00":
                LOG.info('No user-input start_datetime and no restart file, start time arbitrarily 2015-08-16_00:00:00')
            else:
                LOG.info('No user-specified start_datetime, continuing with start time from restart file: %s', t0_str)

        LOG.debug(
            "channel initial states complete in %s seconds."\
            % (time.time() - start_time)
        )

    def build_qlateral_array(self, run, cpu_pool, supernetwork_parameters,):

        start_time = time.time()
        LOG.info("Creating a DataFrame of lateral inflow forcings ...")
        
        # TODO: set default/optional arguments
        qts_subdivisions = run.get("qts_subdivisions", 1)
        nts = run.get("nts", 1)
        qlat_input_folder = run.get("qlat_input_folder", None)
        qlat_input_file = run.get("qlat_input_file", None)

        geo_file_type = supernetwork_parameters.get('geo_file_type')

        if qlat_input_folder:
            qlat_input_folder = pathlib.Path(qlat_input_folder)
            if "qlat_files" in run:
                qlat_files = run.get("qlat_files")
                qlat_files = [qlat_input_folder.joinpath(f) for f in qlat_files]
            elif "qlat_file_pattern_filter" in run:
                qlat_file_pattern_filter = run.get(
                    "qlat_file_pattern_filter", "*CHRT_OUT*"
                )
                qlat_files = sorted(qlat_input_folder.glob(qlat_file_pattern_filter))

            qlat_file_index_col = run.get(
                "qlat_file_index_col", "feature_id"
            )

            if geo_file_type=='NHDNetwork':
                # Parallel reading of qlateral data from CHRTOUT
                with Parallel(n_jobs=cpu_pool) as parallel:
                    jobs = []
                    for f in qlat_files:
                        jobs.append(
                            delayed(nhd_io.get_ql_from_chrtout)
                            #(f, qlat_file_value_col, gw_bucket_col, terrain_ro_col)
                            #delayed(nhd_io.get_ql_from_csv)
                            (f)                    
                        )
                    ql_list = parallel(jobs)

                # get feature_id from a single CHRTOUT file
                with netCDF4.Dataset(qlat_files[0]) as ds:
                    idx = ds.variables[qlat_file_index_col][:].filled()

                # package data into a DataFrame
                qlats_df = pd.DataFrame(
                    np.stack(ql_list).T,
                    index = idx,
                    columns = range(len(qlat_files))
                )

                qlats_df = qlats_df[qlats_df.index.isin(self.segment_index)]
            elif geo_file_type=='HYFeaturesNetwork':
                dfs=[]
                for f in qlat_files:
                    df = read_file(f).set_index(['feature_id']) 
                    dfs.append(df)
                
                # lateral flows [m^3/s] are stored at NEXUS points with NEXUS ids
                nexuses_lateralflows_df = pd.concat(dfs, axis=1)  
                
                # Take flowpath ids entering NEXUS and replace NEXUS ids by the upstream flowpath ids
                qlats_df = pd.concat( (nexuses_lateralflows_df.loc[int(k)].rename(v)
                                    for k,v in self.downstream_flowpath_dict.items() ), axis=1
                                    ).T 
                qlats_df.columns=range(len(qlat_files))
                qlats_df = qlats_df[qlats_df.index.isin(self.segment_index)]

                # The segment_index has the full network set of segments/flowpaths. 
                # Whereas the set of flowpaths that are downstream of nexuses is a 
                # subset of the segment_index. Therefore, all of the segments/flowpaths
                # that are not accounted for in the set of flowpaths downstream of
                # nexuses need to be added to the qlateral dataframe and padded with
                # zeros.
                all_df = pd.DataFrame( np.zeros( (len(self.segment_index), len(qlats_df.columns)) ), index=self.segment_index,
                    columns=qlats_df.columns )
                all_df.loc[ qlats_df.index ] = qlats_df
                qlats_df = all_df.sort_index()

        elif qlat_input_file:
            qlats_df = nhd_io.get_ql_from_csv(qlat_input_file)
        else:
            qlat_const = run.get("qlat_const", 0)
            qlats_df = pd.DataFrame(
                qlat_const,
                index=self.segment_index,
                columns=range(nts // qts_subdivisions),
                dtype="float32",
            )

        # TODO: Make a more sophisticated date-based filter
        max_col = 1 + nts // qts_subdivisions
        if len(qlats_df.columns) > max_col:
            qlats_df.drop(qlats_df.columns[max_col:], axis=1, inplace=True)

        if not self.segment_index.empty:
            qlats_df = qlats_df[qlats_df.index.isin(self.segment_index)]
        
        LOG.debug(
            "lateral inflow DataFrame creation complete in %s seconds." \
                % (time.time() - start_time)
                )

        self._qlateral = qlats_df



def read_file(file_name):
    extension = file_name.suffix
    if extension=='.csv':
        df = pd.read_csv(file_name)
    elif extension=='.parquet':
        df = pq.read_table(file_name).to_pandas().reset_index()
        df.index.name = None
    
    return df