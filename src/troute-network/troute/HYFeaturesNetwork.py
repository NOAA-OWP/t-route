from .AbstractNetwork import AbstractNetwork
import pathlib
import json
import pandas as pd
import numpy as np
import time
import re
import troute.nhd_io as nhd_io #FIXME
from itertools import chain
import geopandas as gpd
from pathlib import Path
import math
import troute.hyfeature_preprocess as hyfeature_prep
from datetime import datetime, timedelta

__verbose__ = False
__showtiming__ = False

def read_qlats(forcing_parameters, segment_index, nexus_to_downstream_flowpath_dict):
    # STEP 5: Read (or set) QLateral Inputs
    if __showtiming__:
        start_time = time.time()
    if __verbose__:
        print("creating qlateral array ...")
    qts_subdivisions = forcing_parameters.get("qts_subdivisions", 1)
    nts = forcing_parameters.get("nts", 1)
    nexus_input_folder = forcing_parameters.get("nexus_input_folder", None)
    nexus_input_folder = pathlib.Path(nexus_input_folder)
    #rt0 = time.time() #FIXME showtiming flag
    if not "nexus_file_pattern_filter" in forcing_parameters:
        raise( RuntimeError("No value for nexus file pattern in config" ) )
    
    nexus_file_pattern_filter = forcing_parameters.get(
        "nexus_file_pattern_filter", "nex-*"
    )
    nexus_files = nexus_input_folder.glob(nexus_file_pattern_filter)
    #TODO Find a way to test that we found some files
    #without consuming the generator...Otherwise, if nexus_files
    #is empty, the following concat raises a ValueError which
    #It may be sufficient to catch this exception and warn that
    #there may not be any pattern matching files in the dir

    nexus_files_list = list(nexus_files)

    if len(nexus_files_list) == 0:
        raise ValueError('No nexus input files found. Recommend checking \
        nexus_input_folder path in YAML configuration.')

    #pd.concat((pd.read_csv(f, index_col=0, usecols=[1,2], header=None, engine='c').rename(columns={2:id_regex.match(f).group(1)}) for f in all_files[0:2]), axis=1).T
    id_regex = re.compile(r".*nex-(\d+)_.*.csv")
    nexuses_flows_df = pd.concat(
            #Read the nexus csv file
            (pd.read_csv(f, index_col=0, usecols=[1,2], header=None, engine='c', skipinitialspace=True, parse_dates=True).rename(
                #Rename the flow column to the id of the nexus
                columns={2:int(id_regex.match(f.name).group(1))})
            for f in nexus_files_list #Build the generator for each required file
            ),  axis=1).T #Have now concatenated a single df (along axis 1).  Transpose it.
    missing = nexuses_flows_df[ nexuses_flows_df.isna().any(axis=1) ]
    if  not missing.empty:
        raise ValueError("The following nexus inputs are incomplete: "+str(missing.index))
    rt1 = time.time()
    #print("Time to build nexus_flows_df: {} seconds".format(rt1-rt0))

    qlat_df = pd.concat( (nexuses_flows_df.loc[int(k)].rename(index={int(k):v})
        for k,v in nexus_to_downstream_flowpath_dict.items() ), axis=1
        ).T
    #qlat_df = pd.concat( (nexuses_flows_df.loc[int(k)].rename(v)
    #    for k,v in nexus_to_downstream_flowpath_dict.items() ), axis=1
    #    ).T
    
    # The segment_index has the full network set of segments/flowpaths. 
    # Whereas the set of flowpaths that are downstream of nexuses is a 
    # subset of the segment_index. Therefore, all of the segments/flowpaths
    # that are not accounted for in the set of flowpaths downstream of
    # nexuses need to be added to the qlateral dataframe and padded with
    # zeros.
    all_df = pd.DataFrame( np.zeros( (len(segment_index), len(qlat_df.columns)) ), index=segment_index,
            columns=qlat_df.columns )
    all_df.loc[ qlat_df.index ] = qlat_df
    qlat_df = all_df.sort_index()

    # Set new nts based upon total nexus inputs
    nts = (qlat_df.shape[1]) * qts_subdivisions
    max_col = 1 + nts // qts_subdivisions
    
    #dt      = 300 # [sec]
    #dt_qlat = 3600 # [sec]
    #nts     = 24 # steps
    #max_col = math.ceil(nts*dt/dt_qlat)

    if len(qlat_df.columns) > max_col:
        qlat_df.drop(qlat_df.columns[max_col:], axis=1, inplace=True)

    if __verbose__:
        print("qlateral array complete")
    if __showtiming__:
        print("... in %s seconds." % (time.time() - start_time))

    return qlat_df

class HYFeaturesNetwork(AbstractNetwork):
    """
    
    """
    __slots__ = ["_flowpath_dict", 
                 "segment_index", 
                 "waterbody_type_specified",
                 "diffusive_network_data", 
                 "topobathy_df", 
                 "refactored_diffusive_domain",
                 "refactored_reaches", 
                 "unrefactored_topobathy_df"]
    def __init__(self, 
                 supernetwork_parameters, 
                 waterbody_parameters=None, 
                 restart_parameters=None, 
                 forcing_parameters=None, 
                 verbose=False, 
                 showtiming=False):
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
        # Load Geo File
        #------------------------------------------------
        (self._dataframe,
         self._flowpath_dict,
         self._waterbody_df,
         self._waterbody_types_df,
        ) = hyfeature_prep.read_geo_file(
            supernetwork_parameters,
            waterbody_parameters,
        )






        #------------------------------------------------
        # Preprocess network attributes
        #------------------------------------------------        
        (self._dataframe,            
         self._flowpath_dict, 
         self._waterbody_types_df,
         self._waterbody_df,
         self.waterbody_type_specified,
         cols,
         terminal_code,
         break_points,
        ) = hyfeature_prep.build_hyfeature_network(
            supernetwork_parameters,
            waterbody_parameters
        )
        
        # called to mainly initialize _waterbody_connections, _connections, _independent_networks,
        # _reverse_network, _reaches_by_tw
        super().__init__(cols, terminal_code, break_points)    
        
        if __verbose__:
            print("supernetwork connections set complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))   
        
   
            
        # list of all segments in the domain (MC + diffusive)
        self.segment_index = self._dataframe.index
        #if self.diffusive_network_data:
        #    for tw in self.diffusive_network_data:
        #        self.segment_index = self.segment_index.append(
        #            pd.Index(self.diffusive_network_data[tw]['mainstem_segs'])
        #        )     
            
        #------------------------------------------------
        #  Handle Channel Initial States
        #------------------------------------------------ 
        if __verbose__:
            print("setting waterbody and channel initial states ...")
        if __showtiming__:
            start_time = time.time()

        (#self._waterbody_df,
         self._q0,
         self._t0,) = hyfeature_prep.hyfeature_initial_warmstate_preprocess(
            #break_network_at_waterbodies,
            restart_parameters,
            #data_assimilation_parameters,
            self.segment_index,
            #self._waterbody_df,
            #self.link_lake_crosswalk,
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
        ) = hyfeature_prep.hyfeature_forcing(
            run, 
            forcing_parameters, 
            hybrid_parameters,
            self._flowpath_dict,
            self.segment_index,
            cpu_pool,
            self._t0,             
            self._coastal_boundary_depth_df,
        )

        #Mask out all non-simulated waterbodies
        self._dataframe['waterbody'] = self.waterbody_null
        
        #This also remaps the initial NHDComID identity to the HY_Features Waterbody ID for the reservoir...
        self._dataframe.loc[self._waterbody_df.index, 'waterbody'] = self._waterbody_df.index.name
        
        #FIXME should waterbody_df and param_df overlap IDS?  Doesn't seem like it should...
        #self._dataframe.drop(self._waterbody_df.index, axis=0, inplace=True)
        #For now, doing it in property waterbody_connections...
        # Remove duplicate rows...but its not really duplicate...
        #self._dataframe = (
        #    self._dataframe.reset_index()
        #    .drop_duplicates(subset="key")
        #    .set_index("key")
        #)
        self._dataframe.sort_index(inplace=True)
        self._waterbody_df.sort_index(inplace=True)

    def extract_waterbody_connections(rows, target_col, waterbody_null=-9999):
        """Extract waterbody mapping from dataframe.
        TODO deprecate in favor of waterbody_connections property"""
        return (
            rows.loc[rows[target_col] != waterbody_null, target_col].astype("int").to_dict()
        )
 
    
    def create_routing_network(self, 
                               conn, 
                               param_df, 
                               wbody_conn, gages, 
                               preprocessing_parameters, 
                               compute_parameters,
                               waterbody_parameters
    ): 

        #--------------------------------------------------------------------------
        # Creation of routing network data objects. Logical ordering of lower-level
        # function calls that build individual network data objects.
        #--------------------------------------------------------------------------        
        (self._independent_networks,
         self._reaches_by_tw,
         self._reverse_network,
         self.diffusive_network_data,
         self.topobathy_df,
         self.refactored_diffusive_domain,
         self.refactored_reaches,
         self.unrefactored_topobathy_df
        ) = hyfeature_prep.hyfeature_hybrid_routing_preprocess(
            conn,
            param_df,
            wbody_conn,
            gages,
            preprocessing_parameters,
            compute_parameters,
            waterbody_parameters, 
        )
        return (self._independent_networks,
                self._reaches_by_tw,
                self._reverse_network,
                self.diffusive_network_data,
                self.topobathy_df,
                self.refactored_diffusive_domain,
                self.refactored_reaches,
                self.unrefactored_topobathy_df)

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
    def downstream_flowpath_dict(self):
        return self._flowpath_dict

    @property
    def waterbody_connections(self):
        """
            A dictionary where the keys are the reach/segment id, and the
            value is the id to look up waterbody parameters
        """
        if( not self._waterbody_connections ):
            #Funny story, NaN != NaN is true!!!!
            #Drop the nan, then check for waterbody_null just in case
            #waterbody_null happens to be NaN
            #FIXME this drops ALL nan, not just `waterbody`
            #waterbody_segments = self._dataframe.dropna().loc[
            #    self._dataframe["waterbody"] != self.waterbody_null, "waterbody"
            #]
            #waterbody_segments = waterbody_segments.loc[self.waterbody_dataframe.index]
            #self._waterbody_connections = waterbody_segments.index\
            #    .to_series(name = waterbody_segments.name)\
            #    .astype("int")\
            #    .to_dict()
            #If we identify as a waterbody, drop from the main dataframe
            #Ugh, but this drops everything that that MIGHT be a "lake"
            #without knowing if it was defined as a lake in the lake params
            #so this should just drop the waterbody_df index, not these segments...
            #In fact, these waterbody connections should probably be entirely reworked
            #with that in mind...
            self._waterbody_connections = self._waterbody_df.index.to_series(name = self._waterbody_df.index.name).astype("int").to_dict()
            #FIXME seems way more appropriate to do this in the constructor so the property doesn't side effect
            #the param df..., but then it breaks down the connection property...so for now, leave it here and fix later
            self._dataframe.drop(self._waterbody_df.index, axis=0, inplace=True)
        return self._waterbody_connections

    @property
    def gages(self):
        """
        FIXME
        """
        if self._gages is None and "gages" in self._dataframe.columns:
            self._gages = nhd_io.build_filtered_gage_df(self._dataframe[["gages"]])
        else:
            self._gages = {}
        return self._gages
    
    @property
    def waterbody_null(self):
        return np.nan #pd.NA

