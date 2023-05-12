from .AbstractNetwork import AbstractNetwork
import pandas as pd
import numpy as np
import geopandas as gpd
import time
import json
from pathlib import Path
import pyarrow.parquet as pq

import troute.nhd_io as nhd_io #FIXME
from troute.nhd_network import reverse_dict, extract_connections

__verbose__ = False
__showtiming__ = False

def read_geopkg(file_path):
    flowpaths = gpd.read_file(file_path, layer="flowpaths")
    attributes = gpd.read_file(file_path, layer="flowpath_attributes").drop('geometry', axis=1)
    #merge all relevant data into a single dataframe
    flowpaths = pd.merge(flowpaths, attributes, on='id')

    return flowpaths

def read_json(file_path, edge_list):
    dfs = []
    with open(edge_list) as edge_file:
        edge_data = json.load(edge_file)
        edge_map = {}
        for id_dict in edge_data:
            edge_map[ id_dict['id'] ] = id_dict['toid']
        with open(file_path) as data_file:
            json_data = json.load(data_file)  
            for key_wb, value_params in json_data.items():
                df = pd.json_normalize(value_params)
                df['id'] = key_wb
                df['toid'] = edge_map[key_wb]
                dfs.append(df)
        df_main = pd.concat(dfs, ignore_index=True)

    return df_main

def numeric_id(flowpath):
    id = flowpath['id'].split('-')[-1]
    toid = flowpath['toid'].split('-')[-1]
    flowpath['id'] = int(id)
    flowpath['toid'] = int(toid)
    return flowpath

def read_ngen_waterbody_df(parm_file, lake_index_field="wb-id", lake_id_mask=None):
    """
    Reads .gpkg or lake.json file and prepares a dataframe, filtered
    to the relevant reservoirs, to provide the parameters
    for level-pool reservoir computation.
    """
    def node_key_func(x):
        return int( x.split('-')[-1] )
    if Path(parm_file).suffix=='.gpkg':
        df = gpd.read_file(parm_file, layer="lake_attributes").set_index('id')
    elif Path(parm_file).suffix=='.json':
        df = pd.read_json(parm_file, orient="index")

    df.index = df.index.map(node_key_func)
    df.index.name = lake_index_field

    if lake_id_mask:
        df = df.loc[lake_id_mask]
    return df

def read_ngen_waterbody_type_df(parm_file, lake_index_field="wb-id", lake_id_mask=None):
    """
    """
    #FIXME: this function is likely not correct. Unclear how we will get 
    # reservoir type from the gpkg files. Information should be in 'crosswalk'
    # layer, but as of now (Nov 22, 2022) there doesn't seem to be a differentiation
    # between USGS reservoirs, USACE reservoirs, or RFC reservoirs...
    def node_key_func(x):
        return int( x.split('-')[-1] )
    
    if Path(parm_file).suffix=='.gpkg':
        df = gpd.read_file(parm_file, layer="crosswalk").set_index('id')
    elif Path(parm_file).suffix=='.json':
        df = pd.read_json(parm_file, orient="index")

    df.index = df.index.map(node_key_func)
    df.index.name = lake_index_field
    if lake_id_mask:
        df = df.loc[lake_id_mask]
        
    return df


class HYFeaturesNetwork(AbstractNetwork):
    """
    
    """
    __slots__ = ["_upstream_terminal"]

    def __init__(self, 
                 supernetwork_parameters, 
                 waterbody_parameters,
                 data_assimilation_parameters,
                 restart_parameters, 
                 compute_parameters,
                 forcing_parameters,
                 hybrid_parameters, 
                 verbose=False, 
                 showtiming=False,
                 from_files=True,
                 value_dict={},
                 segment_attributes=[],
                 waterbody_attributes=[]):
        """
        
        """
        self.supernetwork_parameters = supernetwork_parameters
        self.waterbody_parameters = waterbody_parameters
        self.data_assimilation_parameters = data_assimilation_parameters
        self.restart_parameters = restart_parameters
        self.compute_parameters = compute_parameters
        self.forcing_parameters = forcing_parameters
        self.hybrid_parameters = hybrid_parameters
        self.verbose = verbose
        self.showtiming = showtiming

        if self.verbose:
            print("creating supernetwork connections set")
        if self.showtiming:
            start_time = time.time()
        
        #------------------------------------------------
        # Load Geo File
        #------------------------------------------------
        if from_files:
            self.read_geo_file()
        else:
            self.load_bmi_data(value_dict, segment_attributes, waterbody_attributes)
        
        #TODO Update for waterbodies and DA specific to HYFeatures...
        self._waterbody_connections = {}
        self._waterbody_type_specified = None
        self._gages = None
        self._link_lake_crosswalk = None


        if self.verbose:
            print("supernetwork connections set complete")
        if self.showtiming:
            print("... in %s seconds." % (time.time() - start_time))
            

        super().__init__()   
            
        # Create empty dataframe for coastal_boundary_depth_df. This way we can check if
        # it exists, and only read in SCHISM data during 'assemble_forcings' if it doesn't
        self._coastal_boundary_depth_df = pd.DataFrame()

    def extract_waterbody_connections(rows, target_col, waterbody_null=-9999):
        """Extract waterbody mapping from dataframe.
        TODO deprecate in favor of waterbody_connections property"""
        return (
            rows.loc[rows[target_col] != waterbody_null, target_col].astype("int").to_dict()
        )

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
    
    def read_geo_file(self,):
        
        geo_file_path = self.supernetwork_parameters["geo_file_path"]
        
        file_type = Path(geo_file_path).suffix
        if(  file_type == '.gpkg' ):        
            self._dataframe = read_geopkg(geo_file_path)
        elif( file_type == '.json') :
            edge_list = self.supernetwork_parameters['flowpath_edge_list']
            self._dataframe = read_json(geo_file_path, edge_list) 
        else:
            raise RuntimeError("Unsupported file type: {}".format(file_type))
        
        # Don't need the string prefix anymore, drop it
        mask = ~ self.dataframe['toid'].str.startswith("tnex") 
        self._dataframe = self.dataframe.apply(numeric_id, axis=1)
        
        # make the flowpath linkage, ignore the terminal nexus
        self._flowpath_dict = dict(zip(self.dataframe.loc[mask].toid, self.dataframe.loc[mask].id))
        
        # **********  need to be included in flowpath_attributes  *************
        self._dataframe['alt'] = 1.0 #FIXME get the right value for this... 

        cols = self.supernetwork_parameters.get('columns',None)
        
        if cols:
            self._dataframe = self.dataframe[list(cols.values())]
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
            self._dataframe = self.dataframe.rename(columns=reverse_dict(cols))
            self._dataframe.set_index("key", inplace=True)
            self._dataframe = self.dataframe.sort_index()

        # numeric code used to indicate network terminal segments
        terminal_code = self.supernetwork_parameters.get("terminal_code", 0)

        # There can be an externally determined terminal code -- that's this first value
        self._terminal_codes = set()
        self._terminal_codes.add(terminal_code)
        # ... but there may also be off-domain nodes that are not explicitly identified
        # but which are terminal (i.e., off-domain) as a result of a mask or some other
        # an interior domain truncation that results in a
        # otherwise valid node value being pointed to, but which is masked out or
        # being intentionally separated into another domain.
        self._terminal_codes = self.terminal_codes | set(
            self.dataframe[~self.dataframe["downstream"].isin(self.dataframe.index)]["downstream"].values
        )

        #This is NEARLY redundant to the self.terminal_codes property, but in this case
        #we actually need the mapping of what is upstream of that terminal node as well.
        #we also only want terminals that actually exist based on definition, not user input
        terminal_mask = ~self._dataframe["downstream"].isin(self._dataframe.index)
        terminal = self._dataframe.loc[ terminal_mask ]["downstream"]
        self._upstream_terminal = dict()
        for key, value in terminal.items():
            self._upstream_terminal.setdefault(value, set()).add(key)

        # build connections dictionary
        self._connections = extract_connections(
            self.dataframe, "downstream", terminal_codes=self.terminal_codes
        )

        #Load waterbody/reservoir info
        if self.waterbody_parameters:
            levelpool_params = self.waterbody_parameters.get('level_pool', None)
            if not levelpool_params:
                # FIXME should not be a hard requirement
                raise(RuntimeError("No supplied levelpool parameters in routing config"))
                
            lake_id = levelpool_params.get("level_pool_waterbody_id", "wb-id")
            try:
                self._waterbody_df = read_ngen_waterbody_df(
                            levelpool_params["level_pool_waterbody_parameter_file_path"],
                            lake_id,
                            )
                    
                # Remove duplicate lake_ids and rows
                self._waterbody_df = (
                                self.waterbody_dataframe.reset_index()
                                .drop_duplicates(subset=lake_id)
                                .set_index(lake_id)
                                .sort_index()
                                )
            except ValueError:
                self._waterbody_df = pd.DataFrame()

            try:
                self._waterbody_types_df = read_ngen_waterbody_type_df(
                                        levelpool_params["reservoir_parameter_file"],
                                        lake_id,
                                        #self.waterbody_connections.values(),
                                        )
                # Remove duplicate lake_ids and rows
                self._waterbody_types_df =(
                                    self.waterbody_types_dataframe.reset_index()
                                    .drop_duplicates(subset=lake_id)
                                    .set_index(lake_id)
                                    .sort_index()
                                    )

            except ValueError:
                #FIXME any reservoir operations requires some type
                #So make this default to 1 (levelpool)
                self._waterbody_types_df = pd.DataFrame(index=self.waterbody_dataframe.index)
                self._waterbody_types_df['reservoir_type'] = 1
        else:
            self._waterbody_df = pd.DataFrame()
            self._waterbody_types_df = pd.DataFrame()
    
    def load_bmi_data(self, value_dict, segment_attributes, waterbody_attributes):
        
        self._dataframe = pd.DataFrame(data=None, columns=segment_attributes)
        self._waterbody_df = pd.DataFrame(data=None, columns=waterbody_attributes)

        for var in segment_attributes:
            self._dataframe[var] = value_dict[var]
        
        # make the flowpath linkage, ignore the terminal nexus
        self._flowpath_dict = dict(zip(self.dataframe.segment_toid, self.dataframe.segment_id))
        
        # **********  need to be included in flowpath_attributes  *************
        self._dataframe['alt'] = 1.0 #FIXME get the right value for this... 

        self._dataframe = self.dataframe.rename(columns={'segment_id': 'key',
                                                         'segment_toid': 'downstream'})
        self._dataframe.set_index("key", inplace=True)
        self._dataframe = self.dataframe.sort_index()

        # numeric code used to indicate network terminal segments
        terminal_code = self.supernetwork_parameters.get("terminal_code", 0)

        # There can be an externally determined terminal code -- that's this first value
        self._terminal_codes = set()
        self._terminal_codes.add(terminal_code)
        # ... but there may also be off-domain nodes that are not explicitly identified
        # but which are terminal (i.e., off-domain) as a result of a mask or some other
        # an interior domain truncation that results in a
        # otherwise valid node value being pointed to, but which is masked out or
        # being intentionally separated into another domain.
        self._terminal_codes = self.terminal_codes | set(
            self.dataframe[~self.dataframe["downstream"].isin(self.dataframe.index)]["downstream"].values
        )

        #This is NEARLY redundant to the self.terminal_codes property, but in this case
        #we actually need the mapping of what is upstream of that terminal node as well.
        #we also only want terminals that actually exist based on definition, not user input
        terminal_mask = ~self._dataframe["downstream"].isin(self._dataframe.index)
        terminal = self._dataframe.loc[ terminal_mask ]["downstream"]
        self._upstream_terminal = dict()
        for key, value in terminal.items():
            self._upstream_terminal.setdefault(value, set()).add(key)

        # build connections dictionary
        self._connections = extract_connections(
            self.dataframe, "downstream", terminal_codes=self.terminal_codes
        )
        
        #Load waterbody/reservoir info
        if self.waterbody_parameters:
            for var in waterbody_attributes:
                self._waterbody_df[var] = value_dict[var]

            self._waterbody_df = self._waterbody_df.rename(columns={'waterbody_id': 'wb-id',
                                                                    'waterbody_toid': 'toid'})
            self._waterbody_df.set_index("wb-id", inplace=True)
            self._waterbody_df = self.waterbody_dataframe.sort_index()

            self._waterbody_types_df = pd.DataFrame(self.waterbody_dataframe['reservoir_type'])
            self._waterbody_df.drop('reservoir_type', axis=1, inplace=True)

    
    def build_qlateral_array(self, run,):
        
        # TODO: set default/optional arguments
        qts_subdivisions = run.get("qts_subdivisions", 1)
        nts = run.get("nts", 1)
        qlat_input_folder = run.get("qlat_input_folder", None)
        qlat_input_file = run.get("qlat_input_file", None)

        if qlat_input_folder:
            qlat_input_folder = Path(qlat_input_folder)
            if "qlat_files" in run:
                qlat_files = run.get("qlat_files")
                qlat_files = [qlat_input_folder.joinpath(f) for f in qlat_files]
            elif "qlat_file_pattern_filter" in run:
                qlat_file_pattern_filter = run.get(
                    "qlat_file_pattern_filter", "*CHRT_OUT*"
                )
                qlat_files = sorted(qlat_input_folder.glob(qlat_file_pattern_filter))

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

            '''
            #For a terminal nexus, we want to include the lateral flow from the catchment contributing to that nexus
            #one way to do that is to cheat and put that lateral flow at the upstream...this is probably the simplest way
            #right now.  The other is to create a virtual channel segment downstream to "route" i.e accumulate into
            #but it isn't clear right now how to do that with flow/velocity/depth requirements
            #find the terminal nodes
            for tnx, test_up in self._upstream_terminal.items():
                #first need to ensure there is an upstream location to dump to
                pdb.set_trace()
                for nex in test_up:
                    try:
                        #FIXME if multiple upstreams exist in this case then a choice is to be made as to which it goes into
                        #some cases the choice is easy cause the upstream doesn't exist, but in others, it may not be so simple
                        #in such cases where multiple valid upstream nexuses exist, perhaps the mainstem should be used?
                        pdb.set_trace()
                        qlats_df.loc[up] += nexuses_lateralflows_df.loc[tnx]
                        break #flow added, don't add it again!
                    except KeyError:
                        #this upstream doesn't actually exist on the network (maybe it is a headwater?)
                        #or perhaps the output file doesnt exist?  If this is the case, this isn't a good trap
                        #but for now, add the flow to a known good nexus upstream of the terminal
                        continue
                    #TODO what happens if can't put the qlat anywhere?  Right now this silently ignores the issue...
                qlats_df.drop(tnx, inplace=True)
            '''

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

        self._qlateral = qlats_df

def read_file(file_name):
    extension = file_name.suffix
    if extension=='.csv':
        df = pd.read_csv(file_name)
    elif extension=='.parquet':
        df = pq.read_table(file_name).to_pandas().reset_index()
        df.index.name = None
    
    return df
