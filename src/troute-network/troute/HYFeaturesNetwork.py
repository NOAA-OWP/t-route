from .AbstractNetwork import AbstractNetwork
import pathlib
import json
import pandas as pd
import numpy as np
import time
import re
import troute.nhd_io as nhd_io #FIXME
from itertools import chain

__verbose__ = False
__showtiming__ = False

def node_key_func_nexus(x):
    return int(x[4:])

def node_key_func_wb(x):
    return int(x[3:])

def read_ngen_waterbody_df(parm_file, lake_index_field="wb-id", lake_id_mask=None):
    """
    Reads lake.json file and prepares a dataframe, filtered
    to the relevant reservoirs, to provide the parameters
    for level-pool reservoir computation.
    """
    def node_key_func(x):
        return int(x[3:])
    df = pd.read_json(parm_file, orient="index")

    df.index = df.index.map(node_key_func)
    df.index.name = lake_index_field
    #df = df.set_index(lake_index_field, append=True).reset_index(level=0)
    #df.rename(columns={'level_0':'wb-id'}, inplace=True)
    if lake_id_mask:
        df = df.loc[lake_id_mask]
    return df

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
    #without consuming the generator...Otherwise, if nuxus_files
    #is empty, the following concat raises a ValueError which
    #It may be sufficient to catch this exception and warn that
    #there may not be any pattern matching files in the dir

    #pd.concat((pd.read_csv(f, index_col=0, usecols=[1,2], header=None, engine='c').rename(columns={2:id_regex.match(f).group(1)}) for f in all_files[0:2]), axis=1).T
    id_regex = re.compile(r".*nex-(\d+)_.*.csv")
    nexuses_flows_df = pd.concat(
            #Read the nexus csv file
            (pd.read_csv(f, index_col=0, usecols=[1,2], header=None, engine='c').rename(
                #Rename the flow column to the id of the nexus
                columns={2:int(id_regex.match(f.name).group(1))})
            for f in nexus_files #Build the generator for each required file
            ),  axis=1).T #Have now concatenated a single df (along axis 1).  Transpose it.
    missing = nexuses_flows_df[ nexuses_flows_df.isna().any(axis=1) ]
    if  not missing.empty:
        raise ValueError("The following nexus inputs are incomplete: "+str(missing.index))
    rt1 = time.time()
    #print("Time to build nexus_flows_df: {} seconds".format(rt1-rt0))

    qlat_df = pd.concat( (nexuses_flows_df.loc[int(k)].rename(index={int(k):v})
        for k,v in nexus_to_downstream_flowpath_dict.items() ), axis=1
        ).T
    
    # The segment_index has the full network set of segments/flowpaths. 
    # Whereas the set of flowpaths that are downstream of nexuses is a 
    # subset of the segment_index. Therefore, all of the segments/flowpaths
    # that are not accounted for in the set of flowpaths downstream of
    # nexuses need to be added to the qlateral dataframe and padded with
    # zeros.

    #tq0 = time.time()
    all_df = pd.DataFrame( np.zeros( (len(segment_index), len(qlat_df.columns)) ), index=segment_index,
            columns=qlat_df.columns )

    all_df.loc[ qlat_df.index ] = qlat_df

    #tq1 =  time.time()
    #print("Time to fill null qlats: {}".format(tq1-tq0))
    # Sort qlats
    qlat_df = all_df.sort_index()

    # Set new nts based upon total nexus inputs
    nts = (qlat_df.shape[1]) * qts_subdivisions

    max_col = 1 + nts // qts_subdivisions

    if len(qlat_df.columns) > max_col:
        qlat_df.drop(qlat_df.columns[max_col:], axis=1, inplace=True)

    if __verbose__:
        print("qlateral array complete")
    if __showtiming__:
        print("... in %s seconds." % (time.time() - start_time))

    return qlat_df

def read_nexus_file(nexus_file_path):

    #Currently reading data in format:
    #[
       #{
         #"ID": "wb-44",
         #"toID": "nex-45"
         #},
    with open(nexus_file_path) as data_file:
        json_data_list = json.load(data_file)

    nexus_to_downstream_flowpath_dict_str = {}

    for id_dict in json_data_list:
        if "nex" in id_dict['ID']:
            nexus_to_downstream_flowpath_dict_str[id_dict['ID']] = id_dict['toID']

    # Extract the ID integer values
    nexus_to_downstream_flowpath_dict = {node_key_func_nexus(k): node_key_func_wb(v) for k, v in nexus_to_downstream_flowpath_dict_str.items()}

    return nexus_to_downstream_flowpath_dict

def read_json(file_path):
    dfs = []
    with open(file_path) as data_file:
        json_data = json.load(data_file)

        def node_key_func(x):
            return int(x[3:])
        
        for key_wb, value_params in json_data.items():
            df = pd.json_normalize(value_params)
            df['ID'] = node_key_func(key_wb)
            dfs.append(df)
    df_main = pd.concat(dfs, ignore_index=True)

    return df_main

class HYFeaturesNetwork(AbstractNetwork):
    """
    
    """
    __slots__ = ["_flowpath_dict"]
    def __init__(self, supernetwork_parameters, waterbody_parameters=None, restart_parameters=None, forcing_parameters=None, verbose=False, showtiming=False):
        """
        
        """
        global __verbose__, __showtiming__
        __verbose__ = verbose
        __showtiming__ = showtiming
        if __verbose__:
            print("creating supernetwork connections set")
        if __showtiming__:
            start_time = time.time()
        geo_file_path = supernetwork_parameters["geo_file_path"]
        nexus_file_path = supernetwork_parameters["ngen_nexus_file"]
        cols = supernetwork_parameters["columns"]
        terminal_code = supernetwork_parameters.get("terminal_code", 0)
        break_network_at_waterbodies = supernetwork_parameters.get(
            "break_network_at_waterbodies", False
        )
        break_network_at_gages = supernetwork_parameters.get(
            "break_network_at_gages", False
        )
        break_points = {"break_network_at_waterbodies": break_network_at_waterbodies,
                        "break_network_at_gages": break_network_at_gages}
        #df = pd.read_json(geo_file_path, orient="index")
        #print(df)
        #raise( "BREAK" )
        self._dataframe = read_json(geo_file_path)
        #df["NHDWaterbodyComID"].fillna(-9999, inplace=True)
        #df["NHDWaterbodyComID"] = df["NHDWaterbodyComID"].astype("int64")
        self._flowpath_dict = read_nexus_file(
            pathlib.Path(nexus_file_path)
        )
        self._waterbody_types_df = pd.DataFrame()
        self._waterbody_df = pd.DataFrame()
        #FIXME the base class constructor is finiky
        #as it requires the _dataframe, then sets some 
        #initial default properties...which, at the moment
        #are used by the subclass constructor.
        #So it needs to be called at just the right spot...
        super().__init__(cols, terminal_code, break_points)
        #FIXME once again, order here can hurt....to hack `alt` in, either need to
        #put it as a column in the config, or do this AFTER the super constructor
        #otherwise the alt column gets sliced out...
        self._dataframe['alt'] = 1.0 #FIXME get the right value for this...
        #Load waterbody/reservoir info
        #For ngen HYFeatures, the reservoirs to be simulated
        #are determined by the lake.json file
        #we limit waterbody_connections to only the flowpaths
        #that coincide with a lake listed in this file
        #see `waterbody_connections`

        if waterbody_parameters:
            #FIXME later, DO ALL LAKE PARAMS BETTER
            levelpool_params = waterbody_parameters.get('level_pool', None)
            if not levelpool_params:
                #FIXME should not be a hard requirement
                raise(RuntimeError("No supplied levelpool parameters in routing config"))
            
            lake_id = levelpool_params.get("level_pool_waterbody_id", "wb-id")
            self._waterbody_df = read_ngen_waterbody_df(
                levelpool_params["level_pool_waterbody_parameter_file_path"],
                lake_id,
                #self.waterbody_connections.values()
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
            try: #FIXME for HYFeatures/ngen this will likely be a lot different...
                self._waterbody_types_df = nhd_io.read_reservoir_parameter_file(
                    hybrid_params["reservoir_parameter_file"],
                    lake_id,
                    #self.waterbody_connections.values(),
                )
                # Remove duplicate lake_ids and rows
                self._waterbody_types_df = (
                self._waterbody_types_df.reset_index()
                .drop_duplicates(subset=lake_id)
                .set_index(lake_id)
                )
            except:
                self._waterbody_types_df = pd.DataFrame(index=self._waterbody_df.index)
                self._waterbody_types_df['reservoir_type'] = 1
                #FIXME any reservoir operations requires some type
                #So make this default to 1 (levelpool)
        #At this point, need to adjust some waterbody/channel parameters based on lakes/reservoirs
        adjust = [ zip(x, y) 
                   for x, y in 
                   zip(self._waterbody_df['member_wbs'], self._waterbody_df['partial_length_percent'])
                ]
        #adjust is a generator of a list of list of tuples...use chain to flatten
        for wb, percent in chain.from_iterable(adjust):
            #print(wb, percent)
            #FIXME not sure why some of these are 100%, if that  is the case
            #shouldn't they just not be in the topology???
            wb = node_key_func_wb(wb)
            #Need to adjust waterbodys/channels that  interact with this waterbody
            #print(self._dataframe.loc[wb, 'dx'])
            self._dataframe.loc[wb, 'dx'] = self._dataframe.loc[wb, 'dx'] - self._dataframe.loc[wb, 'dx']*percent
            #print(self._dataframe.loc[wb, 'dx'])

        if __verbose__:
            print("supernetwork connections set complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))   
        #TODO/FIXME reservoir restart load
        #TODO/FIXME channel restart load
        # STEP 4: Handle Channel Initial States
        if __showtiming__:
            start_time = time.time()
        if __verbose__:
            print("setting channel initial states ...")
        #Get channel restarts
        channel_restart_file = restart_parameters.get("channel_restart_file", None)

        if channel_restart_file:
            self._q0 = nhd_io.get_channel_restart_from_csv(channel_restart_file)
            self._q0 = self._q0[self._q0.index.isin(self._dataframe.index)]
            # TODO is this the same???
            #self._q0 = self._q0.loc[self._dataframe.index]
        #TODO/FIXME channel restart t0? self._t0 = ???
        if __verbose__:
          print("channel initial states complete")
        if __showtiming__:
          print("... in %s seconds." % (time.time() - start_time))
        
        self._qlateral = read_qlats(forcing_parameters, self._dataframe.index, self.downstream_flowpath_dict)

        #FIXME should waterbody_df and param_df overlap IDS?  Doesn't seem like it should...
        #self._dataframe.drop(self._waterbody_df.index, axis=0, inplace=True)
        #For now, doing it in property waterbody_connections...
        # Remove duplicate rows...but its not really duplicate...
        #self._dataframe = (
        #    self._dataframe.reset_index()
        #    .drop_duplicates(subset="key")
        #    .set_index("key")
        #)

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
            waterbody_segments = self._dataframe.dropna().loc[
                self._dataframe["waterbody"] != self.waterbody_null, "waterbody"
            ]
            waterbody_segments = waterbody_segments.loc[self.waterbody_dataframe.index]
            self._waterbody_connections = waterbody_segments.index\
                .to_series(name = waterbody_segments.name)\
                .astype("int")\
                .to_dict()
            #If we identify as a waterbody, drop from the main dataframe
            #Ugh, but this drops everything that that MIGHT be a "lake"
            #without knowing if it was defined as a lake in the lake params
            #so this should just drop the waterbody_df index, not these segments...
            #In fact, these waterbody connections should probably be entirely reworked
            #with that in mind...
            self._waterbody_connections = self._waterbody_df.index.to_series(name = self._waterbody_df.index.name).astype("int").to_dict()
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
        return pd.NA
