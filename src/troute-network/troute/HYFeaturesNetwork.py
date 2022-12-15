from .AbstractNetwork import AbstractNetwork
import pandas as pd
import numpy as np
import time
import troute.nhd_io as nhd_io #FIXME
import troute.hyfeature_preprocess as hyfeature_prep

__verbose__ = False
__showtiming__ = False


class HYFeaturesNetwork(AbstractNetwork):
    """
    
    """
    __slots__ = ["_flowpath_dict", 
                 "segment_index",
                 ]
    def __init__(self, 
                 supernetwork_parameters, 
                 waterbody_parameters,
                 data_assimilation_parameters,
                 restart_parameters=None, 
                 compute_parameters=None, 
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
         self._connections,
         self._waterbody_df,
         self._waterbody_types_df,
         self._terminal_codes,
        ) = hyfeature_prep.read_geo_file(
            supernetwork_parameters,
            waterbody_parameters,
        )
        
        if __verbose__:
            print("supernetwork connections set complete")
        if __showtiming__:
            print("... in %s seconds." % (time.time() - start_time))
            
        cols                         = supernetwork_parameters.get('columns',None)
        break_network_at_waterbodies = waterbody_parameters.get("break_network_at_waterbodies", False)        
        streamflow_da = data_assimilation_parameters.get('streamflow_da', False)
        break_network_at_gages       = False       
        if streamflow_da:
            break_network_at_gages   = streamflow_da.get('streamflow_nudging', False)
        break_points                 = {"break_network_at_waterbodies": break_network_at_waterbodies,
                                        "break_network_at_gages": break_network_at_gages}

        super().__init__(
            compute_parameters, 
            waterbody_parameters,
            restart_parameters, 
            cols, 
            break_points,
            verbose=__verbose__,
            showtiming=__showtiming__,
            )   
            
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

