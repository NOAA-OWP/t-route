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


class NHDNetwork(AbstractNetwork):
    """
    
    """
    __slots__ = [  
        "_link_lake_crosswalk", "_usgs_lake_gage_crosswalk", 
        "_usace_lake_gage_crosswalk", "_flowpath_dict"
        ]
    
    def __init__(
                self, 
                supernetwork_parameters, 
                waterbody_parameters=None, 
                restart_parameters=None, 
                forcing_parameters=None, 
                compute_parameters=None, 
                data_assimilation_parameters=None, 
                preprocessing_parameters=None, 
                verbose=False, 
                showtiming=False,
                ):
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
        # Load Geo Data
        #------------------------------------------------

        (
            self._dataframe,
            self._connections,
            self._terminal_codes,
            self._waterbody_df, 
            self._waterbody_types_df,
            self._waterbody_type_specified,
            self._waterbody_connections, 
            self._link_lake_crosswalk,
            self._gages,
            self._usgs_lake_gage_crosswalk, 
            self._usace_lake_gage_crosswalk,
        ) = nhd_prep.read_geo_file(
            supernetwork_parameters,
            waterbody_parameters,
            data_assimilation_parameters,
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
        
        self._flowpath_dict = {}

        super().__init__(
            compute_parameters, 
            waterbody_parameters,
            restart_parameters,
            break_points,
            verbose=__verbose__,
            showtiming=__showtiming__,
            )
        
        # list of all segments in the domain (MC + diffusive)
        self.segment_index = self._dataframe.index
        if self.diffusive_network_data:
            for tw in self.diffusive_network_data:
                self.segment_index = self.segment_index.append(
                    pd.Index(self.diffusive_network_data[tw]['mainstem_segs'])
                ) 

        # Create empty dataframe for coastal_boundary_depth_df. This way we can check if
        # it exists, and only read in SCHISM data during 'assemble_forcings' if it doesn't
        self._coastal_boundary_depth_df = pd.DataFrame()

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

    #@property
    #def wbody_conn(self):
    #    return self._waterbody_connections    
