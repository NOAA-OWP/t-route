from abc import ABC, abstractmethod
import logging
import yaml
import json
import xarray as xr
import pandas as pd
from pathlib import Path
import pyarrow.parquet as pq
from troute.nhd_network import reverse_network
from troute.nhd_network_utilities_v02 import organize_independent_networks, build_refac_connections

LOG = logging.getLogger('')

def read_diffusive_domain(domain_file):
    '''
    Read diffusive domain data from .ymal or .json file.
    
    Arguments
    ---------
    domain_file (str or pathlib.Path): Path of diffusive domain file
    
    Returns
    -------
    data (dict int: [int]): domain tailwater segments: list of segments in domain 
                            (includeing tailwater segment) 
    
    '''
    if domain_file.suffix == ".yaml":
        with open(domain_file) as domain:
            data = yaml.load(domain, Loader=yaml.SafeLoader)
    else:
        with open(domain_file) as domain:
            data = json.load(domain)
            
    return data

def read_netcdf(geo_file_path):
    '''
    Open a netcdf file with xarray and convert to dataframe
    
    Arguments
    ---------
    geo_file_path (str or pathlib.Path): netCDF filepath
    
    Returns
    -------
    ds.to_dataframe() (DataFrame): netCDF contents
    
    Notes
    -----
    - When handling large volumes of netCDF files, xarray is not the most efficient.
    
    '''
    with xr.open_dataset(geo_file_path) as ds:
        return ds.to_dataframe()

def read_parquet(file_path):
    '''
    Open a parquet file with pyarrow read_table and convert to dataframe
    
    Arguments
    ---------
    file_path (str or pathlib.Path): nparquetetCDF filepath
    
    Returns
    -------
    
    ds.to_pandas() (DataFrame): parquet contents 
    
    Notes
    -----
    - It checks if the column 'hy_id' exists in the DataFrame and if there's at least one value in that column that is not numeric (means its like wb-254530)
    

    '''
    tbl = pq.read_table(file_path)
    df = tbl.to_pandas().dropna()
    if 'hy_id' in df.columns and not df['hy_id'].str.isnumeric().all():
        df['hy_id'] = df['hy_id'].apply(lambda x: x.split('-')[-1])
        df = df[df['cs_id'] == 1]
    return df


class AbstractRouting(ABC):
    """
    
    """
    __slots__ = ["hybrid_params", "_diffusive_domain", "_coastal_boundary_depth_df",
                "_diffusive_network_data", "_topobathy_df", "_refactored_diffusive_domain",
                "_refactored_diffusive_network_data", "_refactored_reaches", 
                "_unrefactored_topobathy_df",]
    
    def __init__(self):
        """
        
        """
        self._diffusive_domain                  = None
        self._diffusive_network_data            = None
        self._topobathy_df                      = pd.DataFrame()
        self._unrefactored_topobathy_df         = pd.DataFrame() 
        self._refactored_diffusive_domain       = None
        self._refactored_diffusive_network_data = None   
        self._refactored_reaches                = {}

    @abstractmethod
    def update_routing_domain(self, dataframe, connections):
        pass
    
    @property
    @abstractmethod
    def diffusive_network_data(self):
        pass

    @property
    @abstractmethod
    def topobathy_df(self):
        pass
    
    @property
    @abstractmethod
    def refactored_diffusive_domain(self):
        pass
    
    @property
    @abstractmethod
    def refactored_reaches(self):
        pass
    
    @property
    @abstractmethod
    def unrefactored_topobathy_df(self):
        pass


class MCOnly(AbstractRouting):

    def __init__(self, _):
        self.hybrid_params = None
        
        super().__init__()

    def update_routing_domain(self, dataframe, connections):
        return dataframe, connections

    @property
    def diffusive_network_data(self):
        return self._diffusive_network_data

    @property
    def topobathy_df(self):
        return self._topobathy_df
    
    @property
    def refactored_diffusive_domain(self):
        return self._refactored_diffusive_domain
    
    @property
    def refactored_reaches(self):
        return self._refactored_reaches
    
    @property
    def unrefactored_topobathy_df(self):
        return self._unrefactored_topobathy_df


class MCwithDiffusive(AbstractRouting):

    def __init__(self, hybrid_params):
        self.hybrid_params = hybrid_params
        
        super().__init__()
    
    def update_routing_domain(self, dataframe, connections):
        #==========================================================================
        # build diffusive domain data and edit MC domain data for hybrid simulation
        '''
        # temporary code to build diffusive_domain using given IDs of head segment and tailwater segment of mainstem.
        # Lower Colorado TX mainstem
        #headlink_mainstem = 2421017
        #twlink_mainstem   = 2421105 
        # Guadalupe River TX mainstem
        #headlink_mainstem = 2405287
        #twlink_mainstem   = 2405359  
        # Brazos River TX mainstem
        #headlink_mainstem = 2408552
        #twlink_mainstem   = 2408609
        # Brays Bayou TX mainstem
        #headlink_mainstem = 2398553
        #twlink_mainstem   = 2398446
        # Goose Creek TX mainstem
        #headlink_mainstem = 2404194
        #twlink_mainstem   = 2404199 
        # Cedar Bayou, TX mainstem
        headlink_mainstem = 2404205
        twlink_mainstem   = 2404220          
        
        uslink_mainstem = headlink_mainstem
        dslink_mainstem = 1 # initial value
        mainstem_list =[headlink_mainstem]
        while dslink_mainstem != twlink_mainstem:
            dslink_mainstem = connections[uslink_mainstem][0]
            mainstem_list.append(dslink_mainstem)
            uslink_mainstem = dslink_mainstem   
        diffusive_domain={}
        diffusive_domain[twlink_mainstem] = mainstem_list
        '''
        domain_file = self.hybrid_params.get("diffusive_domain", None)
        self._diffusive_domain = read_diffusive_domain(domain_file)
        self._diffusive_network_data = {}
    
        rconn_diff0 = reverse_network(connections)
        
        for tw in self._diffusive_domain:
            mainstem_segs = self._diffusive_domain[tw]['links']
            # we want mainstem_segs start at a mainstem link right after the upstream boundary mainstem link, which is
            # in turn not under any waterbody. This boundary mainstem link should be turned into a tributary segment.
            upstream_boundary_mainstem_link = self._diffusive_domain[tw]['upstream_boundary_link_mainstem']         
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
            self._diffusive_network_data[tw]['connections'] = {k: connections[k] for k in (mainstem_segs + trib_segs)}       

            # make sure that no downstream link below tw
            self._diffusive_network_data[tw]['connections'][tw] = []
            
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
            self._diffusive_network_data[tw]['param_df'] = dataframe.filter(
                (mainstem_segs + trib_segs),
                axis = 0,
            )
            self._diffusive_network_data[tw]['upstream_boundary_link'] = upstream_boundary_mainstem_link

            # ==== remove diffusive domain segs from MC domain ====        
            # drop indices from param_df
            dataframe = dataframe.drop(mainstem_segs)

            # remove keys from connections dictionary
            for s in mainstem_segs:
                connections.pop(s)

            # update downstream connections of trib segs
            for us in trib_segs:
                connections[us] = []
        
        return dataframe, connections

    @property
    def diffusive_network_data(self):
        return self._diffusive_network_data

    @property
    def topobathy_df(self):
        return self._topobathy_df
    
    @property
    def refactored_diffusive_domain(self):
        return self._refactored_diffusive_domain
    
    @property
    def refactored_reaches(self):
        return self._refactored_reaches
    
    @property
    def unrefactored_topobathy_df(self):
        return self._unrefactored_topobathy_df


class MCwithDiffusiveNatlXSectionNonRefactored(MCwithDiffusive):

    def __init__(self, hybrid_params):

        super().__init__(hybrid_params = hybrid_params)

    @property
    def topobathy_df(self):
        if self._topobathy_df.empty:
            topobathy_file = self.hybrid_params.get("topobathy_domain", None)
            if topobathy_file.suffix == '.nc':
                self._topobathy_df = read_netcdf(topobathy_file).set_index('link')
                self._topobathy_df.index = self._topobathy_df.index.astype(int)

            elif topobathy_file.suffix == '.parquet':
                self._topobathy_df = read_parquet(topobathy_file).set_index('hy_id')
                self._topobathy_df.index = self._topobathy_df.index.astype(int)

                #If any diffusive mainstem segments doesn't have channel bathy date in topobathy_df,
                #estimate one from adjacent segments with available bathy data
                for tw in self._diffusive_domain:
                    for mainstem_segment in self._diffusive_network_data[tw]['mainstem_segs']:
                        if mainstem_segment not in self._topobathy_df.index:
                            # Temp.Solution: when topobaty is not available, use available topobathy of the closest upstream segment
                            upstream_mainstem_segment = set(self._diffusive_network_data[tw]['rconn'][mainstem_segment]).intersection(
                                                                                    self._diffusive_network_data[tw]['mainstem_segs'])
                            upstream_mainstem_segment = list(upstream_mainstem_segment)[0]
                            temp_df = self._topobathy_df[self._topobathy_df.index==upstream_mainstem_segment]
                            new_index = pd.Index([mainstem_segment]*len(temp_df))
                            temp_df.index = new_index
                            cs_id_max = temp_df['cs_id'].max()
                            fill_in_topobathy_df = temp_df[temp_df.cs_id==cs_id_max]
                            combined_df = pd.concat([self._topobathy_df, fill_in_topobathy_df])
                            self._topobathy_df = combined_df

        return self._topobathy_df


class MCwithDiffusiveNatlXSectionRefactored(MCwithDiffusive):

    def __init__(self, hybrid_params):
        
        super().__init__(hybrid_params = hybrid_params)

    @property
    def topobathy_df(self):
        if self._topobathy_df.empty:
            refactored_topobathy_file = self.hybrid_params.get("refactored_topobathy_domain", None)
            self._topobathy_df = read_netcdf(refactored_topobathy_file).set_index('link')
        return self._topobathy_df
    
    @property
    def refactored_diffusive_domain(self):
        if not self._refactored_diffusive_domain:
            refactored_domain_file = self.hybrid_params.get("refactored_domain", None)
            self._refactored_diffusive_domain = read_diffusive_domain(refactored_domain_file)
        return self._refactored_diffusive_domain
    
    @property
    def refactored_reaches(self):
        if not self._refactored_reaches:
            refactored_topobathy_file = self.hybrid_params.get("refactored_topobathy_domain", None)
            diffusive_parameters = {'geo_file_path': refactored_topobathy_file}
            refactored_connections = build_refac_connections(diffusive_parameters)

            for tw in self._diffusive_domain:

                # list of stream segments of a single refactored diffusive domain 
                refac_tw = self.refactored_diffusive_domain[tw]['refac_tw']
                rlinks_tw = self.refactored_diffusive_domain[tw]['rlinks']
                refactored_connections_tw = {}   

                # Subset a connection dictionary (upstream segment as key : downstream segments as values) from refactored_connections
                # for a single refactored diffusive domain defined by a current tw. 
                for k in rlinks_tw:
                    if k in refactored_connections.keys() and k != refac_tw:
                        refactored_connections_tw[k] = refactored_connections[k]
                
                trib_segs = self.diffusive_network_data[tw]['tributary_segments']
                refactored_diffusive_network_data = {}
                refactored_diffusive_network_data[refac_tw] = {}                
                refactored_diffusive_network_data[refac_tw]['tributary_segments'] = trib_segs
                refactored_diffusive_network_data[refac_tw]['connections'] = refactored_connections_tw                 

                for k in trib_segs:
                    refactored_diffusive_network_data[refac_tw]['connections'][k] = [self._refactored_diffusive_domain[tw]['incoming_tribs'][k]]

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
        return self._refactored_reaches
    
    @property
    def unrefactored_topobathy_df(self):
        if self._unrefactored_topobathy_df.empty:
            topobathy_file = self.hybrid_params.get("topobathy_domain",   None)
            self._unrefactored_topobathy_df = read_netcdf(topobathy_file).set_index('link')
            self._unrefactored_topobathy_df.index = self._unrefactored_topobathy_df.index.astype(int)
        return self._unrefactored_topobathy_df


