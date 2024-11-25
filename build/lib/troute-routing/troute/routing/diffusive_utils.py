import numpy as np
from functools import partial, reduce
import troute.nhd_network as nhd_network
from datetime import datetime, timedelta
import pandas as pd
from collections import Counter
import math


def adj_alt1(
    mx_jorder, ordered_reaches, param_df, dbfksegID, z_all
):
    """
    Adjust reach altitude data so that altitude of last node in reach is equal to that of the head segment of
    the neighboring downstream reach. 
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    param_df --(DataFrame) geomorphic parameters
    dbfksegID -- (int) segment ID of fake (ghost) node at network downstream boundary 
    z_all -- (dict) adjusted altitude dictionary with placeholder values to be replaced
    
    Returns
    ----------
    z_all -- (dict) adjusted altitude dictionary
    """

    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            for seg in range(0, ncomp):
                segID = seg_list[seg]
                if seg == ncomp - 1 and seg_list.count(dbfksegID) == 0:
                    # At junction, the altitude of bottom node of an upstream reach
                    # is equal to that of the first segment of the downstream reach

                    # head segment id of downstream reach after a junction
                    dsrchID = reach["downstream_head_segment"]                    
                    z_all[segID]["adj.alt"][0] = param_df.loc[dsrchID, 'alt']

                elif seg == ncomp - 1 and seg_list.count(dbfksegID) > 0:
                    # channel slope-adjusted bottom elevation at the bottom node of TW reach
                    ## AD HOC: need to be corrected later
                    segID2 = seg_list[seg - 1]                    
                    z_all[segID]["adj.alt"][0] = z_all[segID2]["adj.alt"][0] - param_df.loc[segID2, 's0'] * param_df.loc[segID2, 'dx']
                    
                else:                    
                    z_all[segID]["adj.alt"][0] = param_df.loc[segID, 'alt']

    return z_all


def fp_network_map(
                    mainstem_seg_list, 
                    mx_jorder, 
                    ordered_reaches, 
                    rchbottom_reaches, 
                    nrch_g, 
                    frnw_col, 
                    dbfksegID, 
                    pynw,
                    #upstream_boundary_link,
):
    """
    Channel network mapping between Python and Fortran
    
    Parameters
    ----------
    mainstem_seg_list - (int) a list of link IDs of segments of  mainstem reaches solely for diffusive routing 
    mx_jorder -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    rchbottom_reaches -- (dict) reaches and reach metadata keyed by reach tail segment
    nrch_g -- (int) number of reaches in the network
    frnw_col -- (int) number of columns in the fortran network map
    dbfskegID -- (str) segment ID of fake (ghost) node at network downstream boundary 
    pynw -- (dict) ordered reach head segments
    
    Returns
    -------
    frnw_g -- (nparray of int) Fortran-Python network mapping array
    ub_idx_fr -- (int) Fortran reach index marking a link of mainstem at the upstream boundary
    """

    #  Store headwater reach and upstream reaches above a junction
    #  as well as downstream reach after a junction
    #  into python-extension-fortran variables.
    frnw_g = np.zeros((nrch_g, frnw_col), dtype='int32')
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            frnw_g[frj, 0] = ncomp
         
            if not reach["upstream_bottom_segments"]:
                # for diffusive mainstem, this case can be tributaries at the conflunces with the mainstem or mainstem link at its upstream boundary
                frnw_g[frj, 2] = 0  # to mark these tributaries
            else:
                # reaches before a junction
                nusrch = len(reach["upstream_bottom_segments"])
                frnw_g[frj, 2] = nusrch  # the number of upstream reaches
                usrch_bseg_list = list(reach["upstream_bottom_segments"])
                #usrch_hseg_mainstem_j=-100 # ini.value of frj* of reach on mainstem that is just upstream of the current reach of frj
                #frnw_g[frj, 2 + nusrch + 1] = usrch_hseg_mainstem_j + 1
                i = 0
                r = 0
                for usrch in range(0, nusrch):
                    usrch_bseg_id = usrch_bseg_list[
                        usrch
                    ]  # upstream reach's bottom segment
                    usrch_hseg_id = rchbottom_reaches[usrch_bseg_id]["segments_list"][0]
                    # find Fortran js corresponding to individual usrchid
                    for j, sid in pynw.items():
                        if sid == usrch_hseg_id:
                            i = i + 1
                            frnw_g[frj, 2 + i] = j
                            # find if the selected headseg ID of an upstream reach belong to mainstem segment
                            #if mainstem_seg_list.count(usrch_hseg_id) > 0:
                            #    r = r+1
                            #    usrch_hseg_mainstem_j=j
                                # store frj of mainstem headseg's reach in the just upstream of the current reach of frj
                            #    frnw_g[frj, 2 + nusrch + r] = usrch_hseg_mainstem_j + 1
              
            if seg_list.count(dbfksegID) > 0:
                # a reach where downstream boundary condition is set.
                frnw_g[
                    frj, 1
                ] = -100  # head_segment ID that is in terminal downstream reach.
                # That is, -100 indicates the reach of the head segment is
                # terminal downstream reach where ds.bcond. happens.
            else:
                # reach after a junction
                dsrch_hseg_id = reach["downstream_head_segment"]
                # fortran j index equivalent to dsrchID.
                frnw_g[frj, 1] = [
                    j for j, sid in pynw.items() if sid == dsrch_hseg_id[0]
                ][0]

    # Adust frnw_g element values according to Fortran-Python index relationship, that is Python i = Fortran i+1
    for frj in range(0, nrch_g):
        frnw_g[frj, 1] = frnw_g[frj, 1] + 1  # downstream reach index for frj reach
        if frnw_g[frj, 2] > 0:
            nusrch = frnw_g[frj, 2]
            for i in range(0, nusrch):
                frnw_g[frj, 3 + i] = (
                    frnw_g[frj, 3 + i] + 1
                )  # upstream reach indicds for frj reach

    return frnw_g 


def fp_chgeo_map(
    mx_jorder, ordered_reaches, param_df, z_all, mxncomp_g, nrch_g
):
    """
    Channel geometry data mapping between Python and Fortran
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    param_df -- (DataFrame) geomorphic parameters
    z_all -- (dict) adjusted altitude dictionary
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    
    Returns
    -------
    z_ar_g -- (numpy of float64s) altitude (meters)
    bo_ar_g -- (numpy of float64s) bottom width (meters)
    traps_ar_g -- (numpy of float64s) sideslope (m/m)
    tw_ar_g -- (numpy of float64s) top width (meters)
    twcc_ar_g -- (numpy of float64s) top width of compound channel (meters)
    mann_ar_g -- (numpy of float64s) manning's roughness (seconds/meters^1/3)
    mancc_ar_g -- (numpy of float64s) manning's roughness compound channel (seconds/meters^(1/3))
    so_ar_g -- (numpy of float64s) bottom slope (m/m)
    dx_ar_g -- (numpy of float64s) segment length (meters)
    """

    z_ar_g = np.zeros((mxncomp_g, nrch_g))
    bo_ar_g = np.zeros((mxncomp_g, nrch_g))
    traps_ar_g = np.zeros((mxncomp_g, nrch_g))
    tw_ar_g = np.zeros((mxncomp_g, nrch_g))
    twcc_ar_g = np.zeros((mxncomp_g, nrch_g))
    mann_ar_g = np.zeros((mxncomp_g, nrch_g))
    manncc_ar_g = np.zeros((mxncomp_g, nrch_g))
    so_ar_g = np.zeros((mxncomp_g, nrch_g))
    dx_ar_g = np.zeros((mxncomp_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                if seg == ncomp - 1:
                    segID = seg_list[seg - 1]
                else:
                    segID = seg_list[seg]
                
                bo_ar_g[seg, frj] = param_df.loc[segID, 'bw']
                traps_ar_g[seg, frj] = 1/param_df.loc[segID, 'cs']
                tw_ar_g[seg, frj] = param_df.loc[segID, 'tw']
                twcc_ar_g[seg, frj] = param_df.loc[segID, 'twcc']
                mann_ar_g[seg, frj] = param_df.loc[segID, 'n']
                manncc_ar_g[seg, frj] = param_df.loc[segID, 'ncc']
                so_ar_g[seg, frj] = param_df.loc[segID, 's0']
                dx_ar_g[seg, frj] = param_df.loc[segID, 'dx']

                segID1 = seg_list[seg]
                z_ar_g[seg, frj] = z_all[segID1]["adj.alt"][0]

    return (
        z_ar_g,
        bo_ar_g,
        traps_ar_g,
        tw_ar_g,
        twcc_ar_g,
        mann_ar_g,
        manncc_ar_g,
        so_ar_g,
        dx_ar_g,
    )


def fp_qlat_map(
    mx_jorder,
    ordered_reaches,
    nts_ql_g,
    param_df,
    qlat,
    qlat_g,
):
    """
    lateral inflow mapping between Python and Fortran
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    nts_ql_g -- (int) numer of qlateral timesteps
    param_df --(DataFrame) geomorphic parameters
    qlat -- (DataFrame) qlateral data (m3/sec)
    qlat_g -- (ndarray of float32) empty qlateral array to be filled
    
    Returns
    -------
    qlat_g -- (ndarray of float32) qlateral array (m3/sec/m)
    
    Notes
    -----
    data in qlat_g are normalized by segment length with units of m2/sec = m3/sec/m
    """

    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                segID = seg_list[seg]
                for tsi in range(0, nts_ql_g):
                    if seg < ncomp - 1:
                        
                        tlf = qlat.loc[segID, tsi]
                        dx = param_df.loc[segID, 'dx']
                        qlat_g[tsi, seg, frj] = tlf / dx  # [m^2/sec]

                    else:
                        qlat_g[
                            tsi, seg, frj
                        ] = 0.0  # seg=ncomp is actually for bottom node in Fotran code.
                        # And, lateral flow enters between adjacent nodes.

    return qlat_g

def fp_ubcd_map(frnw_g, pynw, nts_ub_g, nrch_g, ds_seg, upstream_inflows):
    """
    Upstream boundary condition mapping between Python and Fortran
    
    Parameters
    ----------
    frnw_g -- (nparray of int) Fortran-Python network mapping array
    pynw -- (dict) ordered reach head segments
    nrch_g -- (int) number of reaches in the network
    ds_seg -- (ndarray of int64) row indices for downstream segments recieving flows in upstream_flows
    upstream_inflows (ndarray of float32) upstream_inflows (m3/sec) to segments in ds_seg
    
    Returns
    -------
    ubcd_g -- (ndarray of float32) upstream boundary data (m3/sec)
    """
    ubcd_g = np.zeros((nts_ub_g, nrch_g))        
    
    return ubcd_g


def fp_dbcd_map(usgsID2tw=None, usgssDT=None, usgseDT=None, usgspCd=None):
    """
    Downststream boundary condition mapping between Python and Fortran using USGS stage observations
    
    Parameters
    ----------
    usgsID2tw -- (str) usgs site ID
    usgssDT -- (str) start data request date (yyyy-mm-dd) 
    usgseDT -- (str) end data request date (yyyy-mm-dd) 
    usgspCd --  (list of str) usgs parameter code 
    
    Returns
    -------
    nts_db_g -- (int) number of timesteps in downstream boundary data
    dbcd_g -- (ndarray of float64) downstream boundary data (meters)
    """

    # ** 1) downstream stage (here, lake elevation) boundary condition
    # from nwis_client.iv import IVDataService
    # install via: pip install hydrotools.nwis_client
    if usgsID2tw:
        try:
            from hydrotools.nwis_client.iv import IVDataService
        except ImportError as err:
            print(err, end="... ")
            print(
                "Please install hydrotools.nwis_client via: `pip install hydrotools.nwis_client`"
            )
            raise  # ensures program exit by re-raising the error.

        # from evaluation_tools.nwis_client.iv import IVDataService
        # Retrieve streamflow and stage data from two sites
        # Note: 1. Retrieved data all are based on UTC time zone (UTC is 4 hours ahead of Eastern Time during
        #          daylight saving time and 5 hours ahead during standard time)
        #       2. Retrieved data are always 1 hour ahead of stated startDT and 1 hour ahead of endDT,
        #          where starDT or endDT equal to yyyy-mm-dd 00:00.
        #       3. Also, retrieved data in 15 min so there are always four more data before startDT
        #          and four less data before endDT.
        #          For example, startDT='2018-08-01' and endDT='2018-09-01', then the retrieved data starts by
        #          2018-07-31-23:00:00
        #          2018-07-31-23:15:00
        #          2018-07-31-23:30:00
        #          2018-07-31-23:45:00
        #          2018-08-01-00:00:00
        #               .......
        #          2018-08-31-22:00:00
        #          2018-08-31-22:15:00
        #          2018-08-31-22:30:00
        #          2018-08-31-22:45:00
        #          2018-08-31-23:00:00
        #       4. '00060' for discharge [ft^3/s]
        #          '00065' for stage [ft]
        #          '62614' for Elevation, lake/res,NGVD29 [ft]
        ivds = IVDataService()
        observations_data = ivds.get(
            sites=usgsID2tw,  # sites='01646500,0208758850',
            startDT=usgssDT,  #'2018-08-01',
            endDT=usgseDT,  #'2020-09-01',
            # parameterCd='62614'
            parameterCd=usgspCd,
        )
        nts_db_g = len(observations_data) - 4
        # ** 4 is added to make data used here has its date time as from startDT 00:00 to (endDT-1day) 23:00, UTC
        # ** usgs data at this site uses NGVD1929 feet datum while 'alt' of RouteLink uses NAD88 meter datum.
        #  -> Has to convert accordingly !!!!
        #
        # source: https://pubs.usgs.gov/sir/2010/5040/section.html
        # Over most USGS study area it is used that NGVD = NAVD88 - 3.6 feet
        dbcd_g = np.zeros(nts_db_g)
        for tsi in range(0, nts_db_g):
            i = tsi + 4
            dmy = observations_data.iloc[i, 4]
            dbcd_g[tsi] = 0.3048 * (
                dmy + 3.6
            )  # accuracy with +-0.5feet for 95 percent of USGS study area.
            # 0.3048 to covert ft to meter. [meter]
    else:
        nts_db_g = 1
        dbcd_g = -np.ones(nts_db_g)
        
    return nts_db_g, dbcd_g

def fp_naturalxsec_map(
                ordered_reaches,
                mainstem_seg_list, 
                topobathy_bytw,
                param_df, 
                mx_jorder, 
                mxncomp_g, 
                nrch_g,
                dbfksegID):
    """
    natural cross section mapping between Python and Fortran using eHydro_ned_cross_sections data
    
    Parameters
    ----------
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    mainstem_seg_list -- (int) a list of link IDs of segs of related mainstem reaches 
    topobathy_bytw -- (DataFrame) natural cross section's x and z values with manning's N   
    param_df --(DataFrame) geomorphic parameters
    mx_jorder -- (int) max junction order
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    dbfksegID -- (int) segment ID of fake node (=bottom node) of TW reach that hosts downstream boundary 
                        condition for a network that is being routed. 
    
    Returns
    -------
    x_bathy_g -- (numpy of float64s) lateral distance of bathy data points
    z_bathy_g -- (numpy of float64s) elevation of bathy data points
    size_bathy_g -- (integer) the nubmer of bathy data points of each cross section
    mxnbathy_g -- (integer) maximum size of bathy data points
    
    Notes
    -----
    - In node-configuration, bottom node takes bathy of the first segment of the downtream reach
      except TW reach where the bottom node bathy is interpolated by bathy of the last segment 
      with so*0.5*dx 
    """  
    if not topobathy_bytw.empty:
        
        # maximum number of stations along a single cross section
        mxnbathy_g = topobathy_bytw.index.value_counts().max()

        # initialize arrays to store cross section data
        x_bathy_g    = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
        z_bathy_g    = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
        mann_bathy_g = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
        size_bathy_g = np.zeros((mxncomp_g, nrch_g), dtype='i4')    

        # loop over reach orders.
        frj = -1
        for x in range(mx_jorder, -1, -1):

            # loop through all reaches of order x
            for head_segment, reach in ordered_reaches[x]:
                frj = frj + 1

                # list of segments in this reach
                seg_list = reach["segments_list"]

                # number of segments in this reach
                ncomp = reach["number_segments"]

                # determine if this reach is part of the mainstem diffusive domain
                if head_segment in mainstem_seg_list: 

                    # loop through segments in mainstem reach
                    for seg, segID in enumerate(seg_list):

                        # identify the index in topobathy dataframe that contains
                        # the data we want for this node.
                        if seg == ncomp-1 and x > 0: 
                            
                            # if last node of a reach, but not the last node in the network
                            # use cross section of downstream neighbor
                            seg_idx = reach["downstream_head_segment"][0]

                        elif segID == dbfksegID:
                            
                            # if last node of reach AND last node in the network,
                            # use cross section of upstream neighbor
                            seg_idx = seg_list[seg-1]
                            
                        else:
                            seg_idx = segID
                            
                        # how many stations are in the node cross section?
                        nstations = len(topobathy_bytw.loc[seg_idx])
                        
                        # populate cross section size (# of stations) array
                        size_bathy_g[seg, frj] = nstations
                        
                        # populate cross section x, z and mannings n arrays
                        x_bathy_g[0:nstations, seg, frj]    = topobathy_bytw.loc[seg_idx].xid_d
                        z_bathy_g[0:nstations, seg, frj]    = topobathy_bytw.loc[seg_idx].z
                        mann_bathy_g[0:nstations, seg, frj] = topobathy_bytw.loc[seg_idx].n
                        
                        # if terminal node of the network, then adjust the cross section z data using
                        # channel slope and length data
                        if segID == dbfksegID:
                            So = param_df.loc[seg_idx].s0
                            dx = param_df.loc[seg_idx].dx
                            z_bathy_g[0:nstations, seg, frj] = z_bathy_g[0:nstations, seg, frj] - So * dx   
   
    else: 
        #if the bathy dataframe is empty, then pass out empty arrays
        x_bathy_g    = np.array([]).reshape(0,0,0)
        z_bathy_g    = np.array([]).reshape(0,0,0)
        mann_bathy_g = np.array([]).reshape(0,0,0)
        size_bathy_g = np.array([], dtype = 'i4').reshape(0,0)
        mxnbathy_g   = int(0)
    
    return x_bathy_g, z_bathy_g, mann_bathy_g, size_bathy_g, mxnbathy_g

def fp_da_map(
    mx_jorder,
    ordered_reaches,
    usgs_df,
    nrch_g,
    t0,
    nsteps,
    dt_da_g,
    t0_g,
    tfin_g):
    """
    Data assimilatoin data mapping between Python and Fortran

    Parameters
    ----------
    mx_jorder       -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    usgs_df         -- (DataFrame) usgs streamflow data interpolated at dt time steps
    t0              -- (datetime) initial date
    nsteps          -- (int) number of simulation time steps
    dt_da_g         -- (int) numer of Data Assimilation timesteps
    t0_g            -- (float) diffusive model's initial simulation time [hr] (by default, zero)
    tfin_g          -- (float) diffusive model's final simulation time [hr] 

    Returns
    -------
    nts_da_g        -- (int) number of DA timesteps
    usgs_da_g       -- (float) usgs oberved streamflow data [cms]
    usgs_da_reach_g -- (int) indices of stream reaches where DA is applied    
    """ 
    nts_da_g    = int((tfin_g - t0_g) * 3600.0 / dt_da_g) + 1  # include initial time 0 to the final time    
    usgs_da_g   = -4444.0*np.ones((nts_da_g, nrch_g))
    usgs_da_reach_g = np.zeros(nrch_g, dtype='i4')
    
    if not usgs_df.empty:    
        dt_timeslice = timedelta(minutes=dt_da_g/60.0)
        tfin         = t0 + dt_timeslice*nsteps
        timestamps   = pd.date_range(t0, tfin, freq=dt_timeslice)
    
        usgs_df_complete = usgs_df.replace(np.nan, -4444.0)
    
        for i in range(len(timestamps)):
            if timestamps[i] not in usgs_df.columns:
                usgs_df_complete.insert(i, timestamps[i], -4444.0*np.ones(len(usgs_df)), allow_duplicates=False)
  
        frj = -1
        for x in range(mx_jorder, -1, -1):
            for head_segment, reach in ordered_reaches[x]:
                seg_list = reach["segments_list"]
                ncomp = reach["number_segments"]
                frj = frj + 1
                for seg in range(0, ncomp):
                    segID = seg_list[seg]
                    if segID in usgs_df_complete.index:
                        usgs_da_g[:,frj] = usgs_df_complete.loc[segID].values[0:nts_da_g]
                        usgs_da_reach_g[frj] = frj + 1  # Fortran-Python index relationship, that is Python i = Fortran i+1
                       
    return nts_da_g, usgs_da_g, usgs_da_reach_g

def fp_refactored_network_map(
    nrch_g,
    frnw_g, 
    pynw, 
    mainstem_seg_list,
    trib_seg_list,
    refactored_diffusive_domain,
    refactored_reaches,
    #upstream_boundary_link
):
    """
    Data assimilatoin data mapping between Python and Fortran
    
    Parameters
    ----------
    nrch_g -- (int) number of reaches in the network
    frnw_g -- (nparray of int) Fortran-Python network mapping array
    pynw -- (dict) ordered reach head segments
    refactored_diffusive_domain -- (dict) geometric relationship information between original and refactored hydrofabrics
    refactored_reaches -- (list of lists) lists of stream segment IDs of diffusive mainstems on refactored hydrofabrics including 
                                          tributaries of original hydrofabric

    Returns
    -------
    rfrnw_g -- (nparray of int) Fortran-Python network mapping array on refactored hydrofabric
    rpynw -- (dict) ordered reach head segments on refactored hydrofabric
    rordered_reaches -- (dict) reaches and reach metadata by junction order on refactored hydrofabric
    """ 
    rfrnw_g = frnw_g.copy()
    rpynw = {}

    for frj in range(0, nrch_g):
        if frnw_g[frj, 2] == 0: # number of upstream reaches = 0 means it's a tributary entering mainstem
            rpynw[frj] = pynw[frj]
        elif frnw_g[frj, 2] > 0: # number of upstream reaches > 0 means it's a mainstem
            rpynw[frj] = refactored_diffusive_domain['incoming_tribs'][pynw[frj-1]]

    head_rlink_frequency = Counter(refactored_diffusive_domain['incoming_tribs'].values())
    head_rlink_multi_trib_links = {key: value for key,value in head_rlink_frequency.items() if value>1}
    
    # make sure multiple tributaries into a single rlink are seperated out by a mainstem link inbetween. Refer to edge cases,p52,RM6
    reversed_pynw = {pynw[k]:k for k in pynw}
    head_rlink_trial = list(head_rlink_multi_trib_links.keys()) 
    head_rlink_2b_removed = []
    
    for i in range(len(head_rlink_trial)):
        rlnk = head_rlink_trial[i]
        trib_link = [k for k, v in refactored_diffusive_domain['incoming_tribs'].items() if v == rlnk]
        # place tributary links in junction order
        trib_link_ordered = []
        for lnk in trib_seg_list:
            if lnk in trib_link:
                trib_link_ordered.append(lnk) 
        
        n_mseg = 0
        for i2 in range(len(trib_link_ordered)-1):
            junct_order1 = reversed_pynw[trib_link_ordered[i2]]
            junct_order2 = reversed_pynw[trib_link_ordered[i2+1]]
            if abs(junct_order1 - junct_order2) == 1:
            # when no mainstem link is between these two tribs. That is, the two tribs both enter a mainstem link at the upper end. 
                head_rlink_2b_removed.append(rlnk)
                #if rlnk in head_rlink_trial:
                #    head_rlink_multi_trib_links.pop(rlnk)
                #    head_rlink_trial.remove(rlnk)
            else:
                # make sure a mainstem link exist inbetween related tribs
                 for odr in range(junct_order1+1, junct_order2):
                    lnk = pynw[odr]
                    if lnk in mainstem_seg_list:
                        n_mseg = n_mseg + 1      
        
        if n_mseg == 0:
            # when no mainstem link between multiple trib links
            head_rlink_2b_removed.append(rlnk)
            #if rlnk in head_rlink_trial:
            #    head_rlink_multi_trib_links.pop(rlnk)            
        
    unique_rlinks = list(set(head_rlink_2b_removed))
    for rlnk in unique_rlinks:
        head_rlink_multi_trib_links.pop(rlnk)

    # make a dictionary of  rlink head (as a key): [rlinks] in the reach (as a value)
    rordered_reaches = {}
    for frj in range(0, nrch_g):
        rreach = [x for x in refactored_reaches if rpynw[frj] in x][0]
        rordered_reaches[rpynw[frj]] = rreach

    # find the number of nodes for refactored hydrofabric
    for frj in range(0, nrch_g):
        rfrnw_g[frj, 0] = len(rordered_reaches[rpynw[frj]]) + 1
  
    frj = -1
    while frj < nrch_g-1:
        frj = frj + 1  
        if rpynw[frj] in head_rlink_multi_trib_links.keys():
            # more than one tributary link enter to a single mainstem rlink
            ntrib_links = head_rlink_frequency[rpynw[frj]] # number of tributary links entering to this mainstem rlink.   
            
            rfrj = frj
            niter = ntrib_links-1
            temp_ntrib_links = ntrib_links
            for iter in range(0, niter):
                # step1: label it invalid reach by the dummy value
                rfrnw_g[rfrj, :] = -33 
                rpynw[rfrj] = -33
                
                # step2: fix downstream frjs of upstream reaches both in mainstem and tributary
                if rfrnw_g[rfrj-1, 1] != -33:
                    rfrnw_g[rfrj-1, 1] = rfrnw_g[rfrj-1, 1] + 2*(temp_ntrib_links - 1)
                if rfrnw_g[rfrj-2, 1] != -33:
                    rfrnw_g[rfrj-2, 1] = rfrnw_g[rfrj-2, 1] + 2*(temp_ntrib_links - 1)
                    
                rfrj = rfrj + 2
                temp_ntrib_links = temp_ntrib_links - 1
                
            # step3:  fix values of frj of final mainstem rlink  
            rfrnw_g[rfrj, 2] = 1 + ntrib_links # number of upstream reaches (1 for upstream mainstem rlink)
            frj_start = rfrj - 2* ntrib_links 
            frj_end   = rfrj - 1
            i=2
            for frj2 in range(frj_start, frj_end + 1):
                if rfrnw_g[frj2, 1] != -33:
                    # frj of upstream reaches
                    i = i + 1
                    rfrnw_g[rfrj, i] = frj2 + 1 # +1 becasue python row index + 1 = fortran row index
            frj = rfrj 

    return rfrnw_g, rpynw, rordered_reaches 

def fp_refactored_qlat_iniq_dx_map(
    mxncomp_g,
    nrch_g,
    mx_jorder,
    ordered_reaches,                                                     
    frnw_g,
    nts_ql_g, 
    qlat_g,   
    initial_conditions,   
    param_df,
    rfrnw_g, 
    rpynw,
    rordered_reaches,            
    refactored_diffusive_domain,
    ): 
                                           
    """
    Data assimilatoin data mapping between Python and Fortran
    
    Parameters
    ----------
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    mx_jorder -- (int) max junction order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    frnw_g -- (nparray of int) Fortran-Python network mapping array 
    nts_ql_g -- (int) numer of qlateral timesteps
    qlat_g -- (ndarray of float32) empty qlateral array to be filled    
    initial_conditions -- (ndarray of float32) initial flow (m3/sec) and depth (m above ch bottom) states for network nodes
    param_df --(DataFrame) geomorphic parameters
    rfrnw_g -- (nparray of int) Fortran-Python network mapping array on refactored hydrofabric
    rpynw -- (dict) ordered reach head segments on refactored hydrofabric
    rordered_reaches -- (dict) reaches and reach metadata by junction order on refactored hydrofabric
    refactored_diffusive_domain -- (dict) geometric relationship information between original and refactored hydrofabrics

    Returns
    -------
    rqlat_g -- (ndarray of float32) empty qlateral array to be filled on refactored hydrofabric
    riniq -- (ndarray of float32) initial flow (m3/sec) for network nodes on refactored hydrofabric
    rdx_ar_g -- (numpy of float64s) segment length on refactored hydrofabric (meters)
    linkid_mainstem_ordered -- (list) segment IDs of mainstem only (excluding tributaries) by junction order
    linkid_mainstem_indices -- (dict) mainstem segment ID : [segment index,  reach index]
    linkid_indices -- (dict) segment ID of both mainstem and tributaries : [segment index,  reach index]    
    """
    linkid_indices={} 
    linkid_mainstem_indices = {}
    rdx_ar_g = np.zeros((mxncomp_g, nrch_g))
    linkid_mainstem_ordered=[]
    
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp-1):                
                linkid = seg_list[seg]
                linkid_indices[linkid] = [seg, frj]
                if  frnw_g[frj, 2] > 0:
                    linkid_mainstem_indices[linkid] = [seg, frj]                
                    linkid_mainstem_ordered.append(linkid)    

    rlinkid_indices= {}
    for frj in range(0, nrch_g):
        if  rfrnw_g[frj, 2] > 0:
            rreach = rordered_reaches[rpynw[frj]]
            for seg in range(0,len(rreach)):
                rlinkid_indices[rreach[seg]] = [seg, frj]

    rqlat_g = np.zeros((nts_ql_g, mxncomp_g, nrch_g))
    riniq   = np.zeros((mxncomp_g, nrch_g))

    for frj in range(0, nrch_g):
        if  rfrnw_g[frj, 2] > 0:
            # mainstem only
            rncomp = rfrnw_g[frj, 0]
            for rseg in range(0, rncomp-1): #qlat at bottom node is not used in diffusive.f90 so set it zero.
                rlinkid = rordered_reaches[rpynw[frj]][rseg]
                linkid_list = refactored_diffusive_domain['lengthMap'][rlinkid].split(',')
                sumql_link      = 0.0
                sumq0_link      = 0.0
                sumlength_rlink = 0.0
                for seg in range(len(linkid_list)):
                    linkid_dxfrac = float(linkid_list[seg])
                    intpart,decimalpart = int(linkid_dxfrac),linkid_dxfrac-int(linkid_dxfrac)
                    if intpart in linkid_indices:
                        # sometimes rlink goes further downstream of TW link so as to include more links that are not in diffusive domain.
                        i_index = linkid_indices[intpart][0]
                        j_index = linkid_indices[intpart][1]
                        qlat_link = qlat_g[:, i_index, j_index] # [m^2/sec]   
                        q0_link   = initial_conditions.loc[intpart,'qu0']  
                        dx_link = param_df.loc[intpart, 'dx']
                        
                        sumql_link      = sumql_link + qlat_link*dx_link*decimalpart*10.0 #[m^3/s] 
                        sumlength_rlink = sumlength_rlink + dx_link*decimalpart*10.0
                        sumq0_link      = sumq0_link + q0_link*dx_link*decimalpart*10.0   #[m^4/s]                
                
                rqlat_g[:,rseg,frj] = sumql_link / sumlength_rlink # [m^2/sec]
                rdx_ar_g[rseg, frj] = sumlength_rlink
                riniq[rseg, frj]    = sumq0_link / sumlength_rlink # weighted averge by a weight of dx_link / dx_rlink
                if rseg == rncomp-2:
                    riniq[rseg+1, frj] = riniq[rseg, frj]
    
    return rqlat_g, riniq, rdx_ar_g, linkid_indices,  linkid_mainstem_indices, linkid_mainstem_ordered
    
def fp_refactored_naturalxsec_map(
    tw,
    mxncomp_g,
    nrch_g,                                                    
    topobathy_bytw, 
    param_df,
    rfrnw_g, 
    rpynw,
    rordered_reaches,            
    refactored_diffusive_domain,
    ):  
            
    """
    Data assimilatoin data mapping between Python and Fortran
    
    Parameters
    ----------
    tw,
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    topobathy_bytw --(DataFrame) natural channel cross section's x and z values with manning's N
    param_df --(DataFrame) geomorphic parameters
    rfrnw_g -- (nparray of int) Fortran-Python network mapping array on refactored hydrofabric
    rpynw -- (dict) ordered reach head segments on refactored hydrofabric
    rordered_reaches -- (dict) reaches and reach metadata by junction order on refactored hydrofabric
    refactored_diffusive_domain -- (dict) geometric relationship information between original and refactored hydrofabrics

    Returns
    -------
    x_bathy_g -- (numpy of float64s) lateral distance of bathy data points on refactored hydrofabric
    z_bathy_g -- (numpy of float64s) elevation of bathy data points on refactored hydrofabric
    mann_bathy_g -- (numpy of float64s) Manning's N value of bathy data points on refactored hydrofabric
    size_bathy_g -- (integer) the nubmer of bathy data points of each cross section on refactored hydrofabric
    mxnbathy_g -- (integer) maximum size of bathy data points on refactored hydrofabric
    """    
    
    if not topobathy_bytw.empty:
        # maximum number of stations along a single cross section
        mxnbathy_g = topobathy_bytw.index.value_counts().max()

        # initialize arrays to store cross section data
        x_bathy_g    = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
        z_bathy_g    = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
        mann_bathy_g = np.zeros((mxnbathy_g, mxncomp_g, nrch_g))
        size_bathy_g = np.zeros((mxncomp_g, nrch_g), dtype='i4')   
    
        rlinkid_tw = refactored_diffusive_domain['refac_tw']

        for frj in range(0, nrch_g):
            if  rfrnw_g[frj, 2] > 0:
                # mainstem only
                rncomp = rfrnw_g[frj, 0]
                rlink_list = rordered_reaches[rpynw[frj]]

                for rseg in range(0, rncomp):
                    if rseg < rncomp-1:
                        segid = rlink_list[rseg]
                    else:
                        segid = rlink_list[rseg-1]
        
                    if rseg == rncomp-1:
                        if segid != rlinkid_tw:
                            # if last node of a reach, but not the last node in TW reach
                            # use cross section of head segment of the downstream reach
                            i = frj + 1
                            while rfrnw_g[i, 2] <= 0:
                                i = i + 1                        
                            rlinkid = rordered_reaches[rpynw[i]][0]
                        else:                            
                            # if last node of reach AND last node in the network,
                            # use cross section of upstream rlink to interpolate toward bottom node of TW reach
                            rlinkid = rlink_list[rncomp-2]                            
                    else:
                        rlinkid = rlink_list[rseg]
                    
                    # how many topo data points in a given cross section?
                    nstations = len(topobathy_bytw.loc[rlinkid])
                        
                    # populate cross section size (# of stations) array
                    size_bathy_g[rseg, frj] = nstations
                        
                    # populate cross section x, z and mannings n arrays
                    x_bathy_g[0:nstations, rseg, frj]    = topobathy_bytw.loc[rlinkid].xid_d
                    z_bathy_g[0:nstations, rseg, frj]    = topobathy_bytw.loc[rlinkid].z
                    mann_bathy_g[0:nstations, rseg, frj] = topobathy_bytw.loc[rlinkid].n
                        
                    # if terminal node of the network, then adjust the cross section z data using
                    # channel slope and length data
                    if rseg == rncomp-1 and rlinkid == rlinkid_tw:
                        So = param_df.loc[tw].s0
                        dx = param_df.loc[tw].dx                    
                        z_bathy_g[0:nstations, rseg, frj] = z_bathy_g[0:nstations, rseg, frj] - So * dx  
    else: 
        #if the bathy dataframe is empty, then pass out empty arrays
        x_bathy_g    = np.array([]).reshape(0,0,0)
        z_bathy_g    = np.array([]).reshape(0,0,0)
        mann_bathy_g = np.array([]).reshape(0,0,0)
        size_bathy_g = np.array([], dtype = 'i4').reshape(0,0)
        mxnbathy_g   = int(0)
   
    return mxnbathy_g, x_bathy_g, z_bathy_g, mann_bathy_g, size_bathy_g
            
def fp_crosswalk_map(
    nrch_g,                                                    
    rfrnw_g, 
    rpynw,
    rordered_reaches,            
    linkid_mainstem_ordered,
    linkid_mainstem_indices,
    linkid_indices,
    refactored_diffusive_domain
    ):  
            
    """
    Data assimilatoin data mapping between Python and Fortran
    
    Parameters
    ----------
    nrch_g -- (int) number of reaches in the network
    rfrnw_g -- (nparray of int) Fortran-Python network mapping array on refactored hydrofabric
    rpynw -- (dict) ordered reach head segments on refactored hydrofabric
    rordered_reaches -- (dict) reaches and reach metadata by junction order on refactored hydrofabric
    linkid_mainstem_ordered -- (list) segment IDs of mainstem only (excluding tributaries) by junction order
    linkid_mainstem_indices -- (dict) mainstem segment ID : [segment index,  reach index]
    linkid_indices -- (dict) segment ID of both mainstem and tributaries : [segment index,  reach index]
    refactored_diffusive_domain -- (dict) geometric relationship information between original and refactored hydrofabrics

    Returns
    -------
    crosswalk_nrow -- (int) number of rows in the crosswalk lookup table
    crosswalk_ncol -- (int) number of columns in the crosswalk lookup table
    crosswalk_g -- (numpy of float64s) crosswalk information between original and refactored hydrofabric, used for mapping 
                                       routing results back to original hydrofabric
    """
    crosswalk_nrow = len(refactored_diffusive_domain['lengthMap'].keys())
    crosswalk_ncol = 15 + 3
    crosswalk_g    = np.zeros((crosswalk_nrow, crosswalk_ncol))

    irow = -1    
    for frj in range(0, nrch_g):
        if  rfrnw_g[frj, 2] > 0:
            # mainstem only
            rncomp = rfrnw_g[frj, 0]
            for rseg in range(0, rncomp-1): #qlat at bottom node is not used in diffusive.f90 so set it zero.
                irow = irow + 1
                rlinkid = rordered_reaches[rpynw[frj]][rseg]
                crosswalk_g[irow, 0] = rseg + 1 # fortran index of [rseg,    ] for refactored 
                crosswalk_g[irow, 1] = frj + 1  # fortran index of [    , frj] for refactored                 
                
                linkid_frac_list = refactored_diffusive_domain['lengthMap'][rlinkid].split(',')         
                linkid_list=[]
                linkid_frac_float_list=[]
                for seg in range(len(linkid_frac_list)):
                    linkid_frac = float(linkid_frac_list[seg])
                    linkid,dxfrac = int(linkid_frac),linkid_frac-int(linkid_frac)
                    linkid_list.append(linkid) 
                    linkid_frac_float_list.append(linkid_frac)

                # order linkid_frac list from upstream to downstream
                linkid_ordered = [x for x in linkid_mainstem_ordered if x in linkid_list] 
                linkid_frac_float_ordered =[]
                for seg in range(len(linkid_frac_list)):
                    linkid_frac = linkid_frac_float_list[seg]
                    if int(linkid_frac) in linkid_indices:     
                        # sometimes rlink goes further downstream of TW link so as to include more links that are not in diffusive domain.
                        idx = linkid_ordered.index(int(linkid_frac))
                        linkid_frac_float_ordered.insert(idx, linkid_frac)                
                crosswalk_g[irow, 2] = len(linkid_frac_float_ordered)  # number of links to cover a rlink     
                
                jcol = 3
                for seg in range(len(linkid_frac_float_ordered)):                     
                    linkid_frac = linkid_frac_float_ordered[seg]
                    linkid, frac = int(linkid_frac),linkid_frac-int(linkid_frac)
                    if linkid in linkid_indices:     
                        # sometimes rlink goes further downstream of TW link so as to include more links that are not in diffusive domain.
                        i_index = linkid_mainstem_indices[linkid][0]
                        j_index = linkid_mainstem_indices[linkid][1]
                                 
                        crosswalk_g[irow, jcol]   = i_index + 1  # +1 for fortran indexing
                        crosswalk_g[irow, jcol+1] = j_index + 1
                        crosswalk_g[irow, jcol+2] = frac*10.0                        
                        jcol = jcol + 3   
    
    return crosswalk_nrow, crosswalk_ncol, crosswalk_g 

def fp_coastal_boundary_input_map(
    tw,
    coastal_boundary_depth_df, 
    nrch_g,
    t0,
    t0_g,
    tfin_g):
    """
    Data assimilatoin data mapping between Python and Fortran
    Parameters
    ----------
    tw -- (int) Tailwater segment ID
    coastal_boundary_depth_df -- (DataFrame) coastal boundary depth data at one hour time steps
    t0              -- (datetime) initial date
    t0_g            -- (float) diffusive model's initial simulation time [hr] (by default, zero)
    tfin_g          -- (float) diffusive model's final simulation time [hr] 
    
    Returns
    -------
    dt_db_g         -- (int)   coastal boundary input data timestep in sec
    dsbd_option     -- (int) 1 or 2 for coastal boundary depth data or normal depth data, respectively
    nts_db_g        -- (int) number of coastal boundary input data timesteps
    dbcd_g          -- (float) coastal boundary input data time series [m]
    """    
  
    if not coastal_boundary_depth_df.empty:   
        date_time_obj1 = datetime.strptime(coastal_boundary_depth_df.columns[1], '%Y-%m-%d %H:%M:%S')
        date_time_obj0 = datetime.strptime(coastal_boundary_depth_df.columns[0], '%Y-%m-%d %H:%M:%S')
        dt_db_g        = (date_time_obj1 - date_time_obj0).total_seconds()    
        nts_db_g       = int((tfin_g - t0_g) * 3600.0 / dt_db_g) + 1  # include initial time 0 to the final time    
        dbcd_g         = np.ones(nts_db_g)

        dt_timeslice = timedelta(minutes=dt_db_g/60.0)
        tfin         = t0 + dt_timeslice*(nts_db_g-1)
        timestamps   = pd.date_range(t0, tfin, freq=dt_timeslice)
 
        timeslice_dbcd_list = []
        tws = coastal_boundary_depth_df.index.values.tolist()
        for i in range(len(timestamps)):
            timeslice = np.full(len(tws), timestamps[i])
            if str(timestamps[i]) in coastal_boundary_depth_df.columns:
                depth = coastal_boundary_depth_df[str(timestamps[i])].tolist()
            else: 
                depth = np.nan

            timeslice_dbcd = (pd.DataFrame({
                                            'stationId' : tws,
                                            'datetime'  : timeslice,
                                            'depth'     : depth 
                                            }).set_index(['stationId', 'datetime']).
                                                unstack(1, fill_value = np.nan)['depth'])
            
            timeslice_dbcd_list.append(timeslice_dbcd)                         

        dbcd_df = pd.concat(timeslice_dbcd_list, axis=1, ignore_index=False)
        
        # replace zero or negative water depth by the smallest postive value
        psmindepth =  dbcd_df.where(dbcd_df>0).min(1).loc[tw]
        dbcd_df.loc[tw, dbcd_df.loc[tw]<= 0] = psmindepth        
        
        # interpolate missing value in NaN (np.nan) inbetween available values. For extrapolation, the result is the same as 
        # a last available value either in the left or right.  
        dbcd_df_interpolated = dbcd_df.interpolate(axis='columns', limit_direction='both', limit=6)

        # if still missing data exists, do not use the coastal depth data as diffusive downstream boundary condition
        if dbcd_df_interpolated.isnull().values.any():
            dsbd_option = 2 # instead, use normal depth as the downstream boundary condition
            dbcd_g[:] = 0.0
        else: 
            dsbd_option = 1 # use the data as prepared
            dbcd_g[:] = dbcd_df_interpolated.loc[tw].values
            
    else:
        dt_db_g     = 3600.0 # by default in sec
        nts_db_g    = int((tfin_g - t0_g) * 3600.0 / dt_db_g) + 1 # include initial time 0 to the final time
        dbcd_g      = np.ones(nts_db_g)
        dsbd_option = 2 # instead, use normal depth as the downstream boundary condition
        dbcd_g[:]   = 0.0               

    return dt_db_g, dsbd_option, nts_db_g, dbcd_g

def fp_thalweg_elev_map(
                ordered_reaches,
                mainstem_seg_list, 
                unrefactored_topobathy_bytw,
                param_df, 
                mx_jorder, 
                mxncomp_g, 
                nrch_g,
                dbfksegID):
    """
    thalweg elevation profile of natural cross sections on unrefactored hydrofabric for crosswalking b/t refactored and unrefactored hydrofabrics
    
    Parameters
    ----------
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    mainstem_seg_list -- (int) a list of link IDs of segs of related mainstem reaches 
    unrefactored_topobathy_bytw -- (DataFrame) natural cross section's x and z values with manning's N on unrefactored hydrofabric   
    param_df --(DataFrame) geomorphic parameters
    mx_jorder -- (int) max junction order
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    dbfksegID -- (int) segment ID of fake node (=bottom node) of TW reach that hosts downstream boundary 
                        condition for a network that is being routed. 
    
    Returns
    -------
    z_thalweg_g -- (numpy of float64s) channel bottom altitude on unrefactored channel networks
    
    Notes
    -----
    - In node-configuration, bottom node takes bathy of the first segment of the downtream reach
      except TW reach where the bottom node bathy is interpolated by bathy of the last segment 
      with so*0.5*dx 
    """  
    if not unrefactored_topobathy_bytw.empty:        
        # initialize arrays to store thalweg elevation values
        z_thalweg_g  = np.zeros((mxncomp_g, nrch_g)) 
        
        # loop over reach orders.
        frj = -1
        for x in range(mx_jorder, -1, -1):
            # loop through all reaches of order x
            for head_segment, reach in ordered_reaches[x]:
                frj = frj + 1
                # list of segments in this reach
                seg_list = reach["segments_list"]
                # number of segments in this reach
                ncomp = reach["number_segments"]

                # determine if this reach is part of the mainstem diffusive domain
                if head_segment in mainstem_seg_list: 
                    # loop through segments in mainstem reach
                    for seg, segID in enumerate(seg_list):
                        # identify the index in topobathy dataframe that contains
                        # the data we want for this node.
                        if seg == ncomp-1 and x > 0:                             
                            # if last node of a reach, but not the last node in the network
                            # use cross section of downstream neighbor
                            seg_idx = reach["downstream_head_segment"][0]
                        
                        elif segID == dbfksegID:                            
                            # if last node of reach AND last node in the network,
                            # use cross section of upstream neighbor
                            seg_idx = seg_list[seg-1]
                            
                        else:
                            seg_idx = segID

                        # find channel bottom's lowest elevation value
                        z_thalweg_g[seg, frj] = unrefactored_topobathy_bytw.loc[seg_idx].z.min()
              
                        # if terminal node of the network, then adjust the thalweg elevation using
                        # channel slope and length data
                        if segID == dbfksegID:
                            So = param_df.loc[seg_idx].s0
                            dx = param_df.loc[seg_idx].dx
                            z_thalweg_g[seg, frj] = z_thalweg_g[seg, frj] - So * dx  

    else: 
        #if the bathy dataframe is empty, then pass out empty arrays
        z_thalweg_g  = np.array([]).reshape(0,0)
    
    return z_thalweg_g

def fp_thalweg_elev_routelink_map(
                                  mx_jorder, 
                                  ordered_reaches, 
                                  z_all, 
                                  mxncomp_g, 
                                  nrch_g,
                                  ):
    """
    thalweg elevation profile of synthetic cross sections of RouteLink.nc for crosswalking b/t refactored and RouteLink hydrofabrics
    
    Parameters
    ----------
    mx_jorder -- (int) maximum network reach order
    ordered_reaches -- (dict) reaches and reach metadata by junction order
    z_all -- (dict) adjusted altitude dictionary
    mxncomp_g -- (int) maximum number of nodes in a reach
    nrch_g -- (int) number of reaches in the network
    
    Returns
    -------
    z_thalweg_g -- (numpy of float64s) channel bottom altitude (meters) from RouteLink.nc 
    """
    z_thalweg_g = np.zeros((mxncomp_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                segID = seg_list[seg]
                z_thalweg_g[seg, frj] = z_all[segID]["adj.alt"][0]

    return z_thalweg_g
     
def diffusive_input_data_v02(
    tw,
    connections,
    rconn,
    reach_list,
    mainstem_seg_list,
    trib_seg_list,
    diffusive_parameters,
    param_df,
    qlat,
    initial_conditions,
    junction_inflows,
    qts_subdivisions,
    t0,
    nsteps,
    dt,
    waterbodies_df,
    topobathy_bytw,
    usgs_df,
    refactored_diffusive_domain,
    refactored_reaches,
    coastal_boundary_depth_df,
    unrefactored_topobathy_bytw,
):

    """
    Build input data objects for diffusive wave model
    
    Parameters
    ----------
    tw -- (int) Tailwater segment ID
    connections -- (dict) donwstream connections for each segment in the network
    rconn -- (dict) upstream connections for each segment in the network
    reach_list -- (list of lists) lists of segments comprising different reaches in the network
    diffusive_parametters -- (dict) Diffusive wave model parameters
    param_df --(DataFrame) geomorphic parameters
    qlat_data -- (ndarray of float32) qlateral data (m3/sec)
    initial_conditions -- (ndarray of float32) initial flow (m3/sec) and depth (m above ch bottom) states for network nodes
    upstream_results -- (dict) with values of 1d arrays upstream flow, velocity, and depth   
    qts_subdivisions -- (int) number of qlateral timestep subdivisions    
    topobathy_bytw --(DataFrame) natural channel cross section data of a channel network draining into a tailwater node
    usgs_df --(DataFrame) observed usgs flow data
    refactored_diffusive_domain -- (dict) geometric relationship information between original and refactored hydrofabrics
    refactored_reaches -- (list of lists) lists of stream segment IDs of diffusive mainstems on refactored hydrofabrics including 
                                          tributaries of original hydrofabric
 
    Returns
    -------
    diff_ins -- (dict) formatted inputs for diffusive wave model
    """

    # lateral inflow timestep (sec). Currently, lateral flow in CHAROUT files at one hour time step
    dt_ql_g = 3600.0
    # upstream boundary condition timestep = MC simulation time step (sec)
    dt_ub_g = 300.0 
    # tributary inflow timestep (sec) - currently, MC simulation time step = 5 min
    dt_qtrib_g = 300.0
    # usgs_df time step used for data assimilation. The original timestep of USGS streamflow is 15 min but later interpolated at every dt in min.
    dt_da_g = dt # [sec]
    # time interval at which flow and depth simulations are written out by Tulane diffusive model
    saveinterval_tu = dt
    # time interval at which depth is written out by cnt model
    #saveinterval_cnt = dt * (qts_subdivisions)
    # time interval at which flow is computed and written out by cnt model
    # initial timestep interval used by Tulane diffusive model
    dtini_g = dt
    t0_g = 0.0  # simulation start hr **set to zero for Fortran computation
    tfin_g = (dt * nsteps)/60/60
    
    # package timestep variables into single array
    timestep_ar_g    = np.zeros(9)
    timestep_ar_g[0] = dtini_g
    timestep_ar_g[1] = t0_g
    timestep_ar_g[2] = tfin_g
    timestep_ar_g[3] = saveinterval_tu 
    timestep_ar_g[4] = dt_ql_g
    timestep_ar_g[5] = dt_ub_g
    # timestep_ar_g[6] =  dt_db_g <- defined after calling fp_coastal_boundary_input_map
    timestep_ar_g[7] = dt_qtrib_g
    timestep_ar_g[8] = dt_da_g

    # CN-mod parameters
    paradim       = 11
    para_ar_g     = np.zeros(paradim)
    para_ar_g[0]  = 0.95    # Courant number (default: 0.95)
    para_ar_g[1]  = 0.5     # lower limit of celerity (default: 0.5)
    para_ar_g[2]  = 10.0    # lower limit of diffusivity (default: 10)
    para_ar_g[3]  = 10000.0  # upper limit of diffusivity (default: 10000)
    para_ar_g[4]  = -15.0   # lower limit of dimensionless diffusivity, used to determine b/t normal depth and diffusive depth
    para_ar_g[5]  = -10.0   #upper limit of dimensionless diffusivity, used to determine b/t normal depth and diffusive depth
    para_ar_g[6]  = 1.0     # 0:run Bisection to compute water level; 1: Newton Raphson (default: 1.0)
    para_ar_g[7]  = 0.02831   # lower limit of discharge (default: 0.02831 cms)
    para_ar_g[8]  = 0.0001    # lower limit of channel bed slope (default: 0.0001)
    para_ar_g[9]  = 1.0     # weight in numerically computing 2nd derivative: 0: explicit, 1: implicit (default: 1.0)
    para_ar_g[10] = 2      # downstream water depth boundary condition: 1: given water depth data, 2: normal depth
    # number of reaches in network
    nrch_g = len(reach_list)

    # maximum number of nodes in a reach
    mxncomp_g = 0
    for r in reach_list:
        nnodes = len(r) + 1
        if nnodes > mxncomp_g:
            mxncomp_g = nnodes

# TODO: How do we plan to utilize upstream boundary condition data object?
#     ds_seg = []
#     offnet_segs = []
#     upstream_flow_array = np.zeros((len(ds_seg), nsteps+1))
#     if upstream_results:

#         # create a list of segments downstream of offnetwork upstreams [ds_seg]
#         # and a list of offnetwork upstream segments [offnet_segs]
#         inv_map = nhd_network.reverse_network(rconn)
#         for seg in upstream_results:
#             ds_seg.append(inv_map[seg][0])
#             offnet_segs.append(seg)

#         # populate an array of upstream flows (boundary condtions)
#         upstream_flow_array = np.zeros((len(set(ds_seg)), nsteps+1))
#         for j, seg in enumerate(set(ds_seg)):

#             # offnetwork-upstream connections
#             us_segs = rconn[seg]

#             # sum upstream flows and initial conditions
#             usq = np.zeros((len(us_segs), nsteps))
#             us_iniq = 0
#             for k, s in enumerate(us_segs):
#                 usq[k] = upstream_results[s]['results'][::3]

#                 if s in lake_segs:
#                     # initial conditions from wbody_param array
#                     idx_segID = np.where(np.asarray(lake_segs) == s)
#                     us_iniq += wbody_params[idx_segID,9]
#                 else:
#                     # initial conditions from initial_conditions array
#                     idx_segID = np.where(geo_index == s)
#                     us_iniq += initial_conditions[idx_segID,0]

#             # write upstream flows to upstream_flow_array
#             upstream_flow_array[j,1:] = np.sum(usq, axis = 0)
#             upstream_flow_array[j,0] = us_iniq

    # Order reaches by junction depth
    path_func = partial(nhd_network.split_at_waterbodies_and_junctions, set(junction_inflows.index.to_list()),rconn)
    tr = nhd_network.dfs_decomposition_depth_tuple(rconn, path_func)    

    jorder_reaches = sorted(tr, key=lambda x: x[0])
    mx_jorder = max(jorder_reaches)[0]  # maximum junction order of subnetwork of TW

    ordered_reaches = {}
    rchhead_reaches = {}
    rchbottom_reaches = {}
    z_all = {}
    for o, rch in jorder_reaches:

        # add one more segment(fake) to the end of a list of segments to account for node configuration.
        fksegID = int(str(rch[-1]) + str(2))
        rch.append(fksegID)

        # additional segment(fake) to upstream bottom segments
        fk_usbseg = [int(str(x) + str(2)) for x in rconn[rch[0]]] 

        if o not in ordered_reaches:
            ordered_reaches.update({o: []})

        # populate the ordered_reaches dictionary with node connection information
        ordered_reaches[o].append(
            [
                rch[0],
                {
                    "number_segments": len(rch),
                    "segments_list": rch,
                    "upstream_bottom_segments": fk_usbseg,
                    "downstream_head_segment": connections[rch[-2]],
                },
            ]
        )

        if rch[0] not in rchhead_reaches:

            # a list of segments for a given reach-head segment
            rchhead_reaches.update(
                {rch[0]: {"number_segments": len(rch), "segments_list": rch}}
            )
            # a list of segments for a given reach-bottom segment
            rchbottom_reaches.update(
                {rch[-1]: {"number_segments": len(rch), "segments_list": rch}}
            )

        # for channel altitude adjustment
        z_all.update({seg: {"adj.alt": np.zeros(1)} for seg in rch})

    # --------------------------------------------------------------------------------------
    #                                 Step 0-3
    #    Adjust altitude so that altitude of the last sement of a reach is equal to that
    #    of the first segment of its downstream reach right after their common junction.
    # --------------------------------------------------------------------------------------
    dbfksegID = int(str(tw) + str(2))

    adj_alt1(
        mx_jorder, ordered_reaches, param_df, dbfksegID, z_all
    )

    # --------------------------------------------------------------------------------------
    #                                 Step 0-4
    #     Make Fortran-Python channel network mapping variables.
    # --------------------------------------------------------------------------------------
    # build a list of head segments in descending reach order [headwater -> tailwater]
    pynw = {}
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            frj = frj + 1
            pynw[frj] = head_segment

    frnw_col = 15
    frnw_g   = fp_network_map(
                              mainstem_seg_list, 
                              mx_jorder, 
                              ordered_reaches, 
                              rchbottom_reaches, 
                              nrch_g, 
                              frnw_col, 
                              dbfksegID, 
                              pynw,
                              #upstream_boundary_link,
                              )

    # covert data type from integer to float for frnw  
    dfrnw_g = frnw_g.astype('float')    
    # ---------------------------------------------------------------------------------
    #                              Step 0-5
    #                  Prepare channel geometry data
    # ---------------------------------------------------------------------------------
    (
        z_ar_g,
        bo_ar_g,
        traps_ar_g,
        tw_ar_g,
        twcc_ar_g,
        mann_ar_g,
        manncc_ar_g,
        so_ar_g,
        dx_ar_g,
    ) = fp_chgeo_map(
        mx_jorder,
        ordered_reaches,
        param_df,
        z_all,
        mxncomp_g,
        nrch_g,
    )

    # ---------------------------------------------------------------------------------
    #                              Step 0-6
    #                  Prepare initial conditions data
    # ---------------------------------------------------------------------------------
    iniq = np.zeros((mxncomp_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            seg_list = reach["segments_list"]
            ncomp = reach["number_segments"]
            frj = frj + 1
            for seg in range(0, ncomp):
                if seg == ncomp - 1:
                    segID = seg_list[seg - 1]
                else:
                    segID = seg_list[seg]

                # retrieve initial condition from initial_conditions DataFrame
                iniq[seg, frj] = initial_conditions.loc[segID, 'qu0']

                # set lower limit on initial flow condition
                if iniq[seg, frj]<0.0001:
                    iniq[seg, frj]=0.0001

    # ---------------------------------------------------------------------------------
    #                              Step 0-7

    #                  Prepare lateral inflow data
    # ---------------------------------------------------------------------------------
    nts_ql_g = (
        math.ceil((tfin_g - t0_g) * 3600.0 / dt_ql_g)
    )  # the number of the entire time steps of lateral flow data

    qlat_g = np.zeros((nts_ql_g, mxncomp_g, nrch_g))

    fp_qlat_map(
        mx_jorder,
        ordered_reaches,
        nts_ql_g,
        param_df,
        qlat,
        qlat_g,
    )
    
    # ---------------------------------------------------------------------------------
    #                              Step 0-8

    #       Prepare upstream boundary (top segments of head basin reaches) data
    # ---------------------------------------------------------------------------------
    # currently all zeroes by default.  
    nts_ub_g = int((tfin_g - t0_g) * 3600.0 / dt_ub_g)  
    ubcd_g = np.zeros((nts_ub_g, nrch_g)) 
    
    # ---------------------------------------------------------------------------------
    #                              Step 0-9

    #       Prepare downstrea boundary (bottom segments of TW reaches) data
    # ---------------------------------------------------------------------------------   
    dt_db_g, dsbd_option, nts_db_g, dbcd_g =  fp_coastal_boundary_input_map(
                                                                        tw,
                                                                        coastal_boundary_depth_df, 
                                                                        nrch_g,
                                                                        t0,
                                                                        t0_g,
                                                                        tfin_g)
    
    timestep_ar_g[6] = dt_db_g
    para_ar_g[10] = dsbd_option  # downstream water depth boundary condition: 1: given water depth data, 2: normal depth    
   
    # ---------------------------------------------------------------------------------------------
    #                              Step 0-9-2

    #       Prepare tributary q time series data generated by MC that flow into a juction boundary  
    # ---------------------------------------------------------------------------------------------
    nts_qtrib_g = int((tfin_g - t0_g) * 3600.0 / dt_qtrib_g) + 1 # Even MC-computed flow start from first 5 min, t0 is coverd by initial_conditions. 
    qtrib_g = np.zeros((nts_qtrib_g, nrch_g))
    frj = -1
    for x in range(mx_jorder, -1, -1):
        for head_segment, reach in ordered_reaches[x]:
            frj = frj + 1
            if head_segment not in mainstem_seg_list:
                qtrib_g[1:,frj] = junction_inflows.loc[head_segment]                
                # TODO - if one of the tributary segments is a waterbody, it's initial conditions
                # will not be in the initial_conditions array, but rather will be in the waterbodies_df array
                qtrib_g[0,frj] = initial_conditions.loc[head_segment, 'qu0']    
         
    # ---------------------------------------------------------------------------------
    #                              Step 0-10

    #                 Prepare cross section bathymetry data
    # ---------------------------------------------------------------------------------    
    if not refactored_diffusive_domain:
        x_bathy_g, z_bathy_g, mann_bathy_g, size_bathy_g, mxnbathy_g = fp_naturalxsec_map(        
                                                                           ordered_reaches,                             
                                                                           mainstem_seg_list, 
                                                                           topobathy_bytw,
                                                                           param_df, 
                                                                           mx_jorder,
                                                                           mxncomp_g, 
                                                                           nrch_g,
                                                                           dbfksegID)
                   
    # ---------------------------------------------------------------------------------------------
    #                              Step 0-11

    #       Prepare interpolated USGS streamflow values at every dt_da_g time step [sec]   
    # ---------------------------------------------------------------------------------------------
    nts_da_g, usgs_da_g, usgs_da_reach_g = fp_da_map(
                                                mx_jorder,
                                                ordered_reaches,
                                                usgs_df,
                                                nrch_g,
                                                t0,
                                                nsteps,
                                                dt_da_g,
                                                t0_g,
                                                tfin_g)
    
    # ---------------------------------------------------------------------------------------------
    #                              Step 0-12-1

    #              Prepare python-fortarn map for refactored hydrofabric   
    #
    #  link : stream segment of unrefactored hydrofabric
    #  rlink: steream segment of refactored hydrofabric (refactored mainly to 
    #                                                    increase min.lenght of segment)
    # ---------------------------------------------------------------------------------------------    
    if refactored_diffusive_domain:
        rfrnw_g, rpynw, rordered_reaches = fp_refactored_network_map(
                                                                    nrch_g,
                                                                    frnw_g, 
                                                                    pynw, 
                                                                    mainstem_seg_list,
                                                                    trib_seg_list,
                                                                    refactored_diffusive_domain,
                                                                    refactored_reaches,
                                                                    )

    # ---------------------------------------------------------------------------------
    #                              Step 0-12-2

    #       Prepare qlateral, iniq, and dx for refactored hydrofabric   
    #
    #  link : stream segment of unrefactored hydrofabric
    #  rlink: steream segment of refactored hydrofabric (refactored mainly to 
    #                                                    increase min.lenght of segment)
    # ---------------------------------------------------------------------------------------------       
    if refactored_diffusive_domain:
        (
                 rqlat_g, 
                 riniq, 
                 rdx_ar_g, 
                 linkid_indices, 
                 linkid_mainstem_indices, 
                 linkid_mainstem_ordered
        ) = fp_refactored_qlat_iniq_dx_map(
                mxncomp_g,
                nrch_g,
                mx_jorder,
                ordered_reaches,
                frnw_g,
                nts_ql_g,
                qlat_g,   
                initial_conditions,   
                param_df,
                rfrnw_g, 
                rpynw,
                rordered_reaches,
                refactored_diffusive_domain,
                )

    # ---------------------------------------------------------------------------------------------
    #                              Step 0-12-3

    #       Prepare topobathy data on refactored hydrofabric   
    #
    #  link : stream segment of unrefactored hydrofabric
    #  rlink: steream segment of refactored hydrofabric (refactored mainly to 
    #                                                    increase min.lenght of segment)
    # ---------------------------------------------------------------------------------------------  
    if refactored_diffusive_domain:
        (
                 mxnbathy_g, 
                 x_bathy_g, 
                 z_bathy_g, 
                 mann_bathy_g, 
                 size_bathy_g, 
        ) = fp_refactored_naturalxsec_map(
                tw,
                mxncomp_g,
                nrch_g,                                                    
                topobathy_bytw, 
                param_df,
                rfrnw_g, 
                rpynw,
                rordered_reaches,            
                refactored_diffusive_domain,                          
                )
        # find thalweg values of natural xsections on unrefactored hydrofabric
        '''
        z_thalweg_g = fp_thalweg_elev_map(
                                          ordered_reaches,
                                          mainstem_seg_list, 
                                          unrefactored_topobathy_bytw,
                                          param_df, 
                                          mx_jorder, 
                                          mxncomp_g, 
                                          nrch_g,
                                          dbfksegID,
                                          ) 
        '''
        # find thalweg values of synthetic xsections of RouteLink.nc
        z_thalweg_g = fp_thalweg_elev_routelink_map(
                                                    mx_jorder, 
                                                    ordered_reaches, 
                                                    z_all, 
                                                    mxncomp_g, 
                                                    nrch_g,
                                                    )
    else: 
        #if refactored diffusive domain is not used, then pass out empty arrays
        z_thalweg_g  = np.array([]).reshape(0,0)

    # ---------------------------------------------------------------------------------------------
    #                              Step 0-12-3

    #       No trapezoidal data when using refactored hydrofabric except dx data from RouteLink.nc
    # ---------------------------------------------------------------------------------------------                     
    if refactored_diffusive_domain:
        rz_ar_g = np.zeros((mxncomp_g, nrch_g))
        rbo_ar_g = np.zeros((mxncomp_g, nrch_g))
        rtraps_ar_g = np.zeros((mxncomp_g, nrch_g))
        rtw_ar_g = np.zeros((mxncomp_g, nrch_g))
        rtwcc_ar_g = np.zeros((mxncomp_g, nrch_g))
        rmann_ar_g = np.zeros((mxncomp_g, nrch_g))
        rmanncc_ar_g = np.zeros((mxncomp_g, nrch_g))
        rso_ar_g = np.zeros((mxncomp_g, nrch_g))
        # rdx_ar_g has already been prepared in fp_refactored_qlat_iniq_dx_map
    else: 
        rdx_ar_g = np.array([]).reshape(0,0) 
    
    # ---------------------------------------------------------------------------------------------
    #                              Step 0-12-4

    #    Prepare crosswalk map between non-refactored and refactored hydrofabric
    # ---------------------------------------------------------------------------------------------   
    if refactored_diffusive_domain:
        (
                crosswalk_nrow, 
                crosswalk_ncol, 
                crosswalk_g 
        ) = fp_crosswalk_map(
                nrch_g,                                                    
                rfrnw_g, 
                rpynw,
                rordered_reaches,   
                linkid_mainstem_ordered,
                linkid_mainstem_indices,
                linkid_indices,
                refactored_diffusive_domain
                ) 
    else:
        crosswalk_nrow = int(0) 
        crosswalk_ncol = int(0)
        crosswalk_g  = np.array([]).reshape(0,0)

    # ---------------------------------------------------------------------------------
    #                              Step 0-13

    #                       Build input dictionary
    # ---------------------------------------------------------------------------------
    ntss_ev_g = int((tfin_g - t0_g) * 3600.0 / dt) + 1

    # build a dictionary of diffusive model inputs and helper variables
    diff_ins = {}
    if not refactored_diffusive_domain:
        # for original hydrofabric
        # model time steps
        diff_ins["timestep_ar_g"] = timestep_ar_g  
        diff_ins["nts_ql_g"]      = nts_ql_g
        diff_ins["nts_ub_g"]      = nts_ub_g
        diff_ins["nts_db_g"]      = nts_db_g
        diff_ins["nts_qtrib_g"]   = nts_qtrib_g
        diff_ins["ntss_ev_g"]     = ntss_ev_g
        diff_ins["nts_da_g"]      = nts_da_g # DA
        # max number of computation nodes of a stream reach and the number of entire stream reaches 
        diff_ins["mxncomp_g"] = mxncomp_g
        diff_ins["nrch_g"] = nrch_g
        # synthetic (trapezoidal) xsection geometry data 
        diff_ins["z_ar_g"] = z_ar_g
        diff_ins["bo_ar_g"] = bo_ar_g
        diff_ins["traps_ar_g"] = traps_ar_g
        diff_ins["tw_ar_g"] = tw_ar_g
        diff_ins["twcc_ar_g"] = twcc_ar_g
        diff_ins["mann_ar_g"] = mann_ar_g
        diff_ins["manncc_ar_g"] = manncc_ar_g
        diff_ins["so_ar_g"] = so_ar_g
        diff_ins["dx_ar_g"] = dx_ar_g
        # python-to-fortran channel network mapping
        diff_ins["frnw_col"] = frnw_col
        diff_ins["frnw_g"] = frnw_g
        # flow forcing data
        diff_ins["qlat_g"] = qlat_g
        diff_ins["ubcd_g"] = ubcd_g
        diff_ins["dbcd_g"] = dbcd_g
        diff_ins["qtrib_g"] = qtrib_g
        # diffusive model internal parameters
        diff_ins["paradim"] = paradim 
        diff_ins["para_ar_g"] = para_ar_g
        # natural xsection bathy data
        diff_ins["mxnbathy_g"] = mxnbathy_g
        diff_ins["x_bathy_g"] = x_bathy_g
        diff_ins["z_bathy_g"] = z_bathy_g
        diff_ins["mann_bathy_g"] = mann_bathy_g
        diff_ins["size_bathy_g"] = size_bathy_g    
        # initial flow value
        diff_ins["iniq"] = iniq
        # python-fortran crosswalk data
        diff_ins["pynw"] = pynw
        diff_ins["ordered_reaches"] = ordered_reaches    
        # Data Assimilation
        diff_ins["usgs_da_g"]   = usgs_da_g 
        diff_ins["usgs_da_reach_g"] = usgs_da_reach_g
        diff_ins["rdx_ar_g"] = rdx_ar_g 
        diff_ins["cwnrow_g"] = crosswalk_nrow
        diff_ins["cwncol_g"] = crosswalk_ncol
        diff_ins["crosswalk_g"] =  crosswalk_g
        diff_ins["z_thalweg_g"] = z_thalweg_g
    else:
        # for refactored hydrofabric
        # model time steps
        diff_ins["timestep_ar_g"] = timestep_ar_g  
        diff_ins["nts_ql_g"]      = nts_ql_g
        diff_ins["nts_ub_g"]      = nts_ub_g
        diff_ins["nts_db_g"]      = nts_db_g
        diff_ins["nts_qtrib_g"]   = nts_qtrib_g
        diff_ins["ntss_ev_g"]     = ntss_ev_g
        diff_ins["nts_da_g"]      = nts_da_g 
        # max number of computation nodes of a stream reach and the number of entire stream reaches 
        diff_ins["mxncomp_g"] = mxncomp_g
        diff_ins["nrch_g"] = nrch_g
        # synthetic (trapezoidal) xsection geometry data 
        diff_ins["z_ar_g"] = rz_ar_g
        diff_ins["bo_ar_g"] = rbo_ar_g
        diff_ins["traps_ar_g"] = rtraps_ar_g
        diff_ins["tw_ar_g"] = rtw_ar_g
        diff_ins["twcc_ar_g"] = rtwcc_ar_g
        diff_ins["mann_ar_g"] = rmann_ar_g
        diff_ins["manncc_ar_g"] = rmanncc_ar_g
        diff_ins["so_ar_g"] = rso_ar_g
        diff_ins["dx_ar_g"] = rdx_ar_g
        # python-to-fortran channel network mapping
        diff_ins["frnw_col"] = frnw_col
        diff_ins["frnw_g"] = rfrnw_g
        # flow forcing data
        diff_ins["qlat_g"] = rqlat_g
        diff_ins["ubcd_g"] = ubcd_g
        diff_ins["dbcd_g"] = dbcd_g
        diff_ins["qtrib_g"] = qtrib_g
        # diffusive model internal parameters
        diff_ins["paradim"] = paradim 
        diff_ins["para_ar_g"] = para_ar_g
        # natural xsection bathy data
        diff_ins["mxnbathy_g"] = mxnbathy_g
        diff_ins["x_bathy_g"] = x_bathy_g
        diff_ins["z_bathy_g"] = z_bathy_g
        diff_ins["mann_bathy_g"] = mann_bathy_g
        diff_ins["size_bathy_g"] = size_bathy_g    
        # initial flow value
        diff_ins["iniq"] = riniq
        # python-fortran crosswalk data
        diff_ins["pynw"] = pynw
        diff_ins["ordered_reaches"] = ordered_reaches    
        # Data Assimilation
        diff_ins["usgs_da_g"]   = usgs_da_g 
        diff_ins["usgs_da_reach_g"] = usgs_da_reach_g
        diff_ins["rdx_ar_g"] = rdx_ar_g 
        diff_ins["cwnrow_g"] = crosswalk_nrow
        diff_ins["cwncol_g"] = crosswalk_ncol
        diff_ins["crosswalk_g"] =  crosswalk_g   
        diff_ins["z_thalweg_g"] = z_thalweg_g
    return diff_ins

def unpack_output(pynw, ordered_reaches, out_q, out_elv):
    """
    Unpack diffusive wave output arrays
    
    Parameters
    ----------
    pynw -- (dict) ordered reach head segments
    ordered_reaches -- 
    out_q -- (diffusive._memoryviewslice of float64) diffusive wave model flow output (m3/sec)
    out_elv -- (diffusive._memoryviewslice of float64) diffusive wave model water surface elevation output (meters)
    
    Returns
    -------
    np.asarray(rch_list, dtype=np.intp) - segment indices 
    np.asarray(dat_all, dtype = 'float32') - flow, velocity, elevation array
    """

    reach_heads = list(pynw.values())
    nts = len(out_q[:, 0, 0])

    i = 1
    rch_list = []
    for o in ordered_reaches.keys():
        for rch in ordered_reaches[o]:

            rch_segs = rch[1]["segments_list"]
            rch_list.extend(rch_segs[:-1])

            j = reach_heads.index(rch[0])

            if i == 1:
                dat_all = np.empty((len(rch_segs)-1, nts * 3))
                dat_all[:] = np.nan
                # flow result
                dat_all[:, ::3] = np.transpose(np.array(out_q[:, 1 : len(rch_segs), j]))
                # elevation result
                dat_all[:, 2::3] = np.transpose(
                    np.array(out_elv[:, 1 : len(rch_segs), j])
                )

            else:
                dat_all_c = np.empty((len(rch_segs)-1, nts * 3))
                dat_all_c[:] = np.nan
                # flow result
                dat_all_c[:, ::3] = np.transpose(
                    np.array(out_q[:, 1 : len(rch_segs), j])
                )
                # elevation result
                dat_all_c[:, 2::3] = np.transpose(
                    np.array(out_elv[:, 1 : len(rch_segs), j])
                )
                # concatenate
                dat_all = np.concatenate((dat_all, dat_all_c))

            i += 1

    return np.asarray(rch_list, dtype=np.intp), np.asarray(dat_all, dtype="float32")
