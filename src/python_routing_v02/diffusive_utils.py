import os
import sys
import numpy as np
import pandas as pd
from itertools import chain, islice
import itertools
from functools import partial

import troute.nhd_network as nhd_network
#from pyuniflowtzlt import uniflow_lookuptable

def adj_alt1(mx_jorder_tw
        , ordered_reaches
        , ch_geo_data_tw 
        , dbfksegID
        , z_all   
        ):  

    # Adjust altitude so that altitude of the last node (=fake segment) of a reach is equal to that of
    # the first segment of its downstream reach right after their common junction.

    for x in range(mx_jorder_tw,-1,-1): 
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= reach['segments_list']
            ncomp= reach['number_segments']
            for seg in range(0,ncomp):
                segID= seg_list[seg]
                if seg==ncomp-1 and seg_list.count(dbfksegID)==0 :
                # At junction, the altitude of fake segment of an upstream reach
                # is equal to that of the first segment of the downstream reach
                    # head segment id of downstream reach from a junction  
                    dsrchID= reach['downstream_head_segment']                    
                    z_all[segID]['adj.alt'][0]=ch_geo_data_tw.loc[dsrchID]["alt"]
                                                      
                elif seg==ncomp-1 and seg_list.count(dbfksegID)>0:
                # Terminal downstream fakesegment                
                    ## AD HOC: need to be corrected later
                    segID2= seg_list[seg-1]
                    So= ch_geo_data_tw.loc[segID2]["s0"]
                    dx= ch_geo_data_tw.loc[segID2]["dx"]
                    z_all[segID]['adj.alt'][0]= z_all[segID2]['adj.alt'][0] - So*dx                    
                else:
                    z_all[segID]['adj.alt'][0]= ch_geo_data_tw.loc[segID]["alt"]
    
    return z_all

#-----------------------------------------------------
# Channel Network mapping between Python and Fortran
#
#-----------------------------------------------------
def fp_network_map(mx_jorder_tw
            , ordered_reaches
            , rchbottom_reaches
            , nrch_g
            , frnw_col
            , dbfksegID
            , pynw
            ):    

    #  Store headwater reach and upstream reaches above a junction
    #  as well as downstream reach after a junction
    #  into python-extension-fortran variables.     
    frnw_g=np.zeros((nrch_g,frnw_col), dtype=int)
    frj=-1
    for x in range(mx_jorder_tw,-1,-1): 
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= reach['segments_list']
            ncomp= reach['number_segments']            
            frj=frj+1
            frnw_g[frj,0]= ncomp       

            if not reach['upstream_bottom_segments']:
            # headwater reach    
                frnw_g[frj,2]= 0 # the number of upstream reaches
            else:
            # reaches before a junction
                nusrch= len(reach['upstream_bottom_segments'])   
                frnw_g[frj,2]= nusrch  # the number of upstream reaches                  
                usrch_bseg_list=list(reach['upstream_bottom_segments']) 
                i=0
                for usrch in range(0, nusrch):
                    usrch_bseg_id=usrch_bseg_list[usrch] #upstream reach's bottom segment
                    usrch_hseg_id=rchbottom_reaches[usrch_bseg_id]["segments_list"][0]                    
                    # find Fortran js corresponding to individual usrchid
                    for j, sid in pynw.items():
                        if sid == usrch_hseg_id:
                            i=i+1
                            frnw_g[frj,2+i]=j
      
            if (seg_list.count(dbfksegID)>0): 
            # a reach where downstream boundary condition is set.
                frnw_g[frj,1]= -100 # head_segment ID that is in terminal downstream reach.
                                    # That is, -100 indicates the reach of the head segment is 
                                    # terminal downstream reach where ds.bcond. happens.
            else:
            # reach after a junction                     
                dsrch_hseg_id= reach['downstream_head_segment'] 
                #fortran j index equivalent to dsrchID.  
                frnw_g[frj,1]=[j for j, sid in pynw.items() if sid==dsrch_hseg_id[0]][0] 
                
    # Adust frnw_g element values according to Fortran-Python index relationship, that is Python i = Fortran i+1
    for frj in range(0,nrch_g):
        frnw_g[frj,1]= frnw_g[frj,1]+1 #downstream reach index for frj reach
        if frnw_g[frj,2]>0:
            nusrch=frnw_g[frj,2]
            for i in range(0, nusrch):
                frnw_g[frj,3+i]= frnw_g[frj,3+i]+1 # upstream reach indicds for frj reach  
    
    
    return frnw_g

#----------------------------------------------------------
# Channel geometry data mapping between Python and Fortran
#
#----------------------------------------------------------
def fp_chgeo_map(mx_jorder_tw
            , ordered_reaches
            , ch_geo_data_tw
            , z_all
            , mxncomp_g     
            , nrch_g
            ):     
        
    z_ar_g=np.zeros((mxncomp_g, nrch_g))
    bo_ar_g=np.zeros((mxncomp_g, nrch_g))
    traps_ar_g=np.zeros((mxncomp_g, nrch_g))
    tw_ar_g=np.zeros((mxncomp_g, nrch_g))
    twcc_ar_g=np.zeros((mxncomp_g, nrch_g))
    mann_ar_g=np.zeros((mxncomp_g, nrch_g))
    manncc_ar_g=np.zeros((mxncomp_g, nrch_g))
    so_ar_g=np.zeros((mxncomp_g, nrch_g))
    dx_ar_g=np.zeros((mxncomp_g, nrch_g))
    frj=-1
    for x in range(mx_jorder_tw,-1,-1):
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= reach['segments_list']
            ncomp= reach['number_segments'] 
            frj=frj+1
            for seg in range(0,ncomp):                    
                if seg==ncomp-1:
                    segID= seg_list[seg-1]
                else:
                    segID= seg_list[seg]
                    
                bo_ar_g[seg, frj]= ch_geo_data_tw.loc[segID]["bw"]         
                traps_ar_g[seg, frj]=ch_geo_data_tw.loc[segID]["cs"]
                tw_ar_g[seg, frj]=ch_geo_data_tw.loc[segID]["tw"]
                twcc_ar_g[seg, frj]=ch_geo_data_tw.loc[segID]["twcc"]
                mann_ar_g[seg, frj]=ch_geo_data_tw.loc[segID]["n"]
                manncc_ar_g[seg, frj]=ch_geo_data_tw.loc[segID]["ncc"]
                so_ar_g[seg, frj]= ch_geo_data_tw.loc[segID]["s0"]
                dx_ar_g[seg, frj]=  ch_geo_data_tw.loc[segID]["dx"]
                segID1= seg_list[seg]  
                z_ar_g[seg, frj]= z_all[segID1]['adj.alt'][0]             
            
    return z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g 
     

#----------------------------------------------------------
# lateral inflow mapping between Python and Fortran
#
#----------------------------------------------------------
def fp_qlat_map(mx_jorder_tw
            , ordered_reaches
            , nts_ql_g
            , qlat_tw            
            , qlat_g
            ):   
    
    frj= -1
    for x in range(mx_jorder_tw,-1,-1):
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= reach['segments_list']
            ncomp= reach['number_segments']          
            frj= frj+1     
            for seg in range(0,ncomp):                               
                segID= seg_list[seg]
                for tsi in range (0,nts_ql_g):
                    if seg<ncomp-1:        
                        qlat_g[tsi,seg,frj]= qlat_tw.loc[segID][tsi]
                    else:
                        qlat_g[tsi,seg,frj]= 0.0 #seg=ncomp is actually for bottom node in Fotran code.
                                                     #And, lateral flow enters between adjacent nodes.
       
    return 

#------------------------------------------------------------------------------------------
# upstream boundary condition mapping between Python and Fortran
# 
# Assumption: top segments of headbasin reaches get their own lateral flows as 
#             upstream boundary conditions while making these qlateral zeroes
#------------------------------------------------------------------------------------------
def fp_ubcd_map(frnw_g
            , pynw
            , nts_ub_g
            , nrch_g
            , ch_geo_data_tw
            , qlat_tw
            , qlat_g
            ):   
    
    ubcd_g=np.zeros((nts_ub_g, nrch_g)) 
    frj= -1
    for frj in range(nrch_g):
        if frnw_g[frj,2]==0: # the number of upstream reaches is zero.
            head_segment=pynw[frj]
            for tsi in range(0,nts_ub_g):
                dx=  ch_geo_data_tw.loc[head_segment]["dx"] # [meter]
                tlf= qlat_tw.loc[head_segment][tsi] # [m^2/s]
                ubcd_g[tsi, frj]= tlf*dx # [m^3/s]
                qlat_g[tsi,0,frj]=0.0                   
    
    return ubcd_g 
                        
#------------------------------------------------------------------------------------------
#       downstream boundary condition mapping between Python and Fortran
# 
#   ** In this mapping, USGS observed stage data at TW is used.
#------------------------------------------------------------------------------------------
def fp_dbcd_map(usgsID2tw
                , usgssDT
                , usgseDT
                , usgspCd
                ):
   
    
    # ** 1) downstream stage (here, lake elevation) boundary condition
    #from nwis_client.iv import IVDataService
    from evaluation_tools.nwis_client.iv import IVDataService
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
    observations_data = IVDataService.get(
                                            sites= usgsID2tw,   #sites='01646500,0208758850', 
                                            startDT=usgssDT,      #'2018-08-01', 
                                            endDT=usgseDT,         #'2020-09-01', 
                                            #parameterCd='62614'
                                            parameterCd=usgspCd)  
    nts_db_g=len(observations_data)-4
    #print(f"nts_db_g:{nts_db_g}")
    #test
    #for tsi in range(0,nts_db_g):
    #    print(f"observatons_data date:{observations_data.iloc[tsi,0]} value:{observations_data.iloc[tsi,4]}")        
    #for tsi in range(0,nts_db_g):
    #    i=tsi+4
    #    print(f"observatons_data shifted date:{observations_data.iloc[i,0]} value:{observations_data.iloc[i,4]}")
    

    # ** 4 is added to make data used here has its date time as from startDT 00:00 to (endDT-1day) 23:00, UTC 
    # ** usgs data at this site uses NGVD1929 feet datum while 'alt' of RouteLink uses NAD88 meter datum.
    #  -> Has to convert accordingly !!!!
    # 
    #source: https://pubs.usgs.gov/sir/2010/5040/section.html
    # Over most USGS study area it is used that NGVD = NAVD88 - 3.6 feet
    dbcd_g= np.zeros(nts_db_g)
    for tsi in range(0,nts_db_g):
        i=tsi+4
        dmy= observations_data.iloc[i,4]   
        dbcd_g[tsi]= 0.3048*(dmy+3.6)  # accuracy with +-0.5feet for 95 percent of USGS study area. 
                                                          # 0.3048 to covert ft to meter. [meter]
    
    return nts_db_g, dbcd_g

def diffusive_input_data_v02(connections
                            , rconn
                            , reaches_bytw
                            , diffusive_parameters
                            , param_df
                            , qlats
                            ):

    usgs_retrievaltool_path= diffusive_parameters.get("usgs_retrievaltool_path",None)
    sys.path.append(usgs_retrievaltool_path)

    # diffusive time steps info.
    dt_ql_g=diffusive_parameters.get("dt_qlat",None) # time step of lateral flow
    dt_ub_g=diffusive_parameters.get("dt_upstream_boundary",None) # time step of us.boundary data
    dt_db_g=diffusive_parameters.get("dt_downstream_boundary",None) # time step of ds.boundary data
    saveinterval_g=diffusive_parameters.get("dt_output",None) # time step for outputting routed results
    saveinterval_ev_g=diffusive_parameters.get("dt_output",None) # time step for evaluating routed results
    dtini_g=diffusive_parameters.get("dt_diffusive",None) # initial simulation time step
    t0_g=0.0 #simulation start hr **set to zero for Fortran computation 
    tfin_g=diffusive_parameters.get("simulation_end_hr",None) # simulation end time

    # USGS data related info.
    usgsID= diffusive_parameters.get("usgsID",None)
    seg2usgsID= diffusive_parameters.get("link2usgsID",None)
    usgssDT= diffusive_parameters.get("usgs_start_date",None)
    usgseDT= diffusive_parameters.get("usgs_end_date",None)
    usgspCd= diffusive_parameters.get("usgs_parameterCd",None)

    # diffusive parameters
    cfl_g= diffusive_parameters.get("courant_number_upper_limit",None)
    theta_g= diffusive_parameters.get("theta_parameter",None)
    tzeq_flag_g= diffusive_parameters.get("chgeo_computation_flag",None)
    y_opt_g= diffusive_parameters.get("water_elevation_computation_flag",None)
    so_llm_g= diffusive_parameters.get("bed_slope_lower_limit",None)

    for twi, (tw, reach_list) in enumerate(reaches_bytw.items(), 1):
    
        # TO DO: once normal depth lower boundary condition capability enabled, remove if statement
        # ??? consider how to use different lower boundary data/methods for different tailwater segs.
        
        if tw==933020089:  #if tw==8777215 or tw==166737669:
            # downstream boundary (tw) segment ID -> make an additional fake tw segment
            dbfksegID= str(tw)+ str(2)
            dbfksegID= int(dbfksegID) 
            ordered_reaches={}
            rchhead_reaches={}
            rchbottom_reaches={}
            z_all={}
    
            flat_list=list(itertools.chain(*reaches_bytw[tw])) # a list of all segments of reaches_bytw for a given tw.
            rconn_tw={key:value for key, value in rconn.items() if key in flat_list} # subset rconn by flat_list
            connections_tw= {key:value for key, value in connections.items() if key in flat_list} # subset connections by flat_list

            path_func = partial(nhd_network.split_at_junction, rconn_tw)
            tr = nhd_network.dfs_decomposition_depth_tuple(rconn_tw, path_func)
            jorder_reaches_tw=sorted(tr, key=lambda x: x[0]) # [ (jorder:[segments]), ... , (jorder:[segments]) ] 

            mx_jorder_tw=max(jorder_reaches_tw)[0] # maximum junction order of subnetwork of TW
            nrch_g=len(jorder_reaches_tw) # the number of reaches        
            maxlist=max(jorder_reaches_tw, key=lambda i:len(i[1]))
            mxncomp_g= len(maxlist[1])+1 # max. number of nodes (segments+one additional segment) within a reach

            for i in jorder_reaches_tw:
                # add one more segment(fake) to the end of a list of segments to account for node configuration.
                fksegID= i[1][len(i[1])-1]
                fksegID= int(str(fksegID) + str(2))
                i[1].append(fksegID)
                # additional segment(fake) to upstream bottom segments
                fk_usbseg=[int(str(x)+str(2)) for x in rconn_tw[i[1][0]]]            

                if i[0] not in ordered_reaches:
                    ordered_reaches.update({i[0]:[]})
                ordered_reaches[i[0]].append([i[1][0],{'number_segments':len(i[1]),\
                                            'segments_list':i[1],\
                                            'upstream_bottom_segments':fk_usbseg,\
                                            'downstream_head_segment':connections_tw[i[1][len(i[1])-2]]}]) 

                if i[1][0] not in rchhead_reaches:    
                    # a list of segments for a given head segment
                    rchhead_reaches.update({i[1][0]:{"number_segments":len(i[1]),\
                                                "segments_list":i[1]}})
                    # a list of segments for a given bottom segment
                    rchbottom_reaches.update({i[1][len(i[1])-1]:{"number_segments":len(i[1]),\
                                                         "segments_list":i[1]}})
                # for channel altitude adjustment
                z_all.update({seg:{'adj.alt':np.zeros(1)}
                                            for seg in i[1]})
            # cahnnel geometry data
            ch_geo_data_tw = param_df.loc[
            flat_list, ["bw", "tw", "twcc", "dx", "n", "ncc", "cs", "s0", "alt"]]
            ch_geo_data_tw[:]["cs"]= 1.0/ch_geo_data_tw[:]["cs"]
        #--------------------------------------------------------------------------------------
        #                                 Step 0-3           

        #    Adjust altitude so that altitude of the last sement of a reach is equal to that 
        #    of the first segment of its downstream reach right after their common junction.
        #--------------------------------------------------------------------------------------
            adj_alt1(mx_jorder_tw
                        , ordered_reaches
                        , ch_geo_data_tw
                        , dbfksegID
                        , z_all
                        ) 
        #--------------------------------------------------------------------------------------
        #                                 Step 0-4           

        #     Make Fortran-Python channel network mapping variables.
        #--------------------------------------------------------------------------------------   
            pynw={}
            frj=-1
            for x in range(mx_jorder_tw,-1,-1): 
                for head_segment, reach in ordered_reaches[x]:
                    frj= frj+1
                    pynw[frj]=head_segment

            #frnw_col=8
            frnw_col= diffusive_parameters.get("fortran_nework_map_col_number",None)
            frnw_g=fp_network_map(mx_jorder_tw
                    , ordered_reaches
                    , rchbottom_reaches
                    , nrch_g
                    , frnw_col
                    , dbfksegID
                    , pynw
                    )  
            #covert data type from integer to float for frnw
            dfrnw_g=np.zeros((nrch_g,frnw_col), dtype=float)
            for j in range(0,nrch_g):
                for col in range(0,frnw_col):
                    dfrnw_g[j,col]=float(frnw_g[j,col])
        #---------------------------------------------------------------------------------
        #                              Step 0-5

        #                  Prepare channel geometry data           
        #---------------------------------------------------------------------------------    
            z_ar_g, bo_ar_g, traps_ar_g, tw_ar_g, twcc_ar_g, mann_ar_g, manncc_ar_g, so_ar_g, dx_ar_g= fp_chgeo_map(mx_jorder_tw
                        , ordered_reaches
                        , ch_geo_data_tw
                        , z_all
                        , mxncomp_g
                        , nrch_g                    
                        )   
        #---------------------------------------------------------------------------------
        #                              Step 0-6

        #                  Prepare lateral inflow data           
        #---------------------------------------------------------------------------------
            segs = list(chain.from_iterable(reach_list))
            qlat_tw = qlats.loc[segs]    
            #tfin_g=len(qlat_tw.columns)-1 #entire simulation period in hrs
            nts_ql_g= int((tfin_g-t0_g)*3600.0/dt_ql_g)+1 # the number of the entire time steps of lateral flow data 

            qlat_g=np.zeros((nts_ql_g, mxncomp_g, nrch_g)) 

            fp_qlat_map(mx_jorder_tw
                , ordered_reaches
                , nts_ql_g
                , qlat_tw            
                , qlat_g
                ) 
        #---------------------------------------------------------------------------------
        #                              Step 0-7

        #       Prepare upstream boundary (top segments of head basin reaches) data            
        #---------------------------------------------------------------------------------
            nts_ub_g= nts_ql_g 
            ubcd_g = fp_ubcd_map(frnw_g
                                    , pynw
                                    , nts_ub_g
                                    , nrch_g
                                    , ch_geo_data_tw
                                    , qlat_tw
                                    , qlat_g
                                    )
        #---------------------------------------------------------------------------------
        #                              Step 0-8

        #       Prepare downstrea boundary (bottom segments of TW reaches) data            
        #---------------------------------------------------------------------------------        
            #import pdb; pdb.set_trace()
            if tw in seg2usgsID:
                ipos= seg2usgsID.index(tw)
                usgsID2tw= usgsID[ipos]         
                nts_db_g, dbcd_g=fp_dbcd_map(usgsID2tw
                            , usgssDT
                            , usgseDT
                            , usgspCd
                            )
            else:
                # no usgs data available at this TW.
                nts_db_g=-1.0
                
        #---------------------------------------------------------------------------------
        #                              Step 0-8

        #                 Prepare uniform flow lookup tables            
        #---------------------------------------------------------------------------------          
            #nhincr_m_g=20
            #nhincr_f_g=20               
            #timesdepth_g=10
            nhincr_m_g= diffusive_parameters.get("normaldepth_lookuptable_main_increment_number",None) 
            nhincr_f_g= diffusive_parameters.get("normaldepth_lookuptable_floodplain_increment_number",None) 
            timesdepth_g= diffusive_parameters.get("normaldepth_lookuptable_depth_multiplier",None)
            ufqlt_m_g= np.zeros((mxncomp_g,nrch_g,nhincr_m_g))                        
            ufhlt_m_g= np.zeros((mxncomp_g,nrch_g,nhincr_m_g))
            ufqlt_f_g= np.zeros((mxncomp_g,nrch_g,nhincr_f_g))
            ufhlt_f_g= np.zeros((mxncomp_g,nrch_g,nhincr_f_g))
            #ufhlt_m_g, ufqlt_m_g, ufhlt_f_g, ufqlt_f_g= uniflow_lookuptable(mxncomp_g 
            #                                                            , nrch_g 
            #                                                            , bo_ar_g 
            #                                                            , traps_ar_g 
            #                                                            , tw_ar_g 
            #                                                            , twcc_ar_g 
            #                                                            , mann_ar_g 
            #                                                            , manncc_ar_g 
            #                                                            , so_ar_g 
            #                                                            , nhincr_m_g 
            #                                                            , nhincr_f_g               
            #                                                            , frnw_col 
            #                                                            , dfrnw_g 
            #                                                            , timesdepth_g)

        #---------------------------------------------------------------------------------
        #                              Step 0-9

        #                       Build input dictionary           
        #---------------------------------------------------------------------------------  
            ntss_ev_g= int((tfin_g - t0_g)*3600.0/saveinterval_ev_g)+1 

            
            # build a dictionary of diffusive model inputs
            diff_ins = {}
            diff_ins["dtini_g"] = dtini_g
            diff_ins["t0_g"] = t0_g
            diff_ins["tfin_g"] = tfin_g
            diff_ins["saveinterval_g"] = saveinterval_g
            diff_ins["saveinterval_ev_g"] = saveinterval_ev_g
            diff_ins["dt_ql_g"] = dt_ql_g
            diff_ins["dt_ub_g"] = dt_ub_g
            diff_ins["dt_db_g"] = dt_db_g
            diff_ins["nts_ql_g"] = nts_ql_g
            diff_ins["nts_ub_g"] = nts_ub_g
            diff_ins["nts_db_g"] = nts_db_g
            diff_ins["mxncomp_g"] = mxncomp_g
            diff_ins["nrch_g"] = nrch_g
            diff_ins["z_ar_g"] = z_ar_g
            diff_ins["bo_ar_g"] = bo_ar_g
            diff_ins["traps_ar_g"] = traps_ar_g
            diff_ins["tw_ar_g"] = tw_ar_g
            diff_ins["twcc_ar_g"] = twcc_ar_g
            diff_ins["mann_ar_g"] = mann_ar_g
            diff_ins["manncc_ar_g"] = manncc_ar_g
            diff_ins["so_ar_g"] = so_ar_g
            diff_ins["dx_ar_g"] = dx_ar_g
            diff_ins["nhincr_m_g"] = nhincr_m_g
            diff_ins["nhincr_f_g"] = nhincr_f_g
            diff_ins["ufhlt_m_g"] = ufhlt_m_g
            diff_ins["ufqlt_m_g"] = ufqlt_m_g
            diff_ins["ufhlt_f_g"] = ufhlt_f_g
            diff_ins["ufqlt_f_g"] = ufqlt_f_g
            diff_ins["frnw_col"] = frnw_col
            diff_ins["frnw_g"] = frnw_g
            diff_ins["qlat_g"] = qlat_g
            diff_ins["ubcd_g"] = ubcd_g
            diff_ins["dbcd_g"] = dbcd_g
            diff_ins["cfl_g"] = cfl_g
            diff_ins["theta_g"] = theta_g
            diff_ins["tzeq_flag_g"] = tzeq_flag_g
            diff_ins["y_opt_g"] = y_opt_g
            diff_ins["so_llm_g"] = so_llm_g
            diff_ins["ntss_ev_g"] = ntss_ev_g
            
#             # save input data as yaml
#             import yaml
#             with open('diff_inputs.yml', 'w') as outfile:
#                 yaml.dump(diff_ins, outfile, default_flow_style=False)

    return diff_ins