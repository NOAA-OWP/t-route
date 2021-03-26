import os
import numpy as np

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


## code storage:

   #output_path="./output"
    #with open(os.path.join(output_path,"Py_ubcd"),'a') as pyubcd:
        #for x in range(mx_jorder_tw,-1,-1):
        #    for head_segment, reach in ordered_reaches[x]:                  
        #        frj=frj+1
        #        if not reach['upstream_bottom_segments']:
        #            for tsi in range(0,nts_ub_g):          
        #                dx=  ch_geo_data_tw.loc[head_segment]["dx"] # [meter]
        #                tlf= qlat_tw.loc[head_segment][tsi] # [m^2/s]
        #                ubcd_g[tsi, frj]= tlf*dx # [m^3/s]
                        # change the existing amount of lateral flow at the head segments of headbasin reaches to zero.
        #                qlat_g[tsi,0,frj]=0.0
                        # test
        #                pyubcd.write("%s %s %s %s %s %s %s\n" %\
        #                        (tsi, frj, head_segment, tlf, dx, ubcd_g[tsi, frj], qlat_g[tsi,0,frj]))        
        #pyubcd.write("\n")
                    # test
                    #pyubcd.write("%s %s %s %s %s %s %s\n" %\
                    #            (tsi, frj, head_segment, tlf, dx, ubcd_g[tsi, frj], qlat_g[tsi,0,frj])) 
                    
                    
                    
    #output_path='./output'
    #with open(os.path.join(output_path,"frnw"),'a') as frnw:
    #    for frj in range(0,nrch_g):
    #        frnw.write("%s %s %s %s %s %s %s %s %s %s\n" % (frj, pynw[frj],frnw_g[frj,0],frnw_g[frj,1],frnw_g[frj,2],\
    #                      frnw_g[frj,3],frnw_g[frj,4],frnw_g[frj,5],frnw_g[frj,6],frnw_g[frj,7]))