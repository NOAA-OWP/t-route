import os
#--------------------------------------------------------------------------------
# Case1: 
#      - set upstream boundary data only at headwater segments using lateral flow  
#      - set downstream boundary data based on USGS observed data 
#---------------------------------------------------------------------------------
def set_bdrydata_usgs1(        
        connections= None  
        , supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , fksegID_dsend= None
        , dbfksegID= None
        , dbseg_usgsID= None
        , nts_ql_g= None
        , nts_ub_g= None
        , nts_db_g= None     
        , ubcond= None
        , dbcond= None
        , usgssDT= None
        , usgseDT= None         
        , qlatral= None
        , output_path= None       
        ):   

    
    # ** 1) downstream stage (here, lake elevation) boundary condition
    #from nwis_client.iv import IVDataService
    from iv import IVDataService
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
                                            sites=dbseg_usgsID,   #sites='01646500,0208758850', 
                                            startDT=usgssDT,      #'2018-08-01', 
                                            endDT=usgseDT,         #'2020-09-01', 
                                            parameterCd='62614')  
 
    # 4 is added to make data used here has its date time as from startDT 00:00 to (endDT-1day) 23:00, UTC 
# *** usgs data at this site uses NGVD1929 feet datum while 'alt' of RouteLink uses NAD88 meter datum.
#  -> Has to convert accordingly !!!!
    # 
    #source: https://pubs.usgs.gov/sir/2010/5040/section.html
    # Over most USGS study area it is used that NGVD = NAVD88 - 3.6 feet
    for tsi in range(0,nts_db_g):
        i=tsi+4
        dmy= observations_data.iloc[i,4]   
        dbcond[dbfksegID]['elev'][tsi]= 0.3048*(dmy+3.6)  # accuracy with +-0.5feet for 95 percent of USGS study area. 
                                                          # 0.3048 to covert ft to meter. [meter]
    
    #import pdb; pdb.set_trace()
   
    # upstream boundary condition
    for x in range(network['maximum_reach_seqorder'],-1,-1):   
        for head_segment, reach in ordered_reaches[x]:      
            if reach['upstream_reaches']=={0}:           
            # upstream boundary condition for headbasin segment is assumed to be zero discharge [m^3/sec]   
                for tsi in range(0,nts_ub_g): 
                    ubcond[head_segment]['discharge'][tsi]= 0.0
                    #ubcond[head_segment]['discharge'][tsi]=qlatral[head_segment]['qlat'][tsi]
                    
    with open(os.path.join(output_path,"dbcond"),'a') as op_db:
        for tsi in range(0,nts_db_g):            
            op_db.write("%s %s %s\n" % (tsi,dbfksegID,dbcond[dbfksegID]['elev'][tsi])) 
 
    with open(os.path.join(output_path,"ubcond"),'a') as op_ub:
        for x in range(network['maximum_reach_seqorder'],-1,-1): 
            for head_segment, reach in ordered_reaches[x]:      
                if reach['upstream_reaches']=={0}:              
                    for tsi in range(0,nts_ub_g):  
                        op_ub.write("%s %s %s\n" % (tsi,head_segment,ubcond[head_segment]['discharge'][tsi])) 
                      
    
    return ubcond
    return dbcond
       
