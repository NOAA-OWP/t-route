
#--------------------------------------------------------------------------------------------------
#   The genearl purpose of this py file is to adjust channel geometries stored in Route_Link.nc 
#   for any hydraulic routing needs.
#--------------------------------------------------------------------------------------------------

# Adjust altitude so that altitude of the last node (=fake segment) of a reach is equal to that of
# the first segment of its downstream reach right after their common junction.
def adj_alt1(mx_jorder_tw
        , ordered_reaches
        , ch_geo_data_tw 
        , dbfksegID
        , z_all   
        ):   

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
               