import numpy as np
import dfcnhi as df

#--------------------------------------------------------------------------------------------------
#   The genearl purpose of this py file is to adjust channel geometries computed in Route_Link.nc 
#   according to any hydraulic routing needs.
#--------------------------------------------------------------------------------------------------

# Adjust altitude so that altitude of the last sement of a reach is equal to that of
# the first segment of its downstream reach right after their common junction.
def adj_alt1(        
        connections= None 
        ,supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , dbfksegID= None
        , z_all= None   
        ):   
    
    for x in range(network['maximum_reach_seqorder'],-1,-1):  
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= seg_list_all[head_segment]                 
            ncomp=ncompall[head_segment] 

            for seg in range(0,ncomp):
                segID= seg_list[seg]
                if seg==ncomp-1 and seg_list.count(dbfksegID)==0 :
                # At junction, the altitude of fake segment of an upstream reach
                # is equal to that of the first segment of the downstream reach
                    # the first segment id of downstream reach from a junction
                    dsrchID= network['reaches'][head_segment]['downstream_reach']                        
                    z_all[segID]['adj.alt'][0]=\
                      connections[dsrchID]['data'][supernetwork_data['alt_col']]
                elif seg==ncomp-1 and seg_list.count(dbfksegID)>0:
                # Terminal downstream fakesegment                
         ## AD HOC: need to be corrected later
                    segID2= seg_list[seg-1]
                    So= connections[segID2]['data'][supernetwork_data['slope_col']]
                    dx= connections[segID2]['data'][supernetwork_data['length_col']]        
                    z_all[segID]['adj.alt'][0]= z_all[segID2]['adj.alt'][0] - So*dx
                else:
                    z_all[segID]['adj.alt'][0]=\
                        connections[segID]['data'][supernetwork_data['alt_col']] 
                #test
                #orgz=-100.0
                #if seg<ncomp-1:
                #    orgz= connections[segID]['data'][supernetwork_data['alt_col']]                
                #print(f"segID:{segID} org.alt:{orgz} adj.alt:{z_all[segID]['adj.alt'][0]}")     
    
    return z_all
               