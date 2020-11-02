import numpy as np
import dfcnhi as df

#--------------------------------------------------------------------------------------------------
#   The genearl purpose of this py file is to create lookup tables for finding water depth or 
#  elevation from given discharge values
#--------------------------------------------------------------------------------------------------
    
    

# ** Create uniform flow-normal depth lookup tables for trapezoidal main + rectangular
#  floodplain using channel geometry equations not channel geometry lookup tables
def ufTR_ltable(        
        connections= None 
        ,supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , nhincr_m= None
        , nhincr_f= None
        , uflt= None
        ):   

    for x in range(network['maximum_reach_seqorder'],-1,-1):  
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= seg_list_all[head_segment]
            ncomp=ncompall[head_segment]                   

            df.var.ncomp=ncomp                   
            #df.var.z_ar=np.zeros(ncomp)
            df.var.bo_ar=np.zeros(ncomp)
            df.var.traps_ar=np.zeros(ncomp)
            df.var.tw_ar=np.zeros(ncomp)
            df.var.twcc_ar=np.zeros(ncomp)
            df.var.mann_ar=np.zeros(ncomp)
            df.var.manncc_ar=np.zeros(ncomp)
            df.var.so_ar=np.zeros(ncomp)                        

            for seg in range(0,ncomp):                    
                if seg==ncomp-1:
                    segID= seg_list[seg-1]
                else:
                    segID= seg_list[seg]                      
                df.var.bo_ar[seg]=connections[segID]['data'][supernetwork_data['bottomwidth_col']]

                traps=connections[segID]['data'][supernetwork_data['ChSlp_col']]
                df.var.traps_ar[seg]=1.0/traps  # s in 1 over s for channel side slope    

                df.var.tw_ar[seg]=connections[segID]['data'][supernetwork_data['topwidth_col']]
                df.var.twcc_ar[seg]=connections[segID]['data'][supernetwork_data['topwidthcc_col']]
                df.var.mann_ar[seg]=connections[segID]['data'][supernetwork_data['manningn_col']]
                df.var.manncc_ar[seg]=connections[segID]['data'][supernetwork_data['manningncc_col']]
                df.var.so_ar[seg]=connections[segID]['data'][supernetwork_data['slope_col']]
                #df.var.dx[seg]=connections[segID]['data'][supernetwork_data['length_col']]
            df.var.ufqlt_m=np.zeros((nhincr_m,ncomp))                        
            df.var.ufhlt_m=np.zeros((nhincr_m,ncomp))
            df.var.ufqlt_f=np.zeros((nhincr_f,ncomp))
            df.var.ufhlt_f=np.zeros((nhincr_f,ncomp))
            df.var.nhincr_m= nhincr_m
            df.var.nhincr_f= nhincr_f

            df.flowlt.uniflowlt_tz()


            for seg in range(0,ncomp):                    
                segID= seg_list[seg]
                for i1 in range(0,nhincr_m):
                    uflt[segID]['q_m'][i1]= df.var.ufqlt_m[i1,seg]
                    uflt[segID]['h_m'][i1]= df.var.ufhlt_m[i1,seg]
                for i1 in range(0,nhincr_f):                        
                    uflt[segID]['q_f'][i1]= df.var.ufqlt_f[i1,seg]
                    uflt[segID]['h_f'][i1]= df.var.ufhlt_f[i1,seg]
    
    return uflt