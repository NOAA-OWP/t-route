import numpy as np
#import dfcnhi as df
import dfcnhi_mt as df  

#--------------------------------------------------------------------------------------------------
#   The genearl purpose of this py file is to create lookup tables for finding water depth or 
#  elevation from given discharge values
#--------------------------------------------------------------------------------------------------
# ** Create uniform flow-normal depth lookup tables for trapezoidal main + rectangular
#  floodplain using channel geometry equations not channel geometry lookup tables
# ** Fortran handles all the reaches with given Fortran-Python mapping code,frnw.

def ufTR_ltable_frnw(        
        connections= None 
        ,supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None
        , mxncomp_g= None    
        , nrch_g= None
        , frnw_g= None
        , timesdepth_g= None
        , nhincr_m_g= None
        , nhincr_f_g= None
        , ufqlt_m_g= None
        , ufhlt_m_g= None
        , ufqlt_f_g= None
        , ufhlt_f_g= None 
        , bo_ar_g= None
        , traps_ar_g= None
        , tw_ar_g= None
        , twcc_ar_g= None
        , mann_ar_g= None
        , manncc_ar_g= None
        , so_ar_g= None
        ):   
    df.var.mxncomp_g= mxncomp_g
    df.var.nrch_g= nrch_g
    
    df.var.frnw_g= np.zeros((nrch_g,8), dtype=int)
    df.var.frnw_g= frnw_g
    
    df.var.timesdepth_g= timesdepth_g   
    
    df.var.nhincr_m_g= nhincr_m_g
    df.var.nhincr_f_g= nhincr_f_g
    df.var.ufqlt_m_g= np.zeros((mxncomp_g,nrch_g,nhincr_m_g))                        
    df.var.ufhlt_m_g= np.zeros((mxncomp_g,nrch_g,nhincr_m_g))
    df.var.ufqlt_f_g= np.zeros((mxncomp_g,nrch_g,nhincr_f_g))
    df.var.ufhlt_f_g= np.zeros((mxncomp_g,nrch_g,nhincr_f_g))

    df.var.bo_ar_g = np.zeros((mxncomp_g,nrch_g))
    df.var.traps_ar_g = np.zeros((mxncomp_g,nrch_g))
    df.var.tw_ar_g = np.zeros((mxncomp_g,nrch_g))
    df.var.twcc_ar_g = np.zeros((mxncomp_g,nrch_g))
    df.var.mann_ar_g = np.zeros((mxncomp_g,nrch_g))
    df.var.manncc_ar_g = np.zeros((mxncomp_g,nrch_g))
    df.var.so_ar_g = np.zeros((mxncomp_g,nrch_g))   
    
    df.var.bo_ar_g = bo_ar_g
    df.var.traps_ar_g = traps_ar_g
    df.var.tw_ar_g = tw_ar_g
    df.var.twcc_ar_g = twcc_ar_g
    df.var.mann_ar_g = mann_ar_g
    df.var.manncc_ar_g = manncc_ar_g
    df.var.so_ar_g = so_ar_g     
    
    df.flowlt.uniflowlt_tz_alrch()

    for frj in range(0,nrch_g):
        for i in range(0,mxncomp_g):
            for i1 in range(0,nhincr_m_g):
                ufqlt_m_g[i,frj,i1] = df.var.ufqlt_m_g[i,frj,i1]                    
                ufhlt_m_g[i,frj,i1] = df.var.ufhlt_m_g[i,frj,i1]
            for i1 in range(0,nhincr_f_g):
                ufqlt_f_g[i,frj,i1] = df.var.ufqlt_f_g[i,frj,i1] 
                ufhlt_f_g[i,frj,i1] = df.var.ufhlt_f_g[i,frj,i1]   
    
    #import pdb; pdb.set_trace()
#    for seg in range(0,ncomp):                    
#        segID= seg_list[seg]
#        for i1 in range(0,nhincr_m):
#            uflt[segID]['q_m'][i1]= df.var.ufqlt_mr[i1,seg]
#            uflt[segID]['h_m'][i1]= df.var.ufhlt_mr[i1,seg]
#        for i1 in range(0,nhincr_f):                        
#            uflt[segID]['q_f'][i1]= df.var.ufqlt_fr[i1,seg]
#            uflt[segID]['h_f'][i1]= df.var.ufhlt_fr[i1,seg]

    return ufqlt_m_g
    return ufhlt_m_g
    return ufqlt_f_g
    return ufhlt_f_g


    

# ** Create uniform flow-normal depth lookup tables for trapezoidal main + rectangular
#  floodplain using channel geometry equations not channel geometry lookup tables
def ufTR_ltable(        
        connections= None 
        ,supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , timesdepth= None
        , nhincr_m= None
        , nhincr_f= None
        , uflt= None
        ):   
    df.var.timesdepth= timesdepth
    
    for x in range(network['maximum_reach_seqorder'],-1,-1):  
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= seg_list_all[head_segment]
            ncomp=ncompall[head_segment]                   

            df.var.ncomp_unif=ncomp                   
            #df.var.z_ar=np.zeros(ncomp)
            df.var.bo_r=np.zeros(ncomp)
            df.var.traps_r=np.zeros(ncomp)
            df.var.tw_r=np.zeros(ncomp)
            df.var.twcc_r=np.zeros(ncomp)
            df.var.mann_r=np.zeros(ncomp)
            df.var.manncc_r=np.zeros(ncomp)
            df.var.so_r=np.zeros(ncomp)                        

            for seg in range(0,ncomp):                    
                if seg==ncomp-1:
                    segID= seg_list[seg-1]
                else:
                    segID= seg_list[seg]                      
                df.var.bo_r[seg]=connections[segID]['data'][supernetwork_data['bottomwidth_col']]

                traps=connections[segID]['data'][supernetwork_data['ChSlp_col']]
                df.var.traps_r[seg]=1.0/traps  # s in 1 over s for channel side slope    

                df.var.tw_r[seg]=connections[segID]['data'][supernetwork_data['topwidth_col']]
                df.var.twcc_r[seg]=connections[segID]['data'][supernetwork_data['topwidthcc_col']]
                df.var.mann_r[seg]=connections[segID]['data'][supernetwork_data['manningn_col']]
                df.var.manncc_r[seg]=connections[segID]['data'][supernetwork_data['manningncc_col']]
                df.var.so_r[seg]=connections[segID]['data'][supernetwork_data['slope_col']]
                #df.var.dx[seg]=connections[segID]['data'][supernetwork_data['length_col']]
            df.var.ufqlt_mr=np.zeros((nhincr_m,ncomp))                        
            df.var.ufhlt_mr=np.zeros((nhincr_m,ncomp))
            df.var.ufqlt_fr=np.zeros((nhincr_f,ncomp))
            df.var.ufhlt_fr=np.zeros((nhincr_f,ncomp))
            df.var.nhincr_m= nhincr_m
            df.var.nhincr_f= nhincr_f

            df.flowlt.uniflowlt_tz()


            for seg in range(0,ncomp):                    
                segID= seg_list[seg]
                for i1 in range(0,nhincr_m):
                    uflt[segID]['q_m'][i1]= df.var.ufqlt_mr[i1,seg]
                    uflt[segID]['h_m'][i1]= df.var.ufhlt_mr[i1,seg]
                for i1 in range(0,nhincr_f):                        
                    uflt[segID]['q_f'][i1]= df.var.ufqlt_fr[i1,seg]
                    uflt[segID]['h_f'][i1]= df.var.ufhlt_fr[i1,seg]
    
    return uflt



