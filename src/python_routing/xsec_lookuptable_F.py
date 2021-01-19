import numpy as np
#import dfcnhi as df
import dfcnhi_mt as df  

#---------------------------------------------------------------------------------------------------

#          Make x-section attributes lookup table for trapezoidal main + rectangular floodplain 
# 
# (p.132,rm4)
#---------------------------------------------------------------------------------------------------
def xsecTR_ltable(        
        connections= None 
        , supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , z_all= None
        , nxsecpt_g= None
        , nel_g= None
        , timesdepth_g= None
        , tzeq_flag_g= None
        , xsec_attr_all= None
        ):   
    
    
    # x-sec attribute lookup table variables    
    #nxsecpt_g=8 #the number of (x,y) points for x-sec geometry
    df.var.nxsecpt_g=nxsecpt_g
    df.var.xcs_g=np.zeros(nxsecpt_g)
    df.var.ycs_g=np.zeros(nxsecpt_g)
    #nel_g=101 # the number of incrementally increasing water depth lines
    df.var.nel_g= nel_g
    #timesdepth: multiplier to bankfull depth to cover potentially largest water depth
    df.var.timesdepth_g=timesdepth_g  
    
    xsec_attr_seg_g=np.zeros((4,nel_g))
    df.var.xsec_attr_seg_g=xsec_attr_seg_g
    arrxsx=np.zeros(nel_g)
    arrxsy=np.zeros(nel_g)      

    for x in range(network['maximum_reach_seqorder'],-1,-1): 
        for head_segment, reach in ordered_reaches[x]:
            seg_list= seg_list_all[head_segment]  
            ncomp=ncompall[head_segment] 

            for seg in range(0,ncomp):
                if seg==ncomp-1:
                # when seg is for fake segment after the last segment of a reach,
                # use the same channel_g attributes for the last segment first and
                # stores ensuing results to fake segmentID of related var.
                    segID= seg_list[seg-1]
                else:
                    segID= seg_list[seg]     
                so=connections[segID]['data'][supernetwork_data['slope_col']]
                bo=connections[segID]['data'][supernetwork_data['bottomwidth_col']]
                traps=connections[segID]['data'][supernetwork_data['ChSlp_col']]
                traps=1.0/traps # s in 1 over s for channel_g side slope
                tw=connections[segID]['data'][supernetwork_data['topwidth_col']]
                twcc=connections[segID]['data'][supernetwork_data['topwidthcc_col']]
                mann= connections[segID]['data'][supernetwork_data['manningn_col']]
                manncc= connections[segID]['data'][supernetwork_data['manningncc_col']]

                segID1= seg_list[seg]  
                z= z_all[segID1]['adj.alt'][0]
              
                #compute (x,y) points for x-sec, p131,rm4
                hbf= (tw-bo)/(2.0*traps)  #bankfull depth
                df.var.xcs_g[0]=0.0
                df.var.xcs_g[1]=df.var.xcs_g[0]
                df.var.xcs_g[2]=(twcc-tw)/2.0
                df.var.xcs_g[3]=df.var.xcs_g[2]+traps*hbf
                df.var.xcs_g[4]=df.var.xcs_g[3]+bo
                df.var.xcs_g[5]=df.var.xcs_g[4]+traps*hbf
                df.var.xcs_g[6]=twcc
                df.var.xcs_g[7]=df.var.xcs_g[6]

                df.var.ycs_g[0]=z+timesdepth_g*hbf
                df.var.ycs_g[1]=z+hbf
                df.var.ycs_g[2]=z+hbf
                df.var.ycs_g[3]=z
                df.var.ycs_g[4]=z
                df.var.ycs_g[5]=z+hbf
                df.var.ycs_g[6]=z+hbf
                df.var.ycs_g[7]=z+timesdepth_g*hbf
                #for Fortran subroutine readXsection
                df.var.mann= mann
                df.var.manncc=manncc                
                # As an alternative, analytically compute channel_g properties
                df.var.tzeq_flag_g= tzeq_flag_g #0 deactivates / 1 activates it
                df.var.z_g=z 
                df.var.so_g=so
                df.var.bo_g=bo
                df.var.traps_g=traps 
                df.var.tw_g=tw 
                df.var.twcc_g=twcc
                df.var.mann_g= mann
                df.var.manncc_g= manncc
                # end of the inputs for the alternative.
                
                #import pdb; pdb.set_trace()
                df.attrtable.read_xsec()

                #if seg==ncomp-1:
                # when seg is for fake segment after the last segment of a reach,
                # stores ensuing results to fake segmentID to related var.   
                #    segID= fksegID_dsend[head_segment]
                #import pdb; pdb.set_trace() 
                segID2= seg_list[seg] 
                
                for il in range(nel_g):
                    xsec_attr_all[segID2]['elev'][il]= df.var.xsec_attr_seg_g[0,il]
                    xsec_attr_all[segID2]['convey'][il]= df.var.xsec_attr_seg_g[1,il]
                    xsec_attr_all[segID2]['topwd'][il]= df.var.xsec_attr_seg_g[2,il]
                    xsec_attr_all[segID2]['uniflow'][il]= df.var.xsec_attr_seg_g[3,il]
                
                            #xsec_attr_all.update({fksegID:{'elev':np.zeros(nel_g)
                            #            ,'convey':np.zeros(nel_g)
                            #            ,'topwd':np.zeros(nel_g)
                            #            ,'uniflow':np.zeros(nel_g)}})
    
    return xsec_attr_all
               