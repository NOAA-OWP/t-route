import numpy as np
import dfcnhi as df

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
        , nxsecpt= None
        , nel= None
        , timesdepth= None
        , tzeq_flag= None
        , xsec_attr_all= None
        ):   
    # x-sec attribute lookup table variables    
    #nxsecpt=8 #the number of (x,y) points for x-sec geometry
    df.var.nxsecpt=nxsecpt
    df.var.xcs=np.zeros(nxsecpt)
    df.var.ycs=np.zeros(nxsecpt)
    #nel=101 # the number of incrementally increasing water depth lines
    df.var.nel= nel
    #timesdepth=5.0  # multiplier to bankfull depth to cover potentially largest water depth
    df.var.timesdepth=timesdepth  
    
    xsec_attr_seg=np.zeros((4,nel))
    df.var.xsec_attr_seg=xsec_attr_seg
    arrxsx=np.zeros(nel)
    arrxsy=np.zeros(nel)      

    for x in range(network['maximum_reach_seqorder'],-1,-1): 
        for head_segment, reach in ordered_reaches[x]:
            seg_list= seg_list_all[head_segment]  
            ncomp=ncompall[head_segment] 

            for seg in range(0,ncomp):
                if seg==ncomp-1:
                # when seg is for fake segment after the last segment of a reach,
                # use the same channel attributes for the last segment first and
                # stores ensuing results to fake segmentID of related var.
                    segID= seg_list[seg-1]
                else:
                    segID= seg_list[seg]     
                so=connections[segID]['data'][supernetwork_data['slope_col']]
                bo=connections[segID]['data'][supernetwork_data['bottomwidth_col']]
                traps=connections[segID]['data'][supernetwork_data['ChSlp_col']]
                traps=1.0/traps # s in 1 over s for channel side slope
                tw=connections[segID]['data'][supernetwork_data['topwidth_col']]
                twcc=connections[segID]['data'][supernetwork_data['topwidthcc_col']]
                mann= connections[segID]['data'][supernetwork_data['manningn_col']]
                manncc= connections[segID]['data'][supernetwork_data['manningncc_col']]

                segID1= seg_list[seg]  
                z= z_all[segID1]['adj.alt'][0]
              
                #compute (x,y) points for x-sec, p131,rm4
                hbf= (tw-bo)/(2.0*traps)  #bankfull depth
                df.var.xcs[0]=0.0
                df.var.xcs[1]=df.var.xcs[0]
                df.var.xcs[2]=(twcc-tw)/2.0
                df.var.xcs[3]=df.var.xcs[2]+traps*hbf
                df.var.xcs[4]=df.var.xcs[3]+bo
                df.var.xcs[5]=df.var.xcs[4]+traps*hbf
                df.var.xcs[6]=twcc
                df.var.xcs[7]=df.var.xcs[6]

                df.var.ycs[0]=z+timesdepth*hbf
                df.var.ycs[1]=z+hbf
                df.var.ycs[2]=z+hbf
                df.var.ycs[3]=z
                df.var.ycs[4]=z
                df.var.ycs[5]=z+hbf
                df.var.ycs[6]=z+hbf
                df.var.ycs[7]=z+timesdepth*hbf
                #for Fortran subroutine readXsection
                df.var.mann= mann
                df.var.manncc=manncc                
                # As an alternative, analytically compute channel properties
                df.var.tzeq_flag= tzeq_flag #0 deactivates / 1 activates it
                df.var.z0=z 
                df.var.so0=so
                df.var.bo0=bo
                df.var.traps0=traps 
                df.var.tw0=tw 
                df.var.twcc0=twcc
                # end of the inputs for the alternative.

                df.attrtable.read_xsec()

                #if seg==ncomp-1:
                # when seg is for fake segment after the last segment of a reach,
                # stores ensuing results to fake segmentID to related var.   
                #    segID= fksegID_dsend[head_segment]
                #import pdb; pdb.set_trace() 
                segID2= seg_list[seg] 

                xsec_attr_all[segID2]['elev']=df.var.xsec_attr_seg[0]
                xsec_attr_all[segID2]['convey']=df.var.xsec_attr_seg[1]
                xsec_attr_all[segID2]['topwd']=df.var.xsec_attr_seg[2]
                xsec_attr_all[segID2]['uniflow']=df.var.xsec_attr_seg[3]
    
    return xsec_attr_all
               