import numpy as np
import dfcnhi as df

#--------------------------------------------------------------------------------------------------
# Create boundary data such as discharge, water depth, lateral flow by either data import or 
# artificial equations.

#--------------------------------------------------------------------------------------------------
def set_bdrydata(        
        connections= None  
        , supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , fksegID_dsend= None
        , dbfksegID= None
        , ntsi= None
        , bdrydata= None         
        , q_bdseg= None
        , h_bdseg= None
        , qlatral= None 
        ):   
        
    #--------------------------------------------------------------------------------
    #               Artificial input data 1
    
    #--------------------------------------------------------------------------------       
    if bdrydata=="artf_data1":  
    # ** 1) headbasin segment's discharge and downstream terminal segment's water depth  
        
        # for headbasin segment q
        vr=3.0
        mtp=20.0
        mn= 10 #40.0
        yintc=1.0
        # for depth at downstream terminal segment
        mtp2=5.0
        mn2= mn+10.0
        vr2=2.0*vr            

        for x in range(network['maximum_reach_seqorder'],-1,-1):   
            for head_segment, reach in ordered_reaches[x]:                
                yintc= yintc + 0.2
                if reach['upstream_reaches']=={0}:              
                    for tsi in range(0,ntsi):
                        n=tsi+1
                        q_bdseg[head_segment]['known discharge'][tsi]=\
                              df.subtools.input_crtor1(mtp, mn, vr, yintc,n)                                              
                        #q_bdseg[head_segment]['known discharge'][tsi]=\
                        #     mtp*np.exp(-0.5*((float(tsi)-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5)+ yintc 
                   
                
                if x==0:    
                    yintc2=1.0
                    for tsi in range(0,ntsi):
                        n=tsi+1
                        h_bdseg[dbfksegID]['known depth'][tsi]=df.subtools.input_crtor1(mtp2, mn2, vr2, yintc2,n)                                        
                        #h_bdseg[dbfksegID]['known depth'][tsi]=\
                        #    mtp2*np.exp(-0.5*((float(tsi)-mn2)/vr2)**2.0)/(vr2*(2.0*3.14)**0.5)+ yintc2                 
                          
    # ** 2) artificial lateral flow
    
    # ** the unit of qlat that is eventually used in the routing SHOULD BE [m^2/sec].
    # ** qlat of fake segment of the last segment of a reach is always ZERO.   
        vr=2.5
        mtp=5.0
        iql=0
        for x in range(network['maximum_reach_seqorder'],-1,-1):   
            for head_segment, reach in ordered_reaches[x]:  
                seg_list= seg_list_all[head_segment]  
                ncomp=ncompall[head_segment] 
                iql=iql+1
                mn= 40.0 + iql
                for seg in range(0,ncomp-1):                               
                    segID= seg_list[seg]
                    dx=connections[segID]['data'][supernetwork_data['length_col']]
                    rdmy=0.0
                    for tsi in range (0,ntsi):
                        n=tsi+1.0
                        yintc=0.1
                        #tlf=mtp*np.exp(-0.5*((rdmy-mn)/vr)**2.0)/(vr*(2.0*3.14)**0.5) + 0.1
                        tlf=df.subtools.input_crtor1(mtp, mn, vr, yintc,n)  
                        tlf= tlf/dx
                        qlatral[segID]['qlat'][tsi]=tlf
                        #with open(os.path.join(output_path,"qlat"),'a') as qlat:
                        #    qlat.write("%s %s %s %s %s\n" % ('qlat', tsi, segID, tlf, dx ))                           

               # qlat of fake segment of the last segment of a reach is always ZERO.
                fksegID= fksegID_dsend[head_segment]
                for tsi in range (0,ntsi):
                    qlatral[fksegID]['qlat'][tsi]=0.0           

    return q_bdseg
    return h_bdseg
    return qlatral 
