import os
import numpy as np
import statistics as stat
import dfcnhi as df

#--------------------------------------------------------------------------------------------------------
# This file constains key components of diffusive model using Crank-Nicolson Hermite Interpolation method 

#--------------------------------------------------------------------------------------------------------


# ****
# ** set initial condition of discharge q for each segment of every reach
# ****
def icond_q(          
            connections= None 
            , supernetwork_data= None
            , network= None
            , ordered_reaches= None
            , fksegID_dsend= None
            , seg_list_all= None
            , ncompall= None 
            , z_all= None
            , q_bdseg= None
            , ydaq= None
            , qlatral= None
            ):   
    
    for x in range(network['maximum_reach_seqorder'],-1,-1):       
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= seg_list_all[head_segment]   
            ncomp=ncompall[head_segment] 
            q=np.zeros(ncomp)

            if reach['upstream_reaches']=={0}:
            # headbasin reach    
                # At headbasin segment, a known discharge is used. 
                q[0]= q_bdseg[head_segment]['known discharge'][0]
            else:
            # reach after a junction
                # add qs of the last segments of all upstream reaches
                qjt=0.0
                usrch_list=list(reach['upstream_reaches'])
                usrchnb= len(reach['upstream_reaches'])                        
                for usrch in range(0,usrchnb):                    
                    uslsegID= fksegID_dsend[usrch_list[usrch]]  #uslsegID= last_segment_reach[usrch_list[usrch]] 
                    qjt= qjt + ydaq[uslsegID]['q'][0]
                # cascade the accumulated q at a junction down to the downream segments
                q[0]= qjt           
            
            segID0= seg_list[0] 
            ydaq[segID0]['q'][0]= q[0]

            for seg in range(1,ncomp):
                segID= seg_list[seg-1]
                dx=connections[segID]['data'][supernetwork_data['length_col']]
                q[seg]= q[0] + qlatral[segID]['qlat'][0]*dx
                
                segID1= seg_list[seg]      
                ydaq[segID1]['q'][0]= q[seg]
    
    return ydaq


# ****
#  ** As moving reach by reach from upstream to downstream, compute the following 
# for each reach:      
#  1) E,F,Ex,Fx from upper end node to lower end node of a given reach, p146,rm4.
#  2) q and qpx at ts+1 from lower end node to upper end node of a given reach.  
# ****
def dfcnhi_qqpx(          
                connections= None 
                , supernetwork_data= None
                , network= None
                , ordered_reaches= None
                , fksegID_dsend= None
                , seg_list_all= None
                , ncompall= None 
                , ntsi= None
                , dtini = None
                , tc= None
                , ts= None
                , theta= None
                , tarri= None 
                , q_bdseg= None
                , ydaq= None
                , dfpara= None
                , qlatral= None
                , output_path= None
                ):       

    df.var.theta= theta
    arri=np.zeros(ntsi)
    
    for x in range(network['maximum_reach_seqorder'],-1,-1):       
        for head_segment, reach in ordered_reaches[x]:                  
            
            seg_list= seg_list_all[head_segment] 
            ncomp=ncompall[head_segment] 
            
            df.var.ncomp=ncomp
            df.var.celty=np.zeros(ncomp)
            df.var.diffty=np.zeros(ncomp)
            df.var.q=np.zeros(ncomp)
            df.var.qpx=np.zeros(ncomp)
            df.var.dx=np.zeros(ncomp)
            df.var.qlatj=np.zeros(ncomp)
            df.var.dtini=dtini

            # input 1: celerity, diffusivity, qpx at current time ts or tc
            if ts==0:
            # initial time:
                dfv_ini=10.0
                for seg in range(0,ncomp):
                    df.var.celty[seg]=1.0
                    df.var.diffty[seg]= dfv_ini
                    df.var.qpx[seg]=0.0
                for seg in range(0,ncomp):
                    segID= seg_list[seg]
                    dfpara[segID]['celerity'][0]=1.0
                    dfpara[segID]['diffusivity'][0]= dfv_ini
                    dfpara[segID]['qpx'][0]= 0.0
            else:
                for seg in range(0,ncomp):
                    segID= seg_list[seg]
                    # celerity,diffusivity, and qpx at current step, ts or tc                                
                    df.var.celty[seg]=dfpara[segID]['celerity'][0] 
                    df.var.diffty[seg]=dfpara[segID]['diffusivity'][0]                                                                
                    df.var.qpx[seg]=dfpara[segID]['qpx'][0]

            # input 2: q, qlat at current time ts or tc & dx
            for seg in range(0,ncomp):                                
                if seg==ncomp-1:
                    #segID= seg_list[seg-1]
                    df.var.dx[seg]=0.0
                else:
                    segID= seg_list[seg]  
                    df.var.dx[seg]=connections[segID]['data'][supernetwork_data['length_col']] 
                
                segID1= seg_list[seg]
                df.var.q[seg]= ydaq[segID1]['q'][ts]
                # estimate lateral flow by interpolating given flow data points
                for tsi in range (0,ntsi):
                    arri[tsi]=qlatral[segID1]['qlat'][tsi] 
                df.var.qlatj[seg]=np.interp(tc,tarri,arri)

            # Interpolate qlat at the last segment before the fake segment at ts+1
            segID2= seg_list[ncomp-2]
            for tsi in range (0,ntsi):                
                arri[tsi]=qlatral[segID2]['qlat'][tsi] 
            tf0= tc + dtini/60.0 
            df.var.qlatf=np.interp(tf0,tarri,arri)

            # call Fortran subroutine ef_calc for outputing q(i), qpx(i) at i=0,ncomp-1 and at ts+1
            df.diff.ef_calc()

            for seg in range(0,ncomp):     
                segID3= seg_list[seg]
                ydaq[segID3]['q'][ts+1]= df.var.q[seg]
                dfpara[segID3]['qpx'][0]= df.var.qpx[seg]

            # update q at the first segment of a reach
            if reach['upstream_reaches']=={0}:
            # At headbasin reach, the first segment takes interpolated discharge value at ts+1    
                for tsi in range (0,ntsi):
                    arri[tsi]=q_bdseg[head_segment]['known discharge'][tsi]
                tf0= tc + dtini/60.0 
                ydaq[head_segment]['q'][ts+1]= np.interp(tf0,tarri,arri)  
            else:
            # At a junction, the first segment take discharge by 
            # adding qs of the last segments of all upstream reaches
                qjt=0.0
                usrch_list=list(reach['upstream_reaches'])
                usrchnb= len(reach['upstream_reaches'])                        
                for usrch in range(0,usrchnb):
                    uslsegID= fksegID_dsend[usrch_list[usrch]] 
                    qjt= qjt + ydaq[uslsegID]['q'][ts+1]
                # cascade the accumulated q at a junction down to the downream segments
                ydaq[head_segment]['q'][ts+1]=qjt      

            with open(os.path.join(output_path,"qqpx"),'a') as qqpx:
                for seg in range(0,ncomp):
                    segID= seg_list[seg]
                    qqpx.write("%s %s %s %s %s %s\n" %\
                               ('qqpx',ts+1, tc+dtini/60.0, segID, ydaq[segID]['q'][ts+1], dfpara[segID]['qpx'][0]))
                #qqpx.write("\n")
    return ydaq
    return dfpara



# **** 
# **As moving reach by reach from downstream to upstream, compute the following 
# for each reach:           

#  1) y, celerity, diffusivity at ts+1 from lower end node to upper end node of 
#  a given reach, p149,rm4.    
# ****
def dfcnhi_elv(          
                connections= None 
                , supernetwork_data= None
                , network= None
                , ordered_reaches= None
                , seg_list_all= None
                , ncompall= None
                , z_all= None
                , ntsi= None
                , dtini = None
                , tc= None
                , ts= None
                , tarri= None 
                , tzeq_flag= None
                , y_opt= None
                , dbfksegID= None
                , h_bdseg= None
                , nel= None
                , xsec_attr_all= None
                , nhincr_m= None
                , nhincr_f= None
                , uflt= None
                , ydaq= None
                , dfpara= None
                , cel_av= None
                , output_path= None
                ):       
    
    rch=0 # count the number of reaches    
    arri=np.zeros(ntsi)
    
    for x in range(0, network['maximum_reach_seqorder']+1):  
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= seg_list_all[head_segment] 
            ncomp=ncompall[head_segment]                   

            df.var.ncomp=ncomp

            df.var.elv=np.zeros(ncomp)
            df.var.celty=np.zeros(ncomp)
            df.var.diffty=np.zeros(ncomp)
            df.var.q=np.zeros(ncomp)  

            df.var.z_ar=np.zeros(ncomp)
            df.var.bo_ar=np.zeros(ncomp)
            df.var.traps_ar=np.zeros(ncomp)
            df.var.tw_ar=np.zeros(ncomp)
            df.var.twcc_ar=np.zeros(ncomp)
            df.var.mann_ar=np.zeros(ncomp)
            df.var.manncc_ar=np.zeros(ncomp)
            df.var.dx=np.zeros(ncomp)

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
                df.var.dx[seg]=connections[segID]['data'][supernetwork_data['length_col']]
                segID1= seg_list[seg]  
                df.var.z_ar[seg]= z_all[segID1]['adj.alt'][0]

            # y downstream boundary value at time ts+1                  
            if (seg_list.count(dbfksegID)>0): 
            # At downstream terminal boundary, use known depth at the boundary segment
                # interpolate y with the given depth data
                for tsi in range (0,ntsi):                               
                    arri[tsi]=h_bdseg[dbfksegID]['known depth'][tsi]
                tf0= tc+dtini/60.0
                elv_fkseg= np.interp(tf0,tarri,arri)  #tc [min], tarri [min]
                elv_fkseg= elv_fkseg + z_all[dbfksegID]['adj.alt'][0]
                df.var.elv[ncomp-1]= elv_fkseg      
            else:
            # At downstream junction boundary, the depth of fake segment of an upstream reach
            # is equal to that of the first segment of the downstream reach
                # the first segment id of downstream reach from a junction
                dsrchID= network['reaches'][head_segment]['downstream_reach']                        
                elv_fkseg= ydaq[dsrchID]['elev'][ts+1] 
                df.var.elv[ncomp-1]= elv_fkseg                         

            df.var.tzeq_flag= tzeq_flag
            df.var.y_opt= y_opt  # 1 for normal depth(kinematic) by lookup table 
                                 # 2 for depth of diffusive wave.
            if y_opt==1:
                if tzeq_flag==0:
                    df.var.xsec_attr_rch=np.zeros((4,nel,ncomp))  
                    for seg in range(0,ncomp):                  
                        segID= seg_list[seg]
                        df.var.q[seg]= ydaq[segID]['q'][ts+1]
                        for i1 in range(0,nel):
                            df.var.xsec_attr_rch[0,i1,seg]= xsec_attr_all[segID]['elev'][i1]
                            df.var.xsec_attr_rch[1,i1,seg]= xsec_attr_all[segID]['convey'][i1]
                            df.var.xsec_attr_rch[2,i1,seg]= xsec_attr_all[segID]['topwd'][i1]
                            df.var.xsec_attr_rch[3,i1,seg]= xsec_attr_all[segID]['uniflow'][i1]
                
                elif tzeq_flag==1:
                    df.var.ufqlt_m= np.zeros((nhincr_m,ncomp))                        
                    df.var.ufhlt_m= np.zeros((nhincr_m,ncomp))
                    df.var.ufqlt_f= np.zeros((nhincr_f,ncomp))
                    df.var.ufhlt_f= np.zeros((nhincr_f,ncomp))
                    for seg in range(0,ncomp):                  
                        segID= seg_list[seg]
                        df.var.q[seg]= ydaq[segID]['q'][ts+1]
                        #if df.var.y_opt==1:
                        for i1 in range(0,nhincr_m):
                            df.var.ufqlt_m[i1,seg]=uflt[segID]['q_m'][i1]
                            df.var.ufhlt_m[i1,seg]=uflt[segID]['h_m'][i1]
                        for i1 in range(0,nhincr_f):
                            df.var.ufqlt_f[i1,seg]=uflt[segID]['q_f'][i1] 
                            df.var.ufhlt_f[i1,seg]=uflt[segID]['h_f'][i1]                

            df.diff.elv_calc()

            for seg in range(0,ncomp):                    
                segID= seg_list[seg]
                ydaq[segID]['elev'][ts+1]= df.var.elv[seg]
                dfpara[segID]['celerity'][0]=df.var.celty[seg]
                dfpara[segID]['diffusivity'][0]= df.var.diffty[seg]

            cel_av[rch]= dfpara[head_segment]['celerity'][0]
            rch= rch+1

            with open(os.path.join(output_path,"elev"),'a') as elev:
                for seg in range(0,ncomp):
                    segID= seg_list[seg]
                    elev.write("%s %s %s %s %s %s %s\n" %\
                        ('elev',ts+1, tc+dtini/60.0, segID, ydaq[segID]['elev'][ts+1],\
                        dfpara[segID]['celerity'][0],dfpara[segID]['diffusivity'][0]))
            #elev.write("\n")  
    
    return ydaq
    return dfpara       
    return cel_av     
    