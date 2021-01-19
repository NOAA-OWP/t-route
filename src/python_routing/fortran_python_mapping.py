import os
import numpy as np
import globalvar as g

def fpmap1(
            connections= None 
            , supernetwork_data= None
            , network= None
            , ordered_reaches= None
            , seg_list_all= None
            , ncompall= None
            , nrch_g= None
            , dbfksegID= None
            , z_all= None
            , pynw= None
            , frnw_g= None
            ):    

    #---------------------------------------------------------------------
    #    Mapping between Python head_segment ID and Fortran j                                     
    #---------------------------------------------------------------------                            
    #pynw={}
    #j=-1
    #for x in range(network['maximum_reach_seqorder'],-1,-1):       
    #    for head_segment, reach in ordered_reaches[x]:
    #        j=j+1
    #        pynw[j]=head_segment
            #print(f"x:{x} j:{j} head_segment:{head_segment} pynw:{pynw}")  

    ##-------------------------------------------------------------------------------
    #  Store headwater reach and upstream reaches above a junction
    #  as well as downstream reach after a junction
    #  into python-extension-fortran variables.
    ##-------------------------------------------------------------------------------           
    frj=-1
    for x in range(network['maximum_reach_seqorder'],-1,-1):       
        for head_segment, reach in ordered_reaches[x]:                  
            seg_list= seg_list_all[head_segment]   
            ncomp=ncompall[head_segment] 
            frj=frj+1
            frnw_g[frj,0]= ncomp       

            # ** headwater reach or upstream reaches
            if reach['upstream_reaches']=={0}:
            # headwater reach    
                frnw_g[frj,2]= 0 # the number of upstream reaches
                #print(f"head_segment:{head_segment} usrch_list:{0}")
            else:
            # reaches before a junction
                nusrch= len(reach['upstream_reaches'])                        
                frnw_g[frj,2]= nusrch  # the number of upstream reaches                  
                usrch_list=list(reach['upstream_reaches']) 
                #print(f"head_segment:{head_segment} usrch_list:{usrch_list}")
                i=0
                for usrch in range(0, nusrch):
                    usrchid=usrch_list[usrch] #upstream reach's head segment ID
                    # find Fortran js corresponding to usrchid
                    for j, sid in pynw.items():
                        if sid == usrchid:
                            j1=j
                            i=i+1
                            frnw_g[frj,2+i]=j1  

            # ** downstream reach            
            if (seg_list.count(dbfksegID)>0): 
            # a reach where downstream boundary condition is set.
                frnw_g[frj,1]= -100 # head_segment ID that is in terminal downstream reach.
                                    # That is, -100 indicates the reach of the head segment is 
                                    # terminal downstream reach where ds.bcond. happens.
            else:
            # reach after a junction
                dsrchID= network['reaches'][head_segment]['downstream_reach']                        
                for j, sid in pynw.items():
                    if sid == dsrchID:
                        j1=j                                
                        frnw_g[frj,1]=j1 # J1 is fortran j index equivalent to dsrchID.  
 
    # Adust frnw_g element values according to Fortran-Python index relationship, that is Python i = Fortran i+1
    for frj in range(0,nrch_g):
        frnw_g[frj,1]= frnw_g[frj,1]+1 #downstream reach index for frj reach
        if frnw_g[frj,2]>0:
            nusrch=frnw_g[frj,2]
            for i in range(0, nusrch):
                frnw_g[frj,3+i]= frnw_g[frj,3+i]+1 # upstream reach indicds for frj reach  
    
    output_path='./output2'
    with open(os.path.join(output_path,"frnw"),'a') as frnw:
        for frj in range(0,nrch_g):
            frnw.write("%s %s %s %s %s %s %s %s\n" % (frnw_g[frj,0],frnw_g[frj,1],\
                       frnw_g[frj,2],frnw_g[frj,3],frnw_g[frj,4],frnw_g[frj,5],frnw_g[frj,6],frnw_g[frj,7]))

    
    #    for x in range(network['maximum_reach_seqorder'],-1,-1): 
    #        for head_segment, reach in ordered_reaches[x]:      
    #            if reach['upstream_reaches']=={0}:              
    #                for tsi in range(0,nts_ub):  
    #                    op_ub.write("%s %s %s\n" % (tsi,head_segment,ubcond[head_segment]['discharge'][tsi]))
    
    
    return frnw_g

    
            
