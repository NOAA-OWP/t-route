import os
import sys
import time
import glob
import nhd_io as nio

#-------------------------------------------------------------------
# The unit of qlateral pulled from .CHRTOUT_DOMAIN1 is m^3/sec.
# As diffusive routing takes the unit of m^2/sece, so it should be
# coverted as such here.
#-------------------------------------------------------------------
def retrieve_qlatral(        
         connections= None  
        , supernetwork_data= None
        , network= None
        , ordered_reaches= None
        , seg_list_all= None
        , ncompall= None 
        , fksegID_dsend= None
        , nts_ql_g= None
        , custom_input_folder = None
        , custom_input_file= None
        , custom_input_file_path= None
        , output_path= None
        , qlatral= None 
        ):
    
    supernetwork_parameters = None
    waterbody_parameters = None
    if custom_input_file:
        (
            supernetwork_parameters,
            waterbody_parameters,
            forcing_parameters,
            restart_parameters,
            output_parameters,
            run_parameters,
        ) = nio.read_custom_input(os.path.join(custom_input_folder,custom_input_file))
        break_network_at_waterbodies = run_parameters.get(
            "break_network_at_waterbodies", None
        )

        verbose = run_parameters.get("verbose", None)
        showtiming = run_parameters.get("showtiming", None)

        qlat_const = forcing_parameters.get("qlat_const", None)
        qlat_input_file = forcing_parameters.get("qlat_input_file", None)
        qlat_input_folder = forcing_parameters.get("qlat_input_folder", None)
        qlat_file_pattern_filter = forcing_parameters.get("qlat_file_pattern_filter", None)
        qlat_file_index_col = forcing_parameters.get("qlat_file_index_col", None)
        qlat_file_value_col = forcing_parameters.get("qlat_file_value_col", None)

    # STEP 5: Read (or set) QLateral Inputs
    if showtiming:
        start_time = time.time()
    if verbose:
        print("creating qlateral array ...")

    # initialize qlateral dict
    qlateral0 = {}

    if qlat_input_folder:
        qlat_files = []
        for pattern in qlat_file_pattern_filter:
            qlat_files.extend(glob.glob(qlat_input_folder + pattern))
        qlat_df = nio.get_ql_from_wrf_hydro(
            qlat_files=qlat_files,
            index_col=qlat_file_index_col,
            value_col=qlat_file_value_col,
        )

    elif qlat_input_file:
        qlat_df = nio.get_ql_from_csv(qlat_input_file)

    else:
        qlat_df = pd.DataFrame(
            qlat_const, index=connections.keys(), columns=range(nts), dtype="float32"
        )
    
    ql_index=[]
    for index, row in qlat_df.iterrows():
        qlateral0[index] = row.tolist()
        ql_index.append(index)
        
    if verbose:
        print("qlateral array complete")

    if showtiming:
        print("... in %s seconds." % (time.time() - start_time))
        start_time = time.time()

    # pass retrieved lateral flow to qlatral dictionary
    for x in range(network['maximum_reach_seqorder'],-1,-1):   
        for head_segment, reach in ordered_reaches[x]:  
            seg_list= seg_list_all[head_segment]  
            ncomp=ncompall[head_segment] 
            for seg in range(0,ncomp-1):                               
                segID= seg_list[seg]
                if ql_index.count(segID)>0:
                    dx=connections[segID]['data'][supernetwork_data['length_col']]
                    qlatral[segID]['qlat']=qlateral0[segID]
                    # convert unit from m^3/s to m^2/s
                    for tsi in range (0,nts_ql_g):
                        tlf= qlatral[segID]['qlat'][tsi]
                        tlf= tlf/dx
                        qlatral[segID]['qlat'][tsi]= tlf

            # qlat of fake segment of the last segment of a reach is always -1.
            fksegID= fksegID_dsend[head_segment]
            for tsi in range (0,nts_ql_g):
                qlatral[fksegID]['qlat'][tsi]=-1.0        
        
        
    #with open(os.path.join(output_path,"qlat"),'a') as op_ql:
    #    for index, row in qlat_df.iterrows():
    #        op_ql.write("%s %s\n" % (index, qlateral0[index])) 
 
    with open(os.path.join(output_path,"qlat_wrfhydro"),'a') as op_ql2:
        for x in range(network['maximum_reach_seqorder'],-1,-1):   
            for head_segment, reach in ordered_reaches[x]:  
                seg_list= seg_list_all[head_segment]  
                ncomp=ncompall[head_segment] 
                for seg in range(0,ncomp-1):                               
                    segID= seg_list[seg]
                    dx=connections[segID]['data'][supernetwork_data['length_col']]
                    for tsi in range (0,nts_ql_g):                    
                        op_ql2.write("%s %s %s\n" % (segID,  qlatral[segID]['qlat'][tsi], dx)) 
    
    return qlatral 