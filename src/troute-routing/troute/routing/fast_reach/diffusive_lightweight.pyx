import cython
import numpy as np
cimport numpy as np

from .fortran_wrappers cimport c_compute_diffusive_couplingtimestep

@cython.boundscheck(False)
cdef void diffusive_couplingtimestep(
                                    double[::1] timestep_ar_g,
                                    int nts_ql_g,
                                    int nts_db_g,
                                    int nts_qtrib_g,
                                    int nts_da_g,
                                    int mxncomp_g,
                                    int nrch_g,
                                    double[::1,:] dx_ar_g,
                                    double[::1,:] iniq,
                                    double[::1,:] inidepth,
                                    int frnw_col, 
                                    int[::1,:] frnw_ar_g, 
                                    double[::1,:,:] qlat_g, 
                                    double[::1] dbcd_g,      
                                    double[::1,:] qtrib_g, 
                                    int paradim, 
                                    double[::1] para_ar_g, 
                                    double[::1,:] usgs_da_g, 
                                    int[::1] usgs_da_reach_g,
                                    int nrow_chxsec_lookuptable, 
                                    double[::1,:,:,:] chxsec_lookuptable, 
                                    double[::1,:] z_adj, 
                                    double t_start, 
                                    double t_end,             
                                    double[:,:] out_q_next_out_time, 
                                    double[:,:] out_elv_next_out_time, 
                                    double[:,:] out_depth_next_out_time
):

    cdef:
        double[::1,:] q_next_out_time = np.empty([mxncomp_g, nrch_g], dtype = np.double, order = 'F')
        double[::1,:] elv_next_out_time = np.empty([mxncomp_g, nrch_g], dtype = np.double, order = 'F')
        double[::1,:] depth_next_out_time = np.empty([mxncomp_g, nrch_g], dtype = np.double, order = 'F')

    c_compute_diffusive_couplingtimestep(&timestep_ar_g[0], 
                                         &nts_ql_g, 
                                         &nts_db_g, 
                                         &nts_qtrib_g, 
                                         &nts_da_g,
                                         &mxncomp_g, 
                                         &nrch_g, 
                                         &dx_ar_g[0,0], 
                                         &iniq[0,0], 
                                         &inidepth[0,0], 
                                         &frnw_col, 
                                         &frnw_ar_g[0,0], 
                                         &qlat_g[0,0,0], 
                                         &dbcd_g[0],      
                                         &qtrib_g[0,0], 
                                         &paradim, 
                                         &para_ar_g[0], 
                                         &usgs_da_g[0,0], 
                                         &usgs_da_reach_g[0],
                                         &nrow_chxsec_lookuptable, 
                                         &chxsec_lookuptable[0,0,0,0], 
                                         &z_adj[0,0], 
                                         &t_start, 
                                         &t_end,             
                                         &q_next_out_time[0,0], 
                                         &elv_next_out_time[0,0], 
                                         &depth_next_out_time[0,0]
                                         )

    # copy data from Fortran to Python memory view
    out_q_next_out_time[:,:]     = q_next_out_time[::1,:]
    out_elv_next_out_time[:,:]   = elv_next_out_time[::1,:]
    out_depth_next_out_time[:,:] = depth_next_out_time[::1,:]

cpdef object compute_diffusive_couplingtimestep(
    dict diff_inputs,
    out_chxsec_lookuptable, 
    out_z_adj,
    couplingtime_start,
    couplingtime_end,
):

    # unpack/declare diffusive input variables
    cdef:
        double[::1] timestep_ar_g = np.asfortranarray(diff_inputs['timestep_ar_g'])
        int nts_ql_g = diff_inputs["nts_ql_g"]
        int nts_db_g = diff_inputs["nts_db_g"]
        int nts_qtrib_g = diff_inputs['nts_qtrib_g']
        int nts_da_g = diff_inputs["nts_da_g"]       
        int mxncomp_g = diff_inputs["mxncomp_g"]
        int nrch_g = diff_inputs["nrch_g"]
        double[::1,:] dx_ar_g = np.asfortranarray(diff_inputs["dx_ar_g"])
        double[::1,:] iniq = np.asfortranarray(diff_inputs["iniq"])
        #double[::1,:] inidepth = np.empty([mxncomp_g,nrch_g], dtype = np.double)
        double[::1,:] inidepth = np.asfortranarray(diff_inputs["inidepth"])
        int frnw_col = diff_inputs["frnw_col"]
        int[::1,:] frnw_ar_g = np.asfortranarray(diff_inputs["frnw_g"])
        double[::1,:,:] qlat_g = np.asfortranarray(diff_inputs["qlat_g"])
        double[::1] dbcd_g = np.asfortranarray(diff_inputs["dbcd_g"])
        double[::1,:] qtrib_g = np.asfortranarray(diff_inputs["qtrib_g"])
        int paradim = diff_inputs['paradim']
        double[::1] para_ar_g = np.asfortranarray(diff_inputs["para_ar_g"])
        double[::1,:] usgs_da_g = np.asfortranarray(diff_inputs["usgs_da_g"])   
        int[::1] usgs_da_reach_g = np.asfortranarray(diff_inputs["usgs_da_reach_g"]) 
        int nrow_chxsec_lookuptable = diff_inputs["nrow_chxsec_lookuptable"]        
        double[::1,:,:,:] chxsec_lookuptable = np.asfortranarray(out_chxsec_lookuptable)
        double[::1,:] z_adj = np.asfortranarray(out_z_adj)
        double t_start = couplingtime_start
        double t_end = couplingtime_end
        double[:,:] out_q_next_out_time = np.empty([mxncomp_g,nrch_g], dtype = np.double)
        double[:,:] out_elv_next_out_time = np.empty([mxncomp_g,nrch_g], dtype = np.double)
        double[:,:] out_depth_next_out_time = np.empty([mxncomp_g,nrch_g], dtype = np.double)

    # call diffusive compute kernel
    diffusive_couplingtimestep(
                                timestep_ar_g,
                                nts_ql_g,
                                nts_db_g,
                                nts_qtrib_g,
                                nts_da_g,
                                mxncomp_g,
                                nrch_g,
                                dx_ar_g,
                                iniq,
                                inidepth,
                                frnw_col, 
                                frnw_ar_g, 
                                qlat_g, 
                                dbcd_g,      
                                qtrib_g, 
                                paradim, 
                                para_ar_g, 
                                usgs_da_g, 
                                usgs_da_reach_g,
                                nrow_chxsec_lookuptable, 
                                chxsec_lookuptable, 
                                z_adj, 
                                t_start, 
                                t_end,             
                                out_q_next_out_time, 
                                out_elv_next_out_time, 
                                out_depth_next_out_time
                                )

    return np.asarray(out_q_next_out_time), np.asarray(out_elv_next_out_time), np.asarray(out_depth_next_out_time)