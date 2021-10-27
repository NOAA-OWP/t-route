import cython
import numpy as np
import pandas as pd
cimport numpy as np
from libc.math cimport isnan, NAN

from .fortran_wrappers cimport c_diffnw
from .. import diffusive_utils as diff_utils
import troute.nhd_network_utilities_v02 as nnu
import troute.nhd_network as nhd_network
from troute.routing.fast_reach.simple_da cimport obs_persist_shift, simple_da_with_decay, simple_da

# TO DO load some example inputs to test the module

@cython.boundscheck(False)
cdef void diffnw(
             double[::1] timestep_ar_g,
             int nts_ql_g,
             int nts_ub_g,
             int nts_db_g,
	     int ntss_ev_g,
             int mxncomp_g,
             int nrch_g,
             double[::1,:] z_ar_g,
             double[::1,:] bo_ar_g,
             double[::1,:] traps_ar_g,
             double[::1,:] tw_ar_g,
             double[::1,:] twcc_ar_g,
             double[::1,:] mann_ar_g,
             double[::1,:] manncc_ar_g,
             double[::1,:] so_ar_g,
             double[::1,:] dx_ar_g,
             double[::1,:] iniq,    
             int frnw_col,
             double[::1,:] frnw_g,
             double[::1,:,:] qlat_g,
             double[::1,:] ubcd_g,
             double[::1] dbcd_g,
             int paradim,
             double[::1] para_ar_g,
             double[:,:,:] out_q,
             double[:,:,:] out_elv):

    cdef:
        double[::1,:,:] q_ev_g = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double, order = 'F')
        double[::1,:,:] elv_ev_g = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double, order = 'F')

    c_diffnw(
            &timestep_ar_g[0],
            &nts_ql_g,
            &nts_ub_g,
            &nts_db_g,
            &ntss_ev_g,
            &mxncomp_g,
            &nrch_g,
            &z_ar_g[0,0],
            &bo_ar_g[0,0],
            &traps_ar_g[0,0],
            &tw_ar_g[0,0],
            &twcc_ar_g[0,0],
            &mann_ar_g[0,0],
            &manncc_ar_g[0,0],
            &so_ar_g[0,0],
            &dx_ar_g[0,0],
            &iniq[0,0],
            &frnw_col,
            &frnw_g[0,0],
            &qlat_g[0,0,0],
            &ubcd_g[0,0],
            &dbcd_g[0],
            &paradim, 
            &para_ar_g[0],	    
            &q_ev_g[0,0,0],
            &elv_ev_g[0,0,0])

    # copy data from Fortran to Python memory view
    out_q[:,:,:] = q_ev_g[::1,:,:]
    out_elv[:,:,:] = elv_ev_g[::1,:,:]

cpdef object compute_diffusive_tst(
    int nsteps,
    float dt,
    int qts_subdivisions,
    list reaches_wTypes, # a list of tuples
    dict rconn,
    const long[:] data_idx,
    object[:] data_cols,
    const float[:,:] data_values,
    const float[:,:] initial_conditions,
    const float[:,:] qlat_values,
    list lake_numbers_col,
    const double[:,:] wbody_cols,
    dict waterbody_parameters,
    const int[:,:] reservoir_types,
    bint reservoir_type_specified,
    str model_start_time,
    const float[:,:] usgs_values,
    const int[:] usgs_positions,
    const int[:] usgs_positions_reach,
    const int[:] usgs_positions_gage,
    const float[:] lastobs_values_init,
    const float[:] time_since_lastobs_init,
    const double da_decay_coefficient,
    dict upstream_results={},
    bint assume_short_ts=False,
    bint return_courant=False,
    dict diffusive_parameters=False,
    ):

    # segment connections dictionary
    connections = nhd_network.reverse_network(rconn)

    # network tailwater
    tw = list(nhd_network.headwaters(rconn))[0]

    # network reaches
    reach_list = []
    for i in reaches_wTypes:
        reach_list.append(i[0])

    # generate diffusive inputs
    diff_inputs = diff_utils.diffusive_input_data_v02(
        tw,
        connections,
        rconn,
        reach_list,
        diffusive_parameters,
        np.asarray(data_cols),
        np.asarray(data_idx),
        np.asarray(data_values),
        np.asarray(qlat_values),
        np.asarray(initial_conditions),
        upstream_results,
        qts_subdivisions,
        nsteps,
        dt,
        np.asarray(lake_numbers_col),
        np.asarray(wbody_cols),
        )

    # unpack/declare diffusive input variables
    cdef:
        double[::1] timestep_ar_g = np.asfortranarray(diff_inputs["timestep_ar_g"])
        int nts_ql_g = diff_inputs["nts_ql_g"]
        int nts_ub_g = diff_inputs["nts_ub_g"]
        int nts_db_g = diff_inputs["nts_db_g"]
        int ntss_ev_g = diff_inputs["ntss_ev_g"]
        int mxncomp_g = diff_inputs["mxncomp_g"]
        int nrch_g = diff_inputs["nrch_g"]
        double[::1,:] z_ar_g = np.asfortranarray(diff_inputs["z_ar_g"])
        double[::1,:] bo_ar_g = np.asfortranarray(diff_inputs["bo_ar_g"])
        double[::1,:] traps_ar_g = np.asfortranarray(diff_inputs["traps_ar_g"])
        double[::1,:] tw_ar_g = np.asfortranarray(diff_inputs["tw_ar_g"])
        double[::1,:] twcc_ar_g = np.asfortranarray(diff_inputs["twcc_ar_g"])
        double[::1,:] mann_ar_g = np.asfortranarray(diff_inputs["mann_ar_g"])
        double[::1,:] manncc_ar_g = np.asfortranarray(diff_inputs["manncc_ar_g"])
        double[::1,:] so_ar_g = np.asfortranarray(diff_inputs["so_ar_g"])
        double[::1,:] dx_ar_g = np.asfortranarray(diff_inputs["dx_ar_g"])
        int frnw_col = diff_inputs["frnw_col"]
        double[::1,:] frnw_g = np.asfortranarray(diff_inputs["frnw_g"], dtype = np.double)
        double[::1,:,:] qlat_g = np.asfortranarray(diff_inputs["qlat_g"])
        double[::1,:] ubcd_g = np.asfortranarray(diff_inputs["ubcd_g"])
        double[::1] dbcd_g = np.asfortranarray(diff_inputs["dbcd_g"])
        int paradim = diff_inputs["paradim"]
        double[::1] para_ar_g = np.asfortranarray(diff_inputs["para_ar_g"])  
        double[::1,:] iniq = np.asfortranarray(diff_inputs["iniq"])
        double[:,:,:] out_q = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)
        double[:,:,:] out_elv = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)

    # call diffusive compute kernel
    diffnw(timestep_ar_g,
           nts_ql_g,
           nts_ub_g,
           nts_db_g,
           ntss_ev_g,
           mxncomp_g,
           nrch_g,
           z_ar_g,
           bo_ar_g,
           traps_ar_g,
           tw_ar_g,
           twcc_ar_g,
           mann_ar_g,
           manncc_ar_g,
           so_ar_g,
           dx_ar_g,
           iniq,
           frnw_col,
           frnw_g,
           qlat_g,
           ubcd_g,
           dbcd_g,
           paradim,
           para_ar_g,
           out_q,
           out_elv)

    # re-format outputs from the diffusive Fortran kernel
    index_array, flowveldepth_unorder = diff_utils.unpack_output(
                                diff_inputs["pynw"],
                                diff_inputs["ordered_reaches"],
                                out_q,
                                out_elv
                                )
        
    # re-index the flowveldepth_unorder array returned by diff_utils
    # TODO return depth not elevation from Tulane model
    flowveldepth_test = np.zeros((data_idx.shape[0], ntss_ev_g*3), dtype='float32')
    flowveldepth_test = (pd.DataFrame(data = flowveldepth_unorder, index = index_array).
                    reindex(index = data_idx).
                    to_numpy(dtype = 'float32'))
    
    cdef float[:,:] flowveldepth = flowveldepth_test
    cdef int gages_size = usgs_positions.shape[0]
    cdef int gage_maxtimestep = usgs_values.shape[1]
    cdef int gage_i, usgs_position_i, timestep
    cdef float[:] lastobs_value, lastobs_time
    cdef float a, da_decay_minutes, da_weighted_shift, replacement_val  # , original_val, lastobs_val,
    cdef (float, float, float, float) da_buf
    cdef int[:] reach_has_gage = np.full(len(reaches_wTypes), np.iinfo(np.int32).min, dtype="int32")
    cdef float[:,:] nudge = np.zeros((1, ntss_ev_g), dtype="float32")
    cdef int qvd_ts_w = 3
    cdef int ts_offset
    cdef float[:,:] flowveldepth_row_fill_buf = np.zeros((1, ntss_ev_g), dtype="float32")
    cdef np.ndarray fill_index_mask = np.ones_like(data_idx, dtype=bool)
    
    lastobs_values = np.array([], dtype="float32")
    lastobs_times = np.array([], dtype="float32")
    if gages_size:
        
        lastobs_times = np.full(gages_size, NAN, dtype="float32")
        lastobs_values = np.full(gages_size, NAN, dtype="float32")

        lastobs_values[0] = lastobs_values_init[0]
        lastobs_times[0] = time_since_lastobs_init[0]
        usgs_positions_i = usgs_positions[0]
                   
        # timestep 0 is the initial condition b/c diffusive wave writes out nsteps + 1 timesteps
        timestep = 0
        while timestep <= ntss_ev_g-1:
            
            ts_offset = timestep * qvd_ts_w

            da_buf = simple_da(
                timestep,
                dt,
                da_decay_coefficient,
                gage_maxtimestep,
                NAN if timestep >= gage_maxtimestep else usgs_values[0,timestep],
                flowveldepth[usgs_positions_i, ts_offset],
                lastobs_times[0],
                lastobs_values[0],
                False,
            )
            
            # fill buffer for all timesteps at gage locations
            flowveldepth_row_fill_buf[0, timestep] = da_buf[0]
            
            # record nudge magnitude and lastobs information
            nudge[0, timestep] = da_buf[1]
            lastobs_times[0] = da_buf[2]
            lastobs_values[0] = da_buf[3]

            timestep += 1
        
        # insert buffer
        flowveldepth[usgs_positions_i,::3] = flowveldepth_row_fill_buf
        
    # create mask array
    for upstream_tw_id in upstream_results:
        tmp = upstream_results[upstream_tw_id]
        fill_index = tmp["position_index"]
        fill_index_mask[fill_index] = False
        
    return np.asarray(data_idx, dtype=np.intp)[fill_index_mask], np.asarray(flowveldepth[:,3:])[fill_index_mask], 0, (np.asarray([data_idx[usgs_position_i] for usgs_position_i in usgs_positions]), np.asarray(lastobs_times), np.asarray(lastobs_values))
