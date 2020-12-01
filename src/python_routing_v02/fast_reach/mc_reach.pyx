# cython: language_level=3, boundscheck=True, wraparound=False, profile=True

import numpy as np
from itertools import chain
from operator import itemgetter
from numpy cimport ndarray
cimport numpy as np
cimport cython
from cython.parallel import prange

#NJF For whatever reason, when cimporting muskingcunge from reach, the linker breaks in weird ways
#the mc_reach.so will have an undefined symbol _muskingcunge, and reach.so will have a ____pyx_f_5reach_muskingcunge
#if you cimport reach, then call explicitly reach.muskingcung, then mc_reach.so maps to the correct module symbol
#____pyx_f_5reach_muskingcunge
#from reach cimport muskingcunge, QVD
cimport reach

@cython.boundscheck(False)
cpdef object binary_find(object arr, object els):
    """
    Find elements in els in arr.
    Args:
        arr: Array to search. Must be sorted
        els:
    Returns:
    """
    cdef long hi = len(arr)
    cdef object idxs = []

    cdef Py_ssize_t L, R, m
    cdef long cand, el
    for el in els:
        L = 0
        R = hi - 1
        m = 0
        while L <= R:
            m = (L + R) // 2
            cand = arr[m]
            if cand < el:
                L = m + 1
            elif cand > el:
                R = m - 1
            else:
                break
        if arr[m] == el:
            idxs.append(m)
        else:
            raise ValueError(f"element {el} not found in {np.asarray(arr)}")
    return idxs


@cython.boundscheck(False)
cdef void compute_reach_kernel(float qup, float quc, int nreach, const float[:,:] input_buf, float[:, :] output_buf, bint assume_short_ts) nogil:
    """
    Kernel to compute reach.
    Input buffer is array matching following description:
    axis 0 is reach
    axis 1 is inputs in th following order:
        qlat, dt, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp
        qup and quc are initial conditions.
    Output buffer matches the same dimsions as input buffer in axis 0
    Input is nxm (n reaches by m variables)
    Ouput is nx3 (n reaches by 3 return values)
        0: current flow, 1: current depth, 2: current velocity
    """
    cdef reach.QVD rv
    cdef reach.QVD *out = &rv

    cdef:
        float dt, qlat, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp
        int i

    for i in range(nreach):
        qlat = input_buf[i, 0] # n x 1
        dt = input_buf[i, 1] # n x 1
        dx = input_buf[i, 2] # n x 1
        bw = input_buf[i, 3]
        tw = input_buf[i, 4]
        twcc =input_buf[i, 5]
        n = input_buf[i, 6]
        ncc = input_buf[i, 7]
        cs = input_buf[i, 8]
        s0 = input_buf[i, 9]
        qdp = input_buf[i, 10]
        velp = input_buf[i, 11]
        depthp = input_buf[i, 12]

        reach.muskingcunge(
                    dt,
                    qup,
                    quc,
                    qdp,
                    qlat,
                    dx,
                    bw,
                    tw,
                    twcc,
                    n,
                    ncc,
                    cs,
                    s0,
                    velp,
                    depthp,
                    out)

#        output_buf[i, 0] = quc = out.qdc # this will ignore short TS assumption at seg-to-set scale?
        output_buf[i, 0] = out.qdc
        output_buf[i, 1] = out.velc
        output_buf[i, 2] = out.depthc
        
        qup = qdp
        
        if assume_short_ts:
            quc = qup
        else:
            quc = out.qdc        

cdef void fill_buffer_column(const Py_ssize_t[:] srows,
    const Py_ssize_t scol,
    const Py_ssize_t[:] drows,
    const Py_ssize_t dcol,
    const float[:, :] src, float[:, ::1] out) nogil:

    cdef Py_ssize_t i
    for i in range(srows.shape[0]):
        out[drows[i], dcol] = src[srows[i], scol]

cpdef object column_mapper(object src_cols):
    """Map source columns to columns expected by algorithm"""
    cdef object index = {}
    cdef object i_label
    for i_label in enumerate(src_cols):
        index[i_label[1]] = i_label[0]

    cdef object rv = []
    cdef object label
    #qlat, dt, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp
    for label in ['dt', 'dx', 'bw', 'tw', 'twcc', 'n', 'ncc', 'cs', 's0']:
        rv.append(index[label])
    return rv


cpdef object compute_network(int nsteps, list reaches, dict connections, 
    const long[:] data_idx, object[:] data_cols, const float[:,:] data_values, 
    const float[:, :] qlat_values, const float[:,:] initial_conditions, 
    # const float[:] wbody_idx, object[:] wbody_cols, const float[:, :] wbody_vals,
    dict upstream_results={},
    bint assume_short_ts=False):
    """
    Compute network
    Args:
        nsteps (int): number of time steps
        reaches (list): List of reaches
        connections (dict): Network
        data_idx (ndarray): a 1D sorted index for data_values
        data_values (ndarray): a 2D array of data inputs (nodes x variables)
        qlats (ndarray): a 2D array of qlat values (nodes x nsteps). The index must be shared with data_values
        initial_conditions (ndarray): an n x 3 array of initial conditions. n = nodes, column 1 = qu0, column 2 = qd0, column 3 = h0
        assume_short_ts (bool): Assume short time steps (quc = qup)
    Notes:
        Array dimensions are checked as a precondition to this method.
        data_idx inc. flowveldepth -- sorted numerically
        Reach_buffer -- sorted topologically
    """
    # Check shapes
    if qlat_values.shape[0] != data_idx.shape[0]:
        raise ValueError(f"Number of rows in Qlat is incorrect: expected ({data_idx.shape[0]}), got ({qlat_values.shape[0]})")
    if qlat_values.shape[1] > nsteps:
        raise ValueError(f"Number of columns (timesteps) in Qlat is incorrect: expected at most ({data_idx.shape[0]}), got ({qlat_values.shape[0]}). The number of columns in Qlat must be equal to or less than the number of routing timesteps")
    if data_values.shape[0] != data_idx.shape[0] or data_values.shape[1] != data_cols.shape[0]:
        raise ValueError(f"data_values shape mismatch")

    # flowveldepth is 2D float array that holds results
    # columns: flow (qdc), velocity (velc), and depth (depthc) for each timestep
    # rows: indexed by data_idx
    cdef float[:,::1] flowveldepth = np.zeros((data_idx.shape[0], nsteps * 3), dtype='float32')

    # Pseudocode: LOOP ON Upstream Inflowers
        # to pre-fill FlowVelDepth
        # fill_index = list_of_all_segments_sorted -- .i.e, data_idx -- .index(upstream_tw_id)
        # # FlowVelDepth[fill_index]['flow'] = UpstreamOutflows[upstream_tw_id]['flow']
        # # FlowVelDepth[fill_index]['depth'] = UpstreamOutflows[upstream_tw_id]['depth']

    # for ts in flowveldepth:
        # print(f"{list(ts)}")
    for upstream_tw_id in upstream_results:
        fill_index = upstream_results[upstream_tw_id]["position_index"]
        # print(f"{upstream_results[upstream_tw_id]['results']}")
        # print(f"filling the {fill_index} row:")
        # print(f"{list(flowveldepth[fill_index])}")
        for idx, val in enumerate(upstream_results[upstream_tw_id]["results"]):
            flowveldepth[fill_index][idx] = val
        # print(f"Now filled, it contains:")
        # print(f"{list(flowveldepth[fill_index])}")

    cdef:
        Py_ssize_t[:] srows  # Source rows indexes
        Py_ssize_t[:] drows_tmp
    
    # Buffers and buffer views
    # These are C-contiguous.
    cdef float[:, ::1] buf, buf_view
    cdef float[:, ::1] out_buf, out_view

    # Source columns
    cdef Py_ssize_t[:] scols = np.array(column_mapper(data_cols), dtype=np.intp)
    
    # hard-coded column. Find a better way to do this
    cdef int buf_cols = 13

    cdef:
        Py_ssize_t i  # Temporary variable
        Py_ssize_t ireach  # current reach index
        Py_ssize_t ireach_cache  # current index of reach cache
        Py_ssize_t iusreach_cache  # current index of upstream reach cache

    # Measure length of all the reaches
    cdef list reach_sizes = list(map(len, reaches))
    # For a given reach, get number of upstream nodes
    cdef list usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    cdef:
        list reach  # Temporary variable
        list bf_results  # Temporary variable

    cdef int reachlen, usreachlen
    cdef Py_ssize_t bidx

    cdef:
        Py_ssize_t[:] reach_cache
        Py_ssize_t[:] usreach_cache

    # reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    # upstream reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)

    ireach_cache = 0
    iusreach_cache = 0
    # copy reaches into an array
    for ireach in range(len(reaches)):
        reachlen = reach_sizes[ireach]
        usreachlen = usreach_sizes[ireach]
        reach = reaches[ireach]

        # set the length (must be negative to indicate reach boundary)
        reach_cache[ireach_cache] = -reachlen
        ireach_cache += 1
        bf_results = binary_find(data_idx, reach)
        for bidx in bf_results:
            reach_cache[ireach_cache] = bidx
            ireach_cache += 1

        usreach_cache[iusreach_cache] = -usreachlen
        iusreach_cache += 1
        if usreachlen > 0:
            for bidx in binary_find(data_idx, connections[reach[0]]):
                usreach_cache[iusreach_cache] = bidx
                iusreach_cache += 1

    cdef int maxreachlen = max(reach_sizes)
    buf = np.empty((maxreachlen, buf_cols), dtype='float32')
    out_buf = np.empty((maxreachlen, 3), dtype='float32')

    drows_tmp = np.arange(maxreachlen, dtype=np.intp)
    cdef Py_ssize_t[:] drows
    cdef float qup, quc
    cdef int timestep = 0
    cdef int ts_offset

# TODO: Split the compute network function so that the part where we set up
# all the indices is separate from the call to loop through them.
# That way, we can refine the functions for preparing the indexes in isolation
# For example, the actual looping function could start about here in the current
# function, and might look like the following Psuedocode
# cpdef(Dataindex):
    # pull indices and put them in arrays
    # minimal validation,
    # Jump straight to nogil.

    with nogil:
        while timestep < nsteps:
            ts_offset = timestep * 3

            ireach_cache = 0
            iusreach_cache = 0
            while ireach_cache < reach_cache.shape[0]:
                reachlen = -reach_cache[ireach_cache]
                usreachlen = -usreach_cache[iusreach_cache]

                ireach_cache += 1
                iusreach_cache += 1
                
                qup = 0.0
                quc = 0.0
                for i in range(usreachlen):
                    
                    '''
                    New logic was added to handle initial conditions:
                    When timestep == 0, the flow from the upstream segments in the previous timestep
                    are equal to the initial conditions. 
                    '''
                        
                    # upstream flow in the current timestep is equal the sum of flows 
                    # in upstream segments, current timestep
                    # Headwater reaches are computed before higher order reaches, so quc can
                    # be evaulated even when the timestep == 0.
                    quc += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset]
                    
                    # upstream flow in the previous timestep is equal to the sum of flows 
                    # in upstream segments, previous timestep
                    if timestep > 0:
                        qup += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset - 3]
                    else:
                        # sum of qd0 (flow out of each segment) over all upstream reaches
                        qup += initial_conditions[usreach_cache[iusreach_cache + i],1]
                        
                buf_view = buf[:reachlen, :]
                out_view = out_buf[:reachlen, :]
                drows = drows_tmp[:reachlen]
                srows = reach_cache[ireach_cache:ireach_cache+reachlen]

                """
                qlat_values may have fewer columns than data_values if qlat data are taken from WRF hydro simulations,
                which are often run at a coarser timestep than routing models. In the fill_buffer_columns call below, 
                the second argument, which defines the column in qlat_values that data should be drawn from, is specified
                such that qlat values are repeated for each of the finer routing timesteps within a WRF hydro timestep. 
                """
                fill_buffer_column(srows, 
                                   int(timestep/(nsteps/qlat_values.shape[1])),  # adjust timestep to WRF-hydro timestep
                                   drows, 
                                   0, 
                                   qlat_values, 
                                   buf_view)
                
                for i in range(scols.shape[0]):
                        fill_buffer_column(srows, scols[i], drows, i + 1, data_values, buf_view)
                        
                # fill buffer with qdp, depthp, velp
                if timestep > 0:
                    fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view)
                    fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view)
                    fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view)
                else:
                    '''
                    Changed made to accomodate initial conditions:
                    when timestep == 0, qdp, and depthp are taken from the initial_conditions array, 
                    using srows to properly index
                    '''
                    for i in range(drows.shape[0]):
                        buf_view[drows[i], 10] = initial_conditions[srows[i],1] #qdp = qd0
                        buf_view[drows[i], 11] = 0.0 # the velp argmument is never used, set to whatever
                        buf_view[drows[i], 12] = initial_conditions[srows[i],2] #hdp = h0

                if assume_short_ts:
                    quc = qup

                compute_reach_kernel(qup, quc, reachlen, buf_view, out_view, assume_short_ts)

                # copy out_buf results back to flowdepthvel
                for i in range(3):
                    fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)

                # Update indexes to point to next reach
                ireach_cache += reachlen
                iusreach_cache += usreachlen
                
            timestep += 1

    return np.asarray(data_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')

#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
cpdef object compute_network_multithread(int nsteps, list reaches, dict connections, 
    const long[:] data_idx, object[:] data_cols, const float[:,:] data_values, 
    const float[:, :] qlat_values, const float[:,:] initial_conditions, 
    const int[:] reach_groups,
    const int[:] reach_group_cache_sizes,
    bint assume_short_ts=False):
    """
    Compute network
    Args:
        nsteps (int): number of time steps
        reaches (list): List of reaches
        connections (dict): Network
        data_idx (ndarray): a 1D sorted index for data_values
        data_values (ndarray): a 2D array of data inputs (nodes x variables)
        qlats (ndarray): a 2D array of qlat values (nodes x nsteps). The index must be shared with data_values
        initial_conditions (ndarray): an n x 3 array of initial conditions. 
        assume_short_ts (bool): Assume short time steps (quc = qup)
    Notes:
        Array dimensions are checked as a precondition to this method.
    """
    # Check shapes
    if qlat_values.shape[0] != data_idx.shape[0]:
        raise ValueError(f"Number of rows in Qlat is incorrect: expected ({data_idx.shape[0]}), got ({qlat_values.shape[0]})")
    if qlat_values.shape[1] > nsteps:
        raise ValueError(f"Number of columns (timesteps) in Qlat is incorrect: expected at most ({data_idx.shape[0]}), got ({qlat_values.shape[0]}). The number of columns in Qlat must be equal to or less than the number of routing timesteps")
    if data_values.shape[0] != data_idx.shape[0] or data_values.shape[1] != data_cols.shape[0]:
        raise ValueError(f"data_values shape mismatch")
    
    cdef float[:,::1] flowveldepth = np.zeros((data_idx.shape[0], nsteps * 3), dtype='float32')
    
    cdef:
        Py_ssize_t[:] srows  # Source rows indexes
        Py_ssize_t[:] drows_tmp
    
    # hard-coded column. Find a better way to do this
    cdef int buf_cols = 13
    
    # Buffers and buffer views
    # These are C-contiguous.
    cdef float[:, ::1] buf, buf_view
    cdef float[:, ::1] out_buf, out_view
    cdef int maxgrouplen = max(reach_group_cache_sizes)
    buf = np.empty((maxgrouplen, buf_cols), dtype='float32')
    out_buf = np.empty((maxgrouplen, 3), dtype='float32')
    
    # Source columns
    cdef Py_ssize_t[:] scols = np.array(column_mapper(data_cols), dtype=np.intp)

    cdef:
        Py_ssize_t i  # Temporary variable
        Py_ssize_t ireach  # current reach index
        Py_ssize_t ireach_cache  # current index of reach cache
        Py_ssize_t iusreach_cache  # current index of upstream reach cache

    # Measure length of all the reaches
    cdef list reach_sizes = list(map(len, reaches))
    # For a given reach, get number of upstream nodes
    cdef list usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    cdef:
        list reach  # Temporary variable
        list bf_results  # Temporary variable

    cdef int reachlen, usreachlen
    cdef int prevreachlen
    cdef Py_ssize_t bidx

    cdef:
        Py_ssize_t[:] reach_cache
        Py_ssize_t[:] usreach_cache
        Py_ssize_t[:] ireach_cache_array
        Py_ssize_t[:] iusreach_cache_array

    # reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    # upstream reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)
    
    # ireach_cache_array
    ireach_cache_array = np.empty(len(reach_sizes), dtype=np.intp)
    iusreach_cache_array = np.empty(len(reach_sizes), dtype=np.intp)

    ireach_cache = 0
    iusreach_cache = 0
    
    ireach_cache_array[0] = 0
    iusreach_cache_array[0] = 0
    # copy reaches into an array
    for ireach in range(len(reaches)):
        reachlen = reach_sizes[ireach]
        usreachlen = usreach_sizes[ireach]
        reach = reaches[ireach]

        # set the length (must be negative to indicate reach boundary)
        reach_cache[ireach_cache] = -reachlen
        ireach_cache += 1
        bf_results = binary_find(data_idx, reach)
        for bidx in bf_results:
            reach_cache[ireach_cache] = bidx
            ireach_cache += 1

        usreach_cache[iusreach_cache] = -usreachlen
        iusreach_cache += 1
        if usreachlen > 0:
            for bidx in binary_find(data_idx, connections[reach[0]]):
                usreach_cache[iusreach_cache] = bidx
                iusreach_cache += 1
                
        if ireach < max(range(len(reaches))):
            ireach_cache_array[ireach+1] = ireach_cache
            iusreach_cache_array[ireach+1] = iusreach_cache

    drows_tmp = np.arange(maxgrouplen, dtype=np.intp)
    
    cdef Py_ssize_t[:] drows
    cdef int timestep = 0
    cdef int ts_offset

    cdef int maxgroupsize = max(reach_groups)
    cdef float[:] qu_buf = np.empty(maxgroupsize, dtype = "float32")
    cdef float[:] quc_view 
    cdef float[:] qup_view
    cdef float quc, qup
    cdef int qu_idx
    cdef int buf_idx
    cdef Py_ssize_t[:] srowsgroup_buf = np.empty(maxgrouplen, dtype = np.intp)
    cdef int srows_idx
    
    cdef:
        Py_ssize_t istart  
        Py_ssize_t iend
        int r
    
    with nogil:
        while timestep < nsteps:
            ts_offset = timestep * 3

            istart = 0
            iend = -1
            for group_i in range(len(reach_group_cache_sizes)):
                
                # index of final reach entry in reach_cache for this group
                iend += reach_groups[group_i]
                
                # prepare group buffers
                buf_view = buf[:reach_group_cache_sizes[group_i],:]
                out_view = out_buf[:reach_group_cache_sizes[group_i],:]
                quc_view = qu_buf[:reach_groups[group_i]]
                qup_view = qu_buf[:reach_groups[group_i]]
                srows = srowsgroup_buf[:reach_group_cache_sizes[group_i]]
                drows = drows_tmp[:reach_group_cache_sizes[group_i]]
                
                # extract upstream flows and populate srows
                qu_idx = 0
                srows_idx = 0
                for r in range(istart,iend+1):
                    
                    ireach_cache = ireach_cache_array[r]
                    iusreach_cache = iusreach_cache_array[r]
                    reachlen = -reach_cache[ireach_cache]
                    usreachlen = -usreach_cache[iusreach_cache]
                    iusreach_cache += 1
                    ireach_cache += 1
                    
                    quc = 0.0
                    qup = 0.0
                    for i in range(usreachlen):
                        quc = quc + flowveldepth[usreach_cache[iusreach_cache + i], ts_offset]
                        if timestep > 0:
                            qup += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset - 3]
                        else:
                            qup += initial_conditions[usreach_cache[iusreach_cache + i],1]
                            
                    quc_view[qu_idx] = quc
                    qup_view[qu_idx] = qup
                    qu_idx += 1
                    
                    # build srows
                    for i in range(reachlen):
                        srows[srows_idx] = reach_cache[ireach_cache + i]
                        srows_idx += 1
                             
                # fill buf_view with qlat, parameter and initial conditions data
                # qlats    
                fill_buffer_column(srows, 
                           int(timestep/(nsteps/qlat_values.shape[1])),
                           drows, 
                           0, 
                           qlat_values, 
                           buf_view)
            
                # parameters
                for i in range(scols.shape[0]):
                    fill_buffer_column(srows, scols[i], drows, i + 1, data_values, buf_view)

                # initial conditions
                if timestep > 0:
                    fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view)
                    fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view)
                    fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view)
                else:
                    for i in range(drows.shape[0]):
                        buf_view[drows[i], 10] = initial_conditions[srows[i],1]
                        buf_view[drows[i], 11] = 0.0
                        buf_view[drows[i], 12] = initial_conditions[srows[i],2]
                
                # ------ !!!!! MULTITHREAD LOOP !!!!! ------ #
                # compute each reach
                for r in prange(istart,iend+1):
                    
                    # reach length for reach r of group
                    ireach_cache = ireach_cache_array[r]
                    reachlen = -reach_cache[ireach_cache]
                    
                    # total reach length of reaches preceeding reach r in group
                    if r > istart:
                        prevreachlen = 0
                        for i in range((r-istart)):
                            prevreachlen = prevreachlen + -reach_cache[ireach_cache_array[(r-(i+1))]]
                    else:
                        prevreachlen = 0
                            
                    # compute reach routing
                    if assume_short_ts:    
                        compute_reach_kernel(qup_view[r-istart], 
                                             quc_view[r-istart], # quc = qup
                                             reachlen, 
                                             buf_view[prevreachlen:prevreachlen+reachlen,:],
                                             out_view[prevreachlen:prevreachlen+reachlen,:],
                                             assume_short_ts)
                        
                    else:
                        compute_reach_kernel(qup_view[r-istart], 
                                             quc_view[r-istart], 
                                             reachlen, 
                                             buf_view[prevreachlen:prevreachlen+reachlen,:],
                                             out_view[prevreachlen:prevreachlen+reachlen,:],
                                             assume_short_ts)
                # END ------ !!!!! MULTITHREAD LOOP !!!!! ------ #
                        
                        
                # place out_view results into flowveldepth
                for i in range(3):
                    fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)
                
                    
                istart = istart + reach_groups[group_i]
                
            timestep += 1

    return np.asarray(data_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')
