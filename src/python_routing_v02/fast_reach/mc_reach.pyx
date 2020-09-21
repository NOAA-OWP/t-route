# cython: language_level=3, boundscheck=True, wraparound=False, profile=True

import numpy as np
from itertools import chain
from operator import itemgetter
from numpy cimport ndarray
cimport numpy as np
cimport cython

from reach cimport muskingcunge, QVD


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
cdef void compute_reach_kernel(float qup, float quc, int nreach, const float[:,:] input_buf, float[:, :] output_buf) nogil:
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
    cdef QVD rv
    cdef QVD *out = &rv

    cdef:
        float dt, qlat, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp
        int i

    # for each segment in the reach, grab parameter data and call the routing model
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

        # call the routing model
        muskingcunge(
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

        # populate outpub_buf with flow, velocity, and depth result
        output_buf[i, 0] = quc = out.qdc
        output_buf[i, 1] = out.velc
        output_buf[i, 2] = out.depthc

        # set qup to qdp before moving to next downstream segment
        qup = qdp

cdef void fill_buffer_column(const Py_ssize_t[:] srows,
    const Py_ssize_t scol,
    const Py_ssize_t[:] drows,
    const Py_ssize_t dcol,
    const float[:, :] src, float[:, ::1] out) nogil:
    """
    Parameters
    ----------
    srows: 
    scol: column index for specified timestep
    drows:
    dcol: column index of buf_view to be edited
    src: data set to grab data from
    out: array to edit and pass back
    """ 

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
<<<<<<< HEAD:src/python_routing_v02/fast_reach/mc_reach.pyx
<<<<<<< HEAD:src/python_routing_v02/fast_reach/mc_reach.pyx
    const float[:, :] qlat_values,
    # const float[:] wbody_idx, object[:] wbody_cols, const float[:, :] wbody_vals,
=======
    const float[:, :] qlat_values, cost float[:,:] initial_conditions, 
>>>>>>> new argument for initial conditions passed, change first three cols of flowveldepth:src/python_routing_v02/mc_reach.pyx
=======
    const float[:, :] qlat_values, const float[:,:] initial_conditions, 
>>>>>>> warm start (initial conditions) allowed):src/python_routing_v02/mc_reach.pyx
    bint assume_short_ts=False):
    """
    Compute network

    Args:
        nsteps (int): 
            number of time steps
        reaches (list): 
            List of reaches
        connections (dict): 
            Network
        data_idx (ndarray): 
            a 1D sorted index for data_values
        data_values (ndarray): 
            a 2D array of data inputs (nodes x variables)
        qlats (ndarray): 
            a 2D array of qlat values (nodes x nsteps of WRF simulation). The index must be shared with data_values.
            NOTE: WRF simulations are often at a coarser timestep than the routing simulation.
        initial_conditions (ndarray): 
            a 2D array of initial states (nodes x [qu0, qd0, h0]). The index must be shared with data_values
                qu0 = initial upstream flow
                qd0 = initial downstream flow
                h0 = initial depth
        assume_short_ts (bool): Assume short time steps (quc = qup)

    Notes:
        Array dimensions are checked as a precondition to this method.
    """
    
    # Check shapes
    """
    EDIT
    Changed the qlateral dim check to allow the number of timesteps (# of columns) to be different than nsteps
    """
    if qlat_values.shape[0] != data_idx.shape[0]:
        raise ValueError(f"Number of rows in Qlat is incorrect: expected ({data_idx.shape[0]}), got ({qlat_values.shape[0]})")
        
    if data_values.shape[0] != data_idx.shape[0] or data_values.shape[1] != data_cols.shape[0]:
        raise ValueError(f"data_values shape mismatch")

    # flowveldepth is 2D float array that holds results
    # columns: flow (qdc), velocity (velc), and depth (depthc) for each timestep
    # rows: indexed by data_idx
    cdef float[:,::1] flowveldepth = np.zeros((data_idx.shape[0], nsteps * 3), dtype='float32')

    cdef:
        Py_ssize_t[:] srows  # Source rows indexes
        Py_ssize_t[:] drows_tmp
        Py_ssize_t[:] usrows # Upstream row indexes 
    
    # Buffers and buffer views
    # These are C-contiguous.
    cdef float[:, ::1] buf, buf_view
    cdef float[:, ::1] out_buf, out_view

    # Source columns
    cdef Py_ssize_t[:] scols = np.array(column_mapper(data_cols), dtype=np.intp)
    
    # hard-coded column. Find a better way to do this
    cdef int buf_cols = 13

<<<<<<< HEAD:src/python_routing_v02/fast_reach/mc_reach.pyx
    cdef:
        Py_ssize_t i  # Temporary variable
        Py_ssize_t ireach  # current reach index
        Py_ssize_t ireach_cache  # current index of reach cache
        Py_ssize_t iusreach_cache  # current index of upstream reach cache

    # Measure length of all the reaches
    cdef list reach_sizes = list(map(len, reaches))
    # For a given reach, get number of upstream nodes
=======
    # reach_sizes - a list of the number of segments in each reach
    cdef list reach_sizes = list(map(len, reaches))
    
    # usreach_sizes - a list of the number of upstream reaches feeding into each reach
    # usreach_size = 0, if reach is a headwater, usreach_size = 2 if reach is below a 2-reach junction
>>>>>>> warm start (initial conditions) allowed):src/python_routing_v02/mc_reach.pyx
    cdef list usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    cdef:
        list reach  # Temporary variable
        list bf_results  # Temporary variable

    cdef int reachlen, usreachlen
    cdef Py_ssize_t bidx
    cdef list buf_cache = []
<<<<<<< HEAD:src/python_routing_v02/fast_reach/mc_reach.pyx

    cdef:
        Py_ssize_t[:] reach_cache
        Py_ssize_t[:] usreach_cache

    # reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    # upstream reach cache is ordered 1D view of reaches
    # [-len, item, item, item, -len, item, item, -len, item, item, ...]
    usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)

=======
    
    # reach_cache and usreach_cache are intializes with np.empty, which creates an array of specified dimensions
    # filled with seemingly random numbers, which are later re-set to specified values. 
    cdef Py_ssize_t[:] reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    cdef Py_ssize_t[:] usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)
    
    # number of timesteps qlat dataframe, from wrf-hydro simulation
    cdef int ntsteps_wrf = qlat_values.shape[1]
    
    
    
    '''
    ***** Populate reach_cache arrays. ******
    cache arrays have length = (sum(reach_sizes) + len(reach_sizes)). 
    Negative values are used to indicate reach boundaries and the number of segments in the reach
    Positive values following (or between) negative values contain indexes for each segment in the reach
    For example: If a reach has 2 segments, and the index values of those segments are 4 and 5, then
    the reach cache array would read [... -2 4 5 ...] for that particular reach.
    
    reach_cache arrays are needed to properly index data that is passed to the muskingumcunge kernel
    '''
    # initialize itterators, reach and usreach
>>>>>>> warm start (initial conditions) allowed):src/python_routing_v02/mc_reach.pyx
    ireach_cache = 0
    iusreach_cache = 0 
    
    # loop through each reach in the network
    for ireach in range(len(reaches)):
        
        reachlen = reach_sizes[ireach] # number of segments in the reach
        usreachlen = usreach_sizes[ireach] # number of segments in the upstream reach
        reach = reaches[ireach] # segment IDs in the reach

        # set the length (must be negative to indicate reach boundary)
        reach_cache[ireach_cache] = -reachlen
        ireach_cache += 1
        
        # find the indices positions associated with the segments in this reach
        bf_results = binary_find(data_idx, reach)
        
        # update reach_cache values with segment indices
        for bidx in bf_results:
            reach_cache[ireach_cache] = bidx
            ireach_cache += 1

        # set the length (must be negative to indicate reach boundary)
        usreach_cache[iusreach_cache] = -usreachlen
        iusreach_cache += 1
        
        # if the reach is NOT a headwater (i.e. usreachlen > 0)
        if usreachlen > 0:
            
            # find the index position associated with segments in the upstream reach
            for bidx in binary_find(data_idx, connections[reach[0]]):
                
                # update usreach_cache values with segment indices
                usreach_cache[iusreach_cache] = bidx
                iusreach_cache += 1

    # initialize buf and out_buf variables using np.empty
    # buf = n x m array, n = number of segments in largest reach, m = 13
    # outbuf = n x 3 array, n = number of segments in the largest reach, 3 b/c model outputs 3 vars: flow, depth, vel
    cdef int maxreachlen = max(reach_sizes)
    buf = np.empty((maxreachlen, buf_cols), dtype='float32')
    out_buf = np.empty((maxreachlen, 3), dtype='float32')

    
    drows_tmp = np.arange(maxreachlen, dtype=np.intp)
    cdef Py_ssize_t[:] drows
    cdef float qup, quc
    cdef int timestep = 0
    cdef int ts_offset
    
    with nogil:
        
        # loop through timesteps
        while timestep < nsteps:
            
            # specify ts_offset as 3 b/c there are three variables we are tracking: flow, velocity, depth
            # why is this inside the time loop? 
            ts_offset = timestep * 3

            # initialize reach cache index
            ireach_cache = 0
            iusreach_cache = 0
            
            # step through each element of reach_cache (number of network segs + number of network reaches)
            while ireach_cache < reach_cache.shape[0]:
               
                # grab the reach length (i.e. number of segments)
                reachlen = -reach_cache[ireach_cache]
                # grab the upstream reach length (i.e. number of segments in the upstream reach)
                usreachlen = -usreach_cache[iusreach_cache]

                # update the reach_cache itterators by 1
                ireach_cache += 1
                iusreach_cache += 1

                # initialize qup and quc as zero
                qup = 0.0
                quc = 0.0
                
                '''
                flows from the terminal segments of ustream reaches (* in schematic below) are 
                upstream boundary conditions for the first segment
                '''
                                # \    /
                                #  \  /
                                #  *\/*
                                #    \
                                #     \
                
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
                 
                # create buf_view and out_view arrays, which contain as many rows as there are segments in the reach
                buf_view = buf[:reachlen, :]
                out_view = out_buf[:reachlen, :]
                
                # create drows array, which contains simple 0 - n index for each segment in the reach
                drows = drows_tmp[:reachlen]
                
                # create the srows array, which clips indexing data from the reach_cache array for this reach
                # use srows to properly index data from initial_conditions
                srows = reach_cache[ireach_cache:ireach_cache+reachlen]
                
                """
                How to handle lateral inflows from WRF-Hydro?
                ---------------------------------------------
                - Lateral inflows from WRF-Hydro are at a shorter timestep than the routing model
                - Need to take the same qlateral value several routing timesteps, before switching to next
                - Lets assume that the total simulation time is equal between WRF and routing sims
                    - if follows that the number of timesteps in the routing simulation, divided by the number
                      of timesteps in the WRF simulation, equals the number of routing timesteps per WRF timestep.
                - For example, there are 100 timesteps in the routing simulation and 10 timesteps in the WRF simulation
                    - Then, there are 10 routing timesteps per WRF timestep. 
                - Need to change src column indexing from scol to     
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
                
                # if NOT initial timestep populate qdp, velp and depthp with simulation results from previous timestep
                if timestep > 0:
                    fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view) # qdp
                    fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view) # velp
                    fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view) # depthp
                
                # if initial timestep, then qdp, velp, and depthp = 0
                else:
                    # fill buffer with constant
                    '''
                    Changed made to accomodate initial conditions:
                    when timestep == 0, qdp, and depthp are taken from the initial_conditions array, 
                    using srows to properly index
                    '''
                    for i in range(drows.shape[0]):
                        buf_view[drows[i], 10] = initial_conditions[srows[i],1] # qdp
                        buf_view[drows[i], 11] = 0 # velp - zero b/c it isn't used in the MC routing model
                        buf_view[drows[i], 12] = initial_conditions[srows[i],2] # depthp

                if assume_short_ts:
                    quc = qup
                
                # call function to compute reach routing for a single timestep
                compute_reach_kernel(qup, quc, reachlen, buf_view, out_view)

                # copy out_buf results back to flowdepthvel
                for i in range(3):
                    fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)

                # Update indexes to point to next reach
                ireach_cache += reachlen
                iusreach_cache += usreachlen
                
            timestep += 1
    return np.asarray(data_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')
