import mc_reach
import numpy as np
import utils

def compute_network(nsteps, reaches, connections, 
    parameter_idx, parameter_cols, parameter_values, 
    qlat_values,
    reach_groups,
    reach_group_cache_sizes,
    assume_short_ts=False,
):
    """
    Compute network

    Args:
        nsteps (int): number of time steps
        reaches (list): List of reaches
    with gil:
        connections (dict): Network
        parameter_idx (ndarray): a 1D sorted index for parameter_values
        parameter_values (ndarray): a 2D array of data inputs (nodes x variables)
        qlats (ndarray): a 2D array of qlat values (nodes x nsteps). The index must be shared with parameter_values
        assume_short_ts (bool): Assume short time steps (quc = qup)
        reach_groups:
        reach_group_cache_sizes:

    Notes:
        Array dimensions are checked as a precondition to this method.
    """
    # Check shapes
    if qlat_values.shape[0] != parameter_idx.shape[0] or qlat_values.shape[1] != nsteps:
        raise ValueError(f"Qlat shape is incorrect: expected ({parameter_idx.shape[0], nsteps}), got ({qlat_values.shape[0], qlat_values.shape[1]})")
    if parameter_values.shape[0] != parameter_idx.shape[0] or parameter_values.shape[1] != parameter_cols.shape[0]:
        raise ValueError(f"parameter_values shape mismatch")

    # flowveldepth is 2D float array that holds results
    # columns: flow (qdc), velocity (velc), and depth (depthc) for each timestep
    # rows: indexed by parameter_idx
    flowveldepth = np.zeros((parameter_idx.shape[0], nsteps * 3), dtype='float32')

    #cdef:
    #    Py_ssize_t[:] srows  # Source rows indexes
    #    Py_ssize_t[:] drows_tmp
    #    Py_ssize_t[:] usrows # Upstream row indexes 
    
    # Buffers and buffer views
    # These are C-contiguous.
    #cdef float[:, ::1] buf, buf_view
    #cdef float[:, ::1] out_buf, out_view

    # Source columns
    scols = np.array(mc_reach.column_mapper(parameter_cols), dtype=np.intp)
    
    # hard-coded column. Find a better way to do this
    buf_cols = 13

    #cdef:
    #    Py_ssize_t i  # Temporary variable
    #    Py_ssize_t ireach  # current reach index
    #    Py_ssize_t ireach_cache  # current index of reach cache
    #    Py_ssize_t ireach_cache_end  # end index of reach cache
    #    Py_ssize_t iusreach_cache  # current index of upstream reach cache

    # Measure length of all the reaches
    reach_sizes = list(map(len, reaches))
    # For a given reach, get number of upstream nodes
    # cdef list usreach_sizes = [0 for reach in reaches]
    usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    #cdef:
    #    list reach  # Temporary variable
    #    list bf_results  # Temporary variable

    #cdef int reachlen, usreachlen
    #cdef Py_ssize_t bidx
    buf_cache = []

    #cdef:
    #    Py_ssize_t[:] reach_cache
    #    Py_ssize_t[:] usreach_cache

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
        print(f"looping: usreachlen {usreachlen}")
        reach = reaches[ireach]
        print(f"looping: reach {reach}")

        # set the length (must be negative to indicate reach boundary)
        reach_cache[ireach_cache] = -reachlen
        ireach_cache += 1
        bf_results = mc_reach.binary_find(parameter_idx, reach)
        for bidx in bf_results:
            reach_cache[ireach_cache] = bidx
            ireach_cache += 1

        usreach_cache[iusreach_cache] = -usreachlen
        iusreach_cache += 1
        if usreachlen > 0:
            for bidx in mc_reach.binary_find(parameter_idx, connections[reach[0]]):
                usreach_cache[iusreach_cache] = bidx
                iusreach_cache += 1
        # print(np.asarray(connections[reach[0]]))
        print(np.asarray(usreach_cache))

    maxreachlen = max(reach_sizes)
    buf = np.empty((maxreachlen, buf_cols), dtype='float32')
    out_buf = np.empty((maxreachlen, 3), dtype='float32')

    drows_tmp = np.arange(maxreachlen, dtype=np.intp)
    #cdef Py_ssize_t[:] drows
    #cdef float qup, quc
    timestep = 0
    #cdef int ts_offset

    ireach_cache = 0
    print(f"reach_cache.shape[0] {reach_cache.shape[0]}")
    print(f"cache_sizes {np.asarray(reach_group_cache_sizes)} reach_groups {np.asarray(reach_groups)}")

    print(f"reach_sizes {reach_sizes}, usreach_sizes {usreach_sizes}")
    print(f"connections {[seg for reach in reaches for seg in reach]}")
    print(f"reaches {reaches}")
    print(f"reach_cache {np.asarray(reach_cache)}")
    print(f"usreach_cache {np.asarray(usreach_cache)}")
    #with nogil:
    if 1 == 1:
        while timestep < nsteps:
            ts_offset = timestep * 3

            ireach_cache = 0
            iusreach_cache = 0
            ireach = 0
            for group_i in range(len(reach_group_cache_sizes)):
                #while ireach_cache < reach_cache.shape[0]:
                ireach_cache_end = ireach_cache + reach_group_cache_sizes[group_i] + reach_groups[group_i]
                while ireach_cache < ireach_cache_end:
                    
                    reachlen = -reach_cache[ireach_cache]
                    usreachlen = -usreach_cache[iusreach_cache]

                    ireach_cache += 1
                    iusreach_cache += 1
                    #print(ireach_cache, iusreach_cache, np.asarray(reach_cache, dtype=np.intp), np.asarray(usreach_cache, dtype=np.intp))

                    qup = 0.0
                    quc = 0.0
                    for i in range(usreachlen):
                        quc += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset]
                        if timestep > 0:
                            qup += flowveldepth[usreach_cache[iusreach_cache + i], ts_offset - 3]

                    buf_view = buf[:reachlen, :]
                    out_view = out_buf[:reachlen, :]
                    drows = drows_tmp[:reachlen]
                    srows = reach_cache[ireach_cache:ireach_cache+reachlen]

                    utils.fill_buffer_column(srows, timestep, drows, 0, qlat_values, buf_view)
                    for i in range(scols.shape[0]):
                        utils.fill_buffer_column(srows, scols[i], drows, i + 1, parameter_values, buf_view)
                        # fill buffer with qdp, depthp, velp
                    if timestep > 0:
                        utils.fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view)
                        utils.fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view)
                        utils.fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view)
                    else:
                        # fill buffer with constant
                        for i in range(drows.shape[0]):
                            buf_view[drows[i], 10] = 0.0
                            buf_view[drows[i], 11] = 0.0
                            buf_view[drows[i], 12] = 0.0

                    if assume_short_ts:
                        quc = qup

                    if timestep < 0:
                        print(f"ts {timestep}, current reach_cache {reach_cache[ireach_cache]}, qup {qup}, quc {quc}, reachlen {reachlen}, buf_view {np.asarray(buf_view)}, out_view {np.asarray(out_view)}")
                    mc_reach.compute_reach_kernel(qup, quc, reachlen, buf_view, out_view)
                    if timestep == 3:
                        print(f"ts {timestep}, reach {reaches[ireach]} segment indexes (reach_cache) {np.asarray(srows)}, qup {qup}, quc {quc}, reachlen {reachlen} ", end="")
                        print(f"upstream segments {[[seg for reach in reaches for seg in reach][r] for r in usreach_cache[iusreach_cache:iusreach_cache + usreachlen]]}")
                        #print(f"{np.asarray(out_view)}")

                    # copy out_buf results back to flowdepthvel
                    for i in range(3):
                        utils.fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)

                    # Update indexes to point to next reach
                    ireach += 1
                    ireach_cache += reachlen
                    iusreach_cache += usreachlen
                    
            timestep += 1

    return np.asarray(parameter_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')
