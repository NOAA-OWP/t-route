# cython: language_level=3, boundscheck=True, wraparound=False, profile=True

import numpy as np
from itertools import chain
from operator import itemgetter
from numpy cimport ndarray
cimport numpy as np
cimport cython


# TODO: Make this QVD
cdef struct QVD:
    float qdc
    float velc
    float depthc

cdef struct QH:
    float resoutflow
    float reslevel

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


cdef extern from "pyResLevelPool.h":
    void c_levelpool_physics(float *dt,
                              float *qi0,
                              float *qi1,
                              float *ql,
                              float *ar,
                              float *we,
                              float *maxh,
                              float *wc,
                              float *wl,
                              float *dl,
                              float *oe,
                              float *oc,
                              float *oa,
                              float *H0,
                              float *H1,
                              float *qo1) nogil;

cdef extern from "pyMCsingleSegStime_NoLoop.h":
    void c_muskingcungenwm(float *dt,
                                  float *qup,
                                  float *quc,
                                  float *qdp,
                                  float *ql,
                                  float *dx,
                                  float *bw,
                                  float *tw,
                                  float *twcc,
                                  float *n,
                                  float *ncc,
                                  float *cs,
                                  float *s0,
                                  float *velp,
                                  float *depthp,
                                  float *qdc,
                                  float *velc,
                                  float *depthc) nogil;


@cython.boundscheck(False)
cdef void levelpool_physics(float dt,
        float qi0,
        float qi1,
        float ql,
        float ar,
        float we,
        float maxh,
        float wc,
        float wl,
        float dl,
        float oe,
        float oc,
        float oa,
        float H0,
        QH *rv) nogil:
    cdef float H1 = 0.0
    cdef float qo1 = 0.0

    c_levelpool_physics(
        &dt,
        &qi0,
        &qi1,
        &ql,
        &ar,
        &we,
        &maxh,
        &wc,
        &wl,
        &dl,
        &oe,
        &oc,
        &oa,
        &H0,
        &H1,
        &qo1)
    rv.reslevel = H1
    rv.resoutflow = qo1


cpdef compute_reservoir_kernel(float dt,
        float qi0,
        float qi1,
        float ql,
        float ar,
        float we,
        float maxh,
        float wc,
        float wl,
        float dl,
        float oe,
        float oc,
        float oa,
        float H0):

    cdef QH rv
    cdef QH *out = &rv

    levelpool_physics(dt,
        qi0,
        qi1,
        ql,
        ar,
        we,
        maxh,
        wc,
        wl,
        dl,
        oe,
        oc,
        oa,
        H0,
        out)

    return rv
    

@cython.boundscheck(False)
cdef void muskingcunge(float dt,
        float qup,
        float quc,
        float qdp,
        float ql,
        float dx,
        float bw,
        float tw,
        float twcc,
        float n,
        float ncc,
        float cs,
        float s0,
        float velp,
        float depthp,
        QVD *rv) nogil:
    cdef float qdc, depthc, velc
    qdc = 0.0
    depthc = 0.0
    velc = 0.0

    c_muskingcungenwm(
        &dt,
        &qup,
        &quc,
        &qdp,
        &ql,
        &dx,
        &bw,
        &tw,
        &twcc,
        &n,
        &ncc,
        &cs,
        &s0,
        &velp,
        &depthp,
        &qdc,
        &velc,
        &depthc)
    rv.qdc = qdc
    rv.depthc = depthc
    rv.velc = velc


@cython.boundscheck(False)
cpdef float[:, :] compute_reach_kernel(float qup, float quc, int nreach, const float[:,:] input_buf, float[:, :] output_buf) nogil:
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

    For computing schemes that require a function return value,
    a reference to output_buf is returned from this function.
    """
    cdef QVD rv
    cdef QVD *out = &rv

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

        output_buf[i, 0] = quc = out.qdc
        output_buf[i, 1] = out.velc
        output_buf[i, 2] = out.depthc

        qup = qdp
    return output_buf

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
    const float[:, :] qlat_values, 
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
        assume_short_ts (bool): Assume short time steps (quc = qup)

    Notes:
        Array dimensions are checked as a precondition to this method.
    """
    # Check shapes
    if qlat_values.shape[0] != data_idx.shape[0] or qlat_values.shape[1] != nsteps:
        raise ValueError(f"Qlat shape is incorrect: expected ({data_idx.shape[0], nsteps}), got ({qlat_values.shape[0], qlat_values.shape[1]})")
    if data_values.shape[0] != data_idx.shape[0] or data_values.shape[1] != data_cols.shape[0]:
        raise ValueError(f"data_values shape mismatch")

    cdef float[:,::1] flowveldepth = np.zeros((data_idx.shape[0], nsteps * 3), dtype='float32')

    cdef Py_ssize_t[:] srows, drows_tmp, usrows
    cdef float[:, ::1] buf, buf_view
    cdef float[:, ::1] out_buf, out_view

    cdef Py_ssize_t[:] scols = np.array(column_mapper(data_cols), dtype=np.intp)
    cdef int buf_cols = 13
    cdef Py_ssize_t i, ireach, ireach_cache, iusreach_cache

    cdef list reach_sizes = list(map(len, reaches))
    cdef list usreach_sizes = [len(connections.get(reach[0], ())) for reach in reaches]

    cdef list reach, bf_results
    cdef int reachlen, usreachlen
    cdef Py_ssize_t bidx
    cdef list buf_cache = []
    cdef Py_ssize_t[:] reach_cache = np.empty(sum(reach_sizes) + len(reach_sizes), dtype=np.intp)
    cdef Py_ssize_t[:] usreach_cache = np.empty(sum(usreach_sizes) + len(usreach_sizes), dtype=np.intp)
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

                fill_buffer_column(srows, timestep, drows, 0, qlat_values, buf_view)
                for i in range(scols.shape[0]):
                        fill_buffer_column(srows, scols[i], drows, i + 1, data_values, buf_view)
                    # fill buffer with qdp, depthp, velp
                if timestep > 0:
                    fill_buffer_column(srows, ts_offset - 3, drows, 10, flowveldepth, buf_view)
                    fill_buffer_column(srows, ts_offset - 2, drows, 11, flowveldepth, buf_view)
                    fill_buffer_column(srows, ts_offset - 1, drows, 12, flowveldepth, buf_view)
                else:
                    # fill buffer with constant
                    for i in range(drows.shape[0]):
                        buf_view[drows[i], 10] = 0.0
                        buf_view[drows[i], 11] = 0.0
                        buf_view[drows[i], 12] = 0.0

                if assume_short_ts:
                    quc = qup

                # TODO: Implement Reach-based ordering
                # for order in orders:
                    # for reach in ordered_reaches[order]:
                        # jobs.append(delayed(
                # TODO: Implement reach-based if statement to segregate waterbody computations 
                # TODO: from regular reaches
                compute_reach_kernel(qup, quc, reachlen, buf_view, out_view)

                # copy out_buf results back to flowdepthvel
                for i in range(3):
                    fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowveldepth)

                # Update indexes to point to next reach
                ireach_cache += reachlen
                iusreach_cache += usreachlen
                
            timestep += 1
    return np.asarray(data_idx, dtype=np.intp), np.asarray(flowveldepth, dtype='float32')
