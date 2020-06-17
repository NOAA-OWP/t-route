# cython: language_level=3, boundscheck=True, wraparound=False, profile=True

import numpy as np
from itertools import chain
from operator import itemgetter
from numpy cimport ndarray
cimport numpy as np
cimport cython


ctypedef fused float:
    float
    double

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
            break
    return idxs



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

cpdef list pymuskingcunge(float dt,
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
        float depthp):
    cdef float qdc, velc, depthc
    qdc = 0.0
    velc = 0.0
    depthc = 0.0

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
    return [qdc, depthc, velc]
    #out[0] = qdc
    #out[1] = depthc
    #out[2] = velc

cdef struct FDV:
    float qdc
    float depthc
    float velc


@cython.boundscheck(False)
@cython.wraparound(False)
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
        FDV *rv) nogil:
    cdef float qdc, depthc, velc
    qdc = 0.0
    depthc = 0.0
    velc = 0.0

    #print(f"muskingcunge({dt},{qup},{quc},{qdp},{ql},{dx},{bw},{tw},{twcc},{n},{ncc},{cs},{s0},{velp},{depthp})")
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
    #print(f"rv=[qdc={qdc}, depthc={depthc}, velc={velc}]")
    rv.qdc = qdc
    rv.depthc = depthc
    rv.velc = velc


@cython.boundscheck(False)
@cython.wraparound(False)
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
    cdef FDV rv
    cdef FDV *out = &rv

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
        depthp = input_buf[i, 11]
        velp = input_buf[i, 12]

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
        output_buf[i, 1] = out.depthc
        output_buf[i, 2] = out.velc

        qup = qdp


cdef void compute_reach(int nreach, const long[:,:] reach, float qup, float quc, float[:, :] flowdepthvel, const float[:, :] data_values, const float[:] qlats, bint assume_short_ts=False):
    cdef:
        long i, current_segment
        int ix
        float dt, bw, tw, twcc, dx, n_manning, n_manning_cc
        float cs, s0
        float qlat, qdp, depthp, velp
        float qdc

    cdef FDV rv
    cdef FDV *out = &rv

    for ix in range(nreach):
        # i: local index for the segment (used for flowdepthvel)
        # current_segment: global index for the segment (used for looking updata_values)
        i = reach[ix, 0]
        current_segment = reach[ix, 1]
        if current_segment == -1:
            print("invalid segment index encoutered not at end of reach!")
            break

        dt = 60.0
        bw = data_values[current_segment, 0]
        tw = data_values[current_segment, 1]
        twcc = data_values[current_segment, 2]
        dx = data_values[current_segment, 3]
        n_manning = data_values[current_segment, 4]
        n_manning_cc = data_values[current_segment, 5]
        cs = data_values[current_segment, 6]
        s0 = data_values[current_segment, 7]
        qlat = qlats[current_segment]

        qdp = flowdepthvel[i, 0]
        depthp = flowdepthvel[i, 1]
        velp = flowdepthvel[i, 2]

        if assume_short_ts:
            quc = qup

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
            n_manning,
            n_manning_cc,
            cs,
            s0,
            velp,
            depthp,
            out
        )

        qdc = out.qdc
        flowdepthvel[i, 3] = qdc
        flowdepthvel[i, 4] = out.depthc  # depthc
        flowdepthvel[i, 5] = out.velc # velc
        #print(f"flowdepthvel[i] = {list(flowdepthvel[i, 3:])}")
        quc = qdc
        qup = qdp


cdef void fill_buffer_column(const Py_ssize_t[:] srows,
    const Py_ssize_t scol,
    const Py_ssize_t[:] drows,
    const Py_ssize_t dcol,
    const float[:, :] src, float[:, ::1] out) nogil:

    cdef Py_ssize_t i
    for i in range(srows.shape[0]):
        out[drows[i], dcol] = src[srows[i], scol]


def column_mapper(object src_cols) -> object:
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

def compute_network(int nsteps, object reaches, object connections, 
    const long[:] data_idx, object[:] data_cols, const float[:,:] data_values, 
    const float[:, :] qlat_values, 
    bint assume_short_ts=False) -> object:
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

    cdef float[:,::1] flowdepthvel = np.zeros((data_idx.shape[0], nsteps * 3), dtype='float32')
    cdef int timestep = 0
    cdef int ts_offset

    cdef Py_ssize_t[:] srows, drows_tmp
    cdef float[:, ::1] buf, buf_view
    cdef float[:, ::1] out_buf, out_view
    cdef object reach  # list
    cdef int reach_length

    cdef Py_ssize_t[:] scols = np.array(column_mapper(data_cols), dtype=np.intp)
    cdef int buf_cols = 13
    cdef Py_ssize_t i, ireach
    cdef float qup, quc

    cdef int maxreachlen = 0
    cdef object buf_cache = []
    for reach in reaches:
        maxreachlen = max(maxreachlen, len(reach))
        srows = np.array(binary_find(data_idx, reach), dtype=np.intp)
        buf_cache.append(srows)

    buf = np.empty((maxreachlen, buf_cols), dtype='float32')
    out_buf = np.empty((maxreachlen, 3), dtype='float32')

    drows_tmp = np.arange(maxreachlen, dtype=np.intp)
    cdef Py_ssize_t[:] drows
    while timestep < nsteps:
        ts_offset = timestep * 3

        #print(f"Timestep: {timestep}".ljust(120, '='))
        for ireach in range(len(reaches)):
            reach = reaches[ireach]
            srows = buf_cache[ireach]
            
            # compute qup/quc
            qup = 0.0
            quc = 0.0
            if reach[0] in connections:
                for i in binary_find(data_idx, connections[reach[0]]):
                    quc += flowdepthvel[i, ts_offset]
                    if timestep > 0:
                        qup += flowdepthvel[i, ts_offset - 3]

            with nogil:
                reach_length = srows.shape[0]
                buf_view = buf[:reach_length, :]
                out_view = out_buf[:reach_length, :]
                drows = drows_tmp[:reach_length]

                # fill the buffer with qlat
                fill_buffer_column(srows, timestep, drows, 0, qlat_values, buf_view)
                for i in range(scols.shape[0]):
                    fill_buffer_column(srows, scols[i], drows, i + 1, data_values, buf_view)
                # fill buffer with qdp, depthp, velp
                if timestep > 0:
                    fill_buffer_column(srows, ts_offset - 3, drows, 10, flowdepthvel, buf_view)
                    fill_buffer_column(srows, ts_offset - 2, drows, 11, flowdepthvel, buf_view)
                    fill_buffer_column(srows, ts_offset - 1, drows, 12, flowdepthvel, buf_view)
                else:
                    # fill buffer with constant
                    for i in range(drows.shape[0]):
                        buf_view[drows[i], 10] = 0.0
                        buf_view[drows[i], 11] = 0.0
                        buf_view[drows[i], 12] = 0.0

                if assume_short_ts:
                    quc = qup

                compute_reach_kernel(qup, quc, reach_length, buf_view, out_view)

                # copy out_buf results back to flowdepthvel
                for i in range(3):
                    fill_buffer_column(drows, i, srows, ts_offset + i, out_view, flowdepthvel)
            
        timestep += 1
    return np.asarray(data_idx, dtype=np.intp), np.asarray(flowdepthvel, dtype='float32')
