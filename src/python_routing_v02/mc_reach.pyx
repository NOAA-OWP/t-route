# cython: language_level=3, boundscheck=False, wraparound=False

import numpy as np
from itertools import chain
from operator import itemgetter
from numpy cimport ndarray
cimport numpy as cnp


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


cpdef void muskingcunge(float dt,
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
        float[:] out) nogil:
    cdef float qdc, velc, depthc

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
    out[0] = qdc
    out[1] = depthc
    out[2] = velc
    #return qdc, depthc, velc


cpdef void compute_reach(Py_ssize_t nreach, const long[:,:] reach, float qup, float quc, float[:, :] flowdepthvel, const float[:, :] data_values, bint assume_short_ts=False):
    cdef:
        long i, current_segment
        Py_ssize_t ix
        float dt, bw, tw, twcc, dx, n_manning, n_manning_cc
        float cs, s0
        float qlat, qdp, depthp, velp
        float qdc, depthc, velc

    cdef float out[3]
    cdef float[:] out_view = out

    with nogil:
        for ix in range(nreach):
            i = reach[ix, 0]
            current_segment = reach[ix, 1]

            dt = 60.0
            bw = data_values[current_segment, 0]
            tw = data_values[current_segment, 1]
            twcc = data_values[current_segment, 2]
            dx = data_values[current_segment, 3]
            n_manning = data_values[current_segment, 4]
            n_manning_cc = data_values[current_segment, 5]
            cs = data_values[current_segment, 6]
            s0 = data_values[current_segment, 7]

            flowdepthvel[i, 7] = 10.0

            qdp = flowdepthvel[i, 0]
            depthp = flowdepthvel[i, 1]
            velp = flowdepthvel[i, 2]
            qlat = flowdepthvel[i, 7]

            flowdepthvel[i, :4] = flowdepthvel[i, 4:]

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
                out_view
            )


            flowdepthvel[i, 4] = qdc = out[0]
            flowdepthvel[i, 5] = out[1]
            flowdepthvel[i, 6] = out[2]
            #flowdepthvel[i, 4:7] = qdc, depthc, velc
            quc = qdc
            qup = qdp


def compute_network(int nsteps, object reaches, object connections, const float[:,:] data_values, bint assume_short_ts=False):

    cdef object local_idx, global_idx, us
    cdef tuple reach

    cdef ndarray[long, ndim=2] findex = np.array(sorted(chain.from_iterable(reaches), key=itemgetter(0)))

    cdef ndarray[float, ndim=2] flowdepthvel = np.zeros((len(findex), 8), dtype='float32')

    cdef long[:] idx_view = findex[:,0]

    cdef list global_idxs = []
    cdef list local_idxs = []
    cdef list upstream_segs = []
    cdef list x, y, us_segs
    cdef object searched, tmp
    cdef dict dtmp
    for reach in reaches:
        x = []
        y = []
        for tmp in reach:
            x.append(tmp[0])
            y.append(tmp[1])
        global_idxs.append(x)
        searched = np.searchsorted(idx_view, x).tolist()
        local_idxs.append(searched)

        # unpack segments
        us_segs = []
        if y[0] in connections:
            dtmp = connections[y[0]]
            for tmp in dtmp['children']:
                us_segs.append(tmp[0])
        searched = np.searchsorted(idx_view, us_segs).tolist()
        upstream_segs.append(searched)


    cdef tuple ri
    cdef long xi
    cdef ndarray[long, ndim=2] reaches_index
    cdef int ts
    cdef float qup, quc
    for ts in range(nsteps):
        # breakpoint()
        for ri in zip(local_idxs, global_idxs, upstream_segs):
            local_idx = ri[0]
            global_idx = ri[1]
            us = ri[2]

            qup = 0.0
            quc = 0.0
            for xi in us:
                qup += flowdepthvel[xi, 0]
                quc += flowdepthvel[xi, 4]
            reach_indexes = np.array(tuple(zip(local_idx, global_idx)), dtype=long)
            compute_reach(reach_indexes.shape[0], reach_indexes, qup, quc, flowdepthvel, data_values,
                                   assume_short_ts=assume_short_ts)
    return idx_view, flowdepthvel


