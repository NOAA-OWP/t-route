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
        float *out) nogil:
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


cdef void compute_reach(Py_ssize_t nreach, const long[:,:] reach, float qup, float quc, float[:, :] flowdepthvel, const float[:, :] data_values, bint assume_short_ts=False) nogil:
    cdef:
        long i, current_segment
        Py_ssize_t ix
        float dt, bw, tw, twcc, dx, n_manning, n_manning_cc
        float cs, s0
        float qlat, qdp, depthp, velp
        float qdc

    cdef float out[3]
    #cdef float[:] out_view = out


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
            out
        )


        flowdepthvel[i, 4] = qdc = out[0]
        flowdepthvel[i, 5] = out[1]
        flowdepthvel[i, 6] = out[2]
        #flowdepthvel[i, 4:7] = qdc, depthc, velc
        quc = qdc
        qup = qdp


def compute_network(int nsteps, list reaches, dict connections, const float[:,:] data_values, bint assume_short_ts=False):

    cdef long[:] us
    cdef long[:, :] reach

    cdef ndarray[long, ndim=2] findex = np.vstack(reaches)
    findex.sort(axis=0)
    cdef long[:] idx_view = findex[:,0]

    cdef float[:, :] flowdepthvel = np.zeros((len(idx_view), 8), dtype='float32')

    cdef list idxs = []
    cdef list upstream_segs = []
    cdef long[:] li_tmp, gi_tmp
    cdef long[:,:] us_segs
    cdef dict dtmp
    for reach in reaches:
        gi_tmp = reach[:,0]
        li_tmp = np.searchsorted(idx_view, gi_tmp)
        idxs.append(np.stack([li_tmp, gi_tmp], axis=1))

        # unpack segments
        if reach[0,1] in connections:
            dtmp = connections[reach[0,1]]
            if dtmp['children'] is not None:
                us_segs = dtmp['children']
                upstream_segs.append(np.searchsorted(idx_view, us_segs[:,0]))
            else:
                upstream_segs.append(np.array((), dtype=long))
        else:
            upstream_segs.append(np.array((), dtype=long))

    cdef long[:,:] reach_indexes
    cdef int ts, ix, ri, us_size
    cdef float qup, quc

    ts = 0
    while ts < nsteps:
        # breakpoint()
        for ri in range(len(idxs)):
            reach_indexes = idxs[ri]
            us = upstream_segs[ri]
            us_size = us.size

            with nogil:
                qup = 0.0
                quc = 0.0
                for ix in range(us_size):
                    qup += flowdepthvel[us[ix], 0]
                    quc += flowdepthvel[us[ix], 4]

                compute_reach(reach_indexes.shape[0], reach_indexes, qup, quc, flowdepthvel, data_values,
                                       assume_short_ts=assume_short_ts)
        ts += 1
    return idx_view, flowdepthvel


