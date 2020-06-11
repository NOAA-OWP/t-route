# cython: language_level=3, boundscheck=False, wraparound=False, profile=True

import numpy as np
from itertools import chain
from operator import itemgetter
from numpy cimport ndarray
cimport numpy as np


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
            idxs.append(-1)
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

cdef void compute_reach(int nreach, const long[:,:] reach, float qup, float quc, float[:, :] flowdepthvel, const float[:, :] data_values, bint assume_short_ts=False) nogil:
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
        if i == -1:
            continue

        dt = 60.0
        bw = data_values[current_segment, 0]
        tw = data_values[current_segment, 1]
        twcc = data_values[current_segment, 2]
        dx = data_values[current_segment, 3]
        n_manning = data_values[current_segment, 4]
        n_manning_cc = data_values[current_segment, 5]
        cs = data_values[current_segment, 6]
        s0 = data_values[current_segment, 7]
        qlat = data_values[current_segment, 8]

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


def compute_network(int nsteps, list reaches, dict connections, const float[:,:] data_values, bint assume_short_ts=False):

    cdef long[:] findex = np.vstack(reaches)[:,0]
    cdef long[:] idx_view = np.sort(findex)

    cdef float[:,:] flowdepthvel = np.zeros((len(idx_view), 6), dtype='float32')

    cdef list idxs = []
    cdef list upstream_segs = []
    cdef long[:] gi_tmp, li_tmp
    cdef long[:,:] us_segs
    cdef long[:,:] reach
    cdef dict dtmp
    for reach in reaches:
        gi_tmp = reach[:,0]
        li_tmp = np.array(binary_find(idx_view, gi_tmp), dtype=long)
        idxs.append(np.stack([li_tmp, gi_tmp], axis=1))

        # unpack segments
        if reach[0,1] in connections:
            dtmp = connections[reach[0,1]]
            if dtmp['children'] is not None:
                us_segs = dtmp['children']
                upstream_segs.append(np.array(binary_find(idx_view, us_segs[:,0]), dtype=long))
            else:
                upstream_segs.append(np.array((), dtype=long))
        else:
            upstream_segs.append(np.array((), dtype=long))

    #print(idxs, upstream_segs)
    cdef long[:,:] reach_indexes
    cdef int ts, ix, ri, us_size
    cdef float qup, quc
    cdef long[:] us
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
                    quc += flowdepthvel[us[ix], 3]

                #print(f"ts={ts}\tqup={qup}\tquc={quc}\nRI={np.asarray(reach_indexes, dtype=long)}")
                    compute_reach(reach_indexes.shape[0], reach_indexes, qup, quc, flowdepthvel, data_values,
                                           assume_short_ts=assume_short_ts)

        # Advance timestep
        #with np.printoptions(precision=6, suppress=True, linewidth=180, edgeitems=5):
        #    print(f"FDV=FDV={np.asarray(flowdepthvel, dtype='float32')}")
        flowdepthvel[:, :3] = flowdepthvel[:, 3:]
        ts += 1

    # coerce back to numpy
    return np.asarray(idx_view[:], dtype=np.long), np.asarray(flowdepthvel[:], dtype='float32')

