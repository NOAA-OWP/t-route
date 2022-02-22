import cython
import numpy as np
from array import array
from numpy cimport ndarray
cimport numpy as np
#from libc.stdio cimport printf

from .fortran_wrappers cimport c_muskingcungenwm, c_reachcompute

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
        QVD *rv):

    cdef:
        float qdc = 0.0
        float depthc = 0.0
        float velc = 0.0
        float ck = 0.0
        float cn = 0.0
        float X = 0.0

    #printf("reach.pyx before %3.9f\t", depthc)
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
        &depthc,
        &ck,
        &cn,
        &X)
    #printf("reach.pyx after %3.9f\t", depthc)

    rv.qdc = qdc
    rv.depthc = depthc
    rv.velc = velc

    # to do: make these additional variable's conditional, somehow
    rv.ck = ck
    rv.cn = cn
    rv.X = X
    
@cython.boundscheck(False)
cpdef void computereach(float dt,
        int nseg,
        int nts,
        float[:] qup,
        float[:] quc,
        float[:] qdp,
        float[::1,:] ql,
        float[:] dx,
        float[:] bw,
        float[:] tw,
        float[:] twcc,
        float[:] n,
        float[:] ncc,
        float[:] cs,
        float[:] s0,
        float[:] velp,
        float[:] depthp,
        float[::1,:] qdc,
        float[::1,:] velc,
        float[::1,:] depthc,
        float[::1,:] cn,
        float[::1,:] ck,
        float[::1,:] X):

    c_reachcompute(
        &dt,
        &nseg,
        &nts,
        &qup[0],
        &quc[0],
        &qdp[0],
        &ql[0,0],
        &dx[0],
        &bw[0],
        &tw[0],
        &twcc[0],
        &n[0],
        &ncc[0],
        &cs[0],
        &s0[0],
        &velp[0],
        &depthp[0],
        &qdc[0,0],
        &velc[0,0],
        &depthc[0,0],
        &ck[0,0],
        &cn[0,0],
        &X[0,0])

