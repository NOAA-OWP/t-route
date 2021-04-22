from reach cimport QVD_double

cdef extern from "muskingumcunge.h" nogil:
    void muskingum_cunge(
        double dt,
        double qup,
        double quc,
        double qdp,
        double ql,
        double dx,
        double bw,
        double tw,
        double twcc,
        double n,
        double ncc,
        double cs,
        double s0,
        double velp,
        double depthp,
        QVD_double *rv
    );