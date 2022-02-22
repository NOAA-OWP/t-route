cdef struct QVD:
    float qdc
    float velc
    float depthc
    float cn
    float ck
    float X


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
        QVD *rv)

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
        float[::1,:] X)

