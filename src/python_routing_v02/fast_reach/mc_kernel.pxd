cdef struct QVD:
    float qdc
    float velc
    float depthc
    float cn
    float ck
    float X


cdef void cython_muskingcunge(
    float dt,
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
    QVD *rv,
) nogil
