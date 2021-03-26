cdef extern from "pydiffusive.h":
    void c_diffnw(int *mxncomp_g,
                int *nrch_g,
                double *z_ar_g,
                int *ntss_ev_g,
                double *q_ev_g,
                double *elv_ev_g) nogil;
