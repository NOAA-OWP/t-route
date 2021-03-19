from numpy import empty
from numpy cimport ndarray as ar
cimport numpy as np

cdef extern from "pyuniflowtzlt.h" nogil:
    void c_uniflowtzlt(int* mxncomp_g, 
                       int* nrch_g, 
                       double* bo_ar_g, 
                       double* traps_ar_g, 
                       double* tw_ar_g, 
                       double* twcc_ar_g, 
                       double* mann_ar_g, 
                       double* manncc_ar_g, 
                       double* so_ar_g, 
                       int* nhincr_m_g, 
                       int* nhincr_f_g, 
                       double* ufhlt_m_g, 
                       double* ufqlt_m_g, 
                       double* ufhlt_f_g, 
                       double* ufqlt_f_g, 
                       int* frnw_col, 
                       double* dfrnw_g, 
                       double* timesdepth_g)

def uniflow_lookuptable(int mxncomp_g, 
                        int nrch_g, 
                        bo_ar_g, 
                        traps_ar_g, 
                        tw_ar_g, 
                        twcc_ar_g, 
                        mann_ar_g, 
                        manncc_ar_g, 
                        so_ar_g, 
                        int nhincr_m_g, 
                        int nhincr_f_g,
                        int frnw_col, 
                        dfrnw_g, 
                        double timesdepth_g): 
    cdef:
        ar[double,ndim=2] bo = bo_ar_g
        ar[double,ndim=2] traps = traps_ar_g
        ar[double,ndim=2] tw = tw_ar_g
        ar[double,ndim=2] twcc = twcc_ar_g
        ar[double,ndim=2] mann = mann_ar_g
        ar[double,ndim=2] manncc = manncc_ar_g
        ar[double,ndim=2] so = so_ar_g
        ar[double,ndim=2] dfrnw= dfrnw_g  
        ar[double,ndim=3] ufhlt_m_g = empty((mxncomp_g, nrch_g, nhincr_m_g), order='F')
        ar[double,ndim=3] ufqlt_m_g = empty((mxncomp_g, nrch_g, nhincr_m_g), order='F')
        ar[double,ndim=3] ufhlt_f_g = empty((mxncomp_g, nrch_g, nhincr_f_g), order='F')
        ar[double,ndim=3] ufqlt_f_g = empty((mxncomp_g, nrch_g, nhincr_f_g), order='F')  
 
    with nogil:
        c_uniflowtzlt(&mxncomp_g, 
                      &nrch_g, 
                      <double*> bo.data, 
                      <double*> traps.data, 
                      <double*> tw.data,
                      <double*> twcc.data,
                      <double*> mann.data,
                      <double*> manncc.data,
                      <double*> so.data,
                      &nhincr_m_g, 
                      &nhincr_f_g,
                      <double*> ufhlt_m_g.data,
                      <double*> ufqlt_m_g.data,
                      <double*> ufhlt_f_g.data,
                      <double*> ufqlt_f_g.data,
                      &frnw_col,
                      <double*> dfrnw.data,
                      &timesdepth_g)
    return ufhlt_m_g, ufqlt_m_g, ufhlt_f_g, ufqlt_f_g

