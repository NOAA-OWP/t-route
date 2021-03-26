import cython
import numpy as np
cimport numpy as np

from fortran_wrappers cimport c_diffnw


# TO DO load some example inputs to test the module

@cython.boundscheck(False)
cdef void diffnw(
                int mxncomp_g,
                int nrch_g,
                double[::1,:] z_ar_g,
                int ntss_ev_g,
                double[:,:,:] out_q,
                double[:,:,:] out_elv
            ):
    
    cdef:
        double[::1,:,:] q_ev_g = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double, order = 'F') 
        double[::1,:,:] elv_ev_g = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double, order = 'F') 
        
    c_diffnw(
            &mxncomp_g,
            &nrch_g,
            &z_ar_g[0,0],
            &ntss_ev_g,
            &q_ev_g[0,0,0],
            &elv_ev_g[0,0,0]
        )
    
    # copy data from Fortran to Python memory view
    out_q[:,:,:] = q_ev_g[::1,:,:]
    out_elv[:,:,:] = elv_ev_g[::1,:,:]
    
cpdef object compute_diffusive(
                 int mxncomp_g,
                 int nrch_g,
                 double[::1,:] z_ar_g,
                 int ntss_ev_g
                ):

    cdef:
        double[:,:,:] out_q = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)
        double[:,:,:] out_elv = np.empty([ntss_ev_g,mxncomp_g,nrch_g], dtype = np.double)
    
    diffnw(
         mxncomp_g,
         nrch_g,
         z_ar_g,
         ntss_ev_g,
         out_q,
         out_elv
        )
    
    return np.asarray(out_q), np.asarray(out_elv)