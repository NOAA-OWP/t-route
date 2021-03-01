from numpy cimport ndarray
from numpy import empty

cdef extern nogil:
    void c_meshexp(double* r_min , double* r_max, double* a, int* N, double* mesh)

def mesh_exp(double r_min, double r_max, double a, int N):
    cdef ndarray[double, mode="c"] mesh = empty(N+1, dtype="double")
    with nogil:
        c_meshexp(&r_min, &r_max, &a, &N, &mesh[0])
    return mesh
