from numpy import linspace, empty
from numpy cimport ndarray as ar
cimport numpy as np

cdef extern from "pygfunc.h" nogil:
    void c_gfunc(double* x, int* n, int* m, double* a, double* b, double* c) 

def f2cytest(double x, double a, double b, int n): 
    cdef:
        ar[double] ax = linspace(a, b, n)
        ar[double,ndim=2] c = empty((n, n), order='F')
    with nogil:
        c_gfunc(&x, &n, &n, <double*> ax.data, <double*> ax.data, <double*> c.data)
    return c


