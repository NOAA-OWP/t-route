cpdef Py_ssize_t bisect_left_long(const long[:] arr, long el) nogil
cpdef Py_ssize_t bisect_left(object arr, object el)
cpdef list binary_find(object arr, object els)
cpdef Py_ssize_t[:,:] build_upstream_graph(object rconnections, const long[:] data_index)