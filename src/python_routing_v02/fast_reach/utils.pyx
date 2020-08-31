# cython: language_level=3, boundscheck=True, wraparound=False

import numpy as np
cimport cython
cimport numpy as np

@cython.boundscheck(False)
cpdef Py_ssize_t bisect_left_long(const long[:] arr, long el) nogil:
    cdef Py_ssize_t L = 0
    cdef Py_ssize_t R = len(arr) - 1
    cdef Py_ssize_t m
    cdef long cand

    while L <= R:
        m = (L + R) // 2
        cand = arr[m]
        if cand < el:
            L = m + 1
        elif cand > el:
            R = m - 1
        else:
            return m
    return -1

@cython.boundscheck(False)
cpdef Py_ssize_t bisect_left(object arr, object el):
    cdef Py_ssize_t L = 0
    cdef Py_ssize_t R = len(arr) - 1
    cdef Py_ssize_t m
    cdef object cand

    while L <= R:
        m = (L + R) // 2
        cand = arr[m]
        if cand < el:
            L = m + 1
        elif cand > el:
            R = m - 1
        else:
            return m
    return -1


@cython.boundscheck(False)
cpdef list binary_find(object arr, object els):
    """
    Find elements in els in arr.
    Args:
        arr: Array to search. Must be sorted
        els:

    Returns:

    """
    cdef list idxs = []
    cdef Py_ssize_t r

    for el in els:
        r = bisect_left(arr, el)
        if r > 0:
            idxs.append(r)
        else:
            raise ValueError(f"element {el} not found in {np.asarray(arr)}")
    return idxs


@cython.boundscheck(False)
cpdef Py_ssize_t[:,:] build_upstream_graph(object rconnections, const long[:] data_index):
    """
    Build an array representation of upstream connections.

    Returns a 2D array of shape 2xn. The first row are the destination nodes in sorted order
    The second row are the index of the source nodes in the data_index.
    Data_index is assumed to be sorted, and a binary search is used to find occurrence of src.
    This function assumes that node ids are integers.

    Example:
        2 2 2 3 3 5 -> node ids
        0 0 0 1 2 3 -> index of upstream nodes in data_index!
    """
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t ind, k, len_v
    cdef object v
    cdef object kv
    #cdef const long[:] v_view
    cdef Py_ssize_t[:, :] arr_view

    # Get the cumulative size of all values
    cdef long total = 0
    for _v in rconnections.values():
        total += len(_v)

    arr_view = np.empty((2, total), dtype=np.intp)
    for kv in sorted(rconnections.items()):
        k = kv[0]
        v = kv[1]
        len_v = len(v)
        for ind in range(i, i+len_v):
            arr_view[0, ind] = k
            arr_view[1, ind] = bisect_left_long(data_index, v[ind-i])
        i += len_v
    return arr_view
