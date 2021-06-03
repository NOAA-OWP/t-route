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
    cdef Py_ssize_t ind, k, len_v, kl
    cdef object v
    cdef object kv
    cdef Py_ssize_t[:, :] arr_view

    #Preprocess connections to be able to release gil below
    # Get the cumulative size of all values
    cdef long total = 0
    cdef long[:] keys_lens = np.empty(len(rconnections)*2, dtype=np.long)
    for v in rconnections.values():
        total += len(v)

    cdef long[:] values = np.empty(total, dtype=np.long)
    cdef long vpos = 0
    i = 0
    for kv in sorted(rconnections.items()):
        keys_lens[i] = kv[0]
        v = kv[1]
        len_v = len(v)
        keys_lens[i+1] = len_v
        i += 2

        # copy values
        for ind in range(len_v):
            values[vpos+ind] = v[ind]
        vpos += len_v

    arr_view = np.empty((2, total), dtype=np.intp)
    vpos = 0
    with nogil:
        for kl in range(0, keys_lens.shape[0]//2, 2):
            k = keys_lens[kl]
            len_v = keys_lens[kl+1]
            for ind in range(vpos, vpos+len_v):
                arr_view[0, ind] = k
                arr_view[1, ind] = bisect_left_long(data_index, values[ind-vpos])
            vpos += len_v

    return arr_view
