import cython

from fortran_wrappers cimport c_levelpool_physics


@cython.boundscheck(False)
cdef void levelpool_physics(float dt,
        float qi0,
        float qi1,
        float ql,
        float ar,
        float we,
        float maxh,
        float wc,
        float wl,
        float dl,
        float oe,
        float oc,
        float oa,
        float H0,
        QH *rv) nogil:
    cdef:
        float H1 = 0.0
        float qo1 = 0.0

    c_levelpool_physics(
        &dt,
        &qi0,
        &qi1,
        &ql,
        &ar,
        &we,
        &maxh,
        &wc,
        &wl,
        &dl,
        &oe,
        &oc,
        &oa,
        &H0,
        &H1,
        &qo1)
    rv.reslevel = H1
    rv.resoutflow = qo1


cpdef dict compute_reservoir_kernel(float dt,
        float qi0,
        float qi1,
        float ql,
        float ar,
        float we,
        float maxh,
        float wc,
        float wl,
        float dl,
        float oe,
        float oc,
        float oa,
        float H0):

    cdef QH rv
    cdef QH *out = &rv

    levelpool_physics(dt,
        qi0,
        qi1,
        ql,
        ar,
        we,
        maxh,
        wc,
        wl,
        dl,
        oe,
        oc,
        oa,
        H0,
        out)

    return rv


cpdef long boundary_shape() nogil:
    return 2

cpdef long previous_state_cols() nogil:
    return output_buffer_cols()

cpdef long parameter_inputs_cols() nogil:
    return 12

cpdef long output_buffer_cols() nogil:
    return 3


@cython.boundscheck(False)
cpdef float[:,:] compute_reservoir(const float[:] boundary,
                                    const float[:,:] previous_state,
                                    const float[:,:] parameter_inputs,
                                    float[:,:] output_buffer,
                                    Py_ssize_t size=0) nogil:
    """
    Compute a reservoir

    Arguments:
        boundary: [qi0, qi1]
        previous_state: Previous state for each node in the reach [qdp, velp, depthp]
        parameter_inputs: Parameterization of the reach at node.
            dt,ql,ar,we,maxh,wc,wl,dl,oe,oc,oa,H0
        output_buffer: Current state [H1, qo0]
    """
    cdef QH rv
    cdef QH *out = &rv

    cdef:
        float dt, ql, ar, we, maxh, wc, wl, dl, oe, oc, oa, H0
        Py_ssize_t i, rows

    # check that previous state, parameter_inputs and output_buffer all have same axis 0
    if size > 0:
        rows = size
        if (parameter_inputs.shape[0] < rows
                or output_buffer.shape[0] < rows
                or previous_state.shape[0] < rows):
            raise ValueError(f"axis 0 is not long enough for {size}")
    else:
        rows = previous_state.shape[0]
        if rows != parameter_inputs.shape[0] or rows != output_buffer.shape[0]:
            raise ValueError("axis 0 of input arguments do not agree")

    # check bounds
    if boundary.shape[0] < 2:
        raise IndexError
    if parameter_inputs.shape[1] < 12:
        raise IndexError
    if output_buffer.shape[1] < 2:
        raise IndexError

    cdef float qi0 = boundary[0]
    cdef float qi1 = boundary[1]

    for i in range(rows):
        dt = parameter_inputs[i, 0]
        ql = parameter_inputs[i, 1]
        ar = parameter_inputs[i, 2]
        we = parameter_inputs[i, 3]
        maxh = parameter_inputs[i, 4]
        wc = parameter_inputs[i, 5]
        wl = parameter_inputs[i, 6]
        dl = parameter_inputs[i, 7]
        oe = parameter_inputs[i, 8]
        oc = parameter_inputs[i, 9]
        oa = parameter_inputs[i, 10]
        H0 = parameter_inputs[i, 11]

        levelpool_physics(dt,
            qi0,
            qi1,
            ql,
            ar,
            we,
            maxh,
            wc,
            wl,
            dl,
            oe,
            oc,
            oa,
            H0,
            out)

        output_buffer[i, 0] = out.reslevel
        output_buffer[i, 1] = out.resoutflow

    return output_buffer
