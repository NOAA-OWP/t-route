import cython
#from libc.stdio cimport printf

from .fortran_wrappers cimport c_muskingcungenwm

@cython.boundscheck(False)
cdef void muskingcunge(float dt,
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
        QVD *rv) nogil:

    cdef:
        float qdc = 0.0
        float depthc = 0.0
        float velc = 0.0
        float ck = 0.0
        float cn = 0.0
        float X = 0.0

    #printf("reach.pyx before %3.9f\t", depthc)
    c_muskingcungenwm(
        &dt,
        &qup,
        &quc,
        &qdp,
        &ql,
        &dx,
        &bw,
        &tw,
        &twcc,
        &n,
        &ncc,
        &cs,
        &s0,
        &velp,
        &depthp,
        &qdc,
        &velc,
        &depthc,
        &ck,
        &cn,
        &X)
    #printf("reach.pyx after %3.9f\t", depthc)

    rv.qdc = qdc
    rv.depthc = depthc
    rv.velc = velc

    # to do: make these additional variable's conditional, somehow
    rv.ck = ck
    rv.cn = cn
    rv.X = X

cpdef dict compute_reach_kernel(float dt,
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
        float depthp):

    cdef QVD rv
    cdef QVD *out = &rv

    muskingcunge(
        dt,
        qup,
        quc,
        qdp,
        ql,
        dx,
        bw,
        tw,
        twcc,
        n,
        ncc,
        cs,
        s0,
        velp,
        depthp,
        out)

    return rv


cpdef long boundary_shape() nogil:
    return 2

cpdef long previous_state_cols() nogil:
    return output_buffer_cols()

cpdef long parameter_inputs_cols() nogil:
    return 13

cpdef long output_buffer_cols() nogil:
    return 3

@cython.boundscheck(False)
cpdef float[:,:] compute_reach(const float[:] boundary,
                                const float[:,:] previous_state,
                                const float[:,:] parameter_inputs,
                                float[:,:] output_buffer,
                                Py_ssize_t size=0) nogil:
    """
    Compute a reach

    Arguments:
        boundary: [qup, quc]
        previous_state: Previous state for each node in the reach [qdp, velp, depthp]
        parameter_inputs: Parameterization of the reach at node.
            qlat, dt, dx, bw, tw, twcc, n, ncc, cs, s0
        output_buffer: Current state [qdc, velc, depthc]

    """
    cdef QVD rv
    cdef QVD *out = &rv

    cdef:
        float dt, qlat, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp
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
    if parameter_inputs.shape[1] < 10:
        raise IndexError
    if output_buffer.shape[1] < 3:
        raise IndexError
    if previous_state.shape[1] < 3:
        raise IndexError

    cdef float qup = boundary[0]
    cdef float quc = boundary[1]

    for i in range(rows):
        qlat = parameter_inputs[i, 0]
        dt = parameter_inputs[i, 1]
        dx = parameter_inputs[i, 2]
        bw = parameter_inputs[i, 3]
        tw = parameter_inputs[i, 4]
        twcc = parameter_inputs[i, 5]
        n = parameter_inputs[i, 6]
        ncc = parameter_inputs[i, 7]
        cs = parameter_inputs[i, 8]
        s0 = parameter_inputs[i, 9]

        qdp = previous_state[i, 0]
        velp = previous_state[i, 1]
        depthp = previous_state[i, 2]

        muskingcunge(
                    dt,
                    qup,
                    quc,
                    qdp,
                    qlat,
                    dx,
                    bw,
                    tw,
                    twcc,
                    n,
                    ncc,
                    cs,
                    s0,
                    velp,
                    depthp,
                    out)

        output_buffer[i, 0] = quc = out.qdc
        output_buffer[i, 1] = out.velc
        output_buffer[i, 2] = out.depthc

        qup = qdp
    return output_buffer
