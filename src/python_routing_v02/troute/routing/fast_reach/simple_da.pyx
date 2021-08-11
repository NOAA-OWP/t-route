from libc.math cimport exp, isnan
# from libc.stdio cimport printf

cdef float simple_da_with_decay(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil:
    """
    Given a modeled value, last valid observation,
    time since that observation, and an exponential
    decay_coefficient, compute the 'nudged' value
    """

    cdef float da_weight, da_shift, da_weighted_shift
    da_weight = exp(minutes_since_last_valid/-decay_coeff)  # TODO: This could be pre-calculated knowing when obs finish relative to simulation time
    # TODO: we need to be able to export these values to compute the 'Nudge'
    # One possibility would be to return only the nudge from this function...
    da_shift = last_valid_obs - model_val
    da_weighted_shift = da_shift * da_weight
    # printf("%g %g %g %g %g --> %g\n", minutes_since_last_valid, decay_coeff, da_shift, da_weighted_shift, model_val, model_val + da_weighted_shift)
    return model_val + da_weighted_shift
