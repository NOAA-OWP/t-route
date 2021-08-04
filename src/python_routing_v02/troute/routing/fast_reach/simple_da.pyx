from libc.math cimport exp, isnan

cdef float simple_da_with_decay(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil:

    cdef float da_weight, da_shift, da_weighted_shift, replacement_value
    da_weight = exp(minutes_since_last_valid/-decay_coeff)  # TODO: This could be pre-calculated knowing when obs finish relative to simulation time
    # replacement_value = f(lastobs_value, da_weight)  # TODO: we need to be able to export these values to compute the 'Nudge'
    da_shift = last_valid_obs - model_val
    da_weighted_shift = da_shift * da_weight
    return model_val + da_weighted_shift
