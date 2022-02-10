cdef (float, float, float, float) simple_da(
    const float timestep,
    const float routing_period,
    const float decay_coeff,
    const float gage_maxtimestep,
    const float target_val,
    const float model_val,
    float lastobs_time,
    float lastobs_val,
    bint da_check_gage=*,
) nogil


cdef float simple_da_with_decay(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil


cdef float obs_persist_shift(
    const float last_valid_obs,
    const float model_val,
    const float minutes_since_last_valid,
    const float decay_coeff,
) nogil
