cimport numpy as np
"""
FIXME add some significant inline documentation
"""
from troute.network.reach cimport Reach, compute_type

############ Other Reservoir Interface ############
cdef void run(_Reach* reach, float inflow, float lateral_inflow, float routing_period, float* outflow,  float* water_elevation) nogil

cdef extern from "rfc_structs.h":
  ctypedef struct _MC_RFC:
    int lake_number
    float dam_length, area, max_depth
    float orifice_area, orifice_coefficient, orifice_elevation
    float weir_coefficient, weir_elevation, weir_length
    float initial_fractional_depth, water_elevation
    int reservoir_type
    bytes reservoir_parameter_file
    bytes start_date
    bytes time_series_path
    int forecast_lookback_hours
  ctypedef struct _Reach:
    pass

cdef class MC_RFC(Reach):
  """
  """
  cpdef (float,float) run(self, float inflow, float lateral_inflow, float routing_period)