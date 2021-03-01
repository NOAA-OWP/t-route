cimport numpy as np
"""
FIXME add some significant inline documentation
"""
from troute.network.reach cimport Reach, compute_type

############ Other Reservoir Interface ############
cdef void run(_Reach* reach, float routing_period, float inflow, float lateral_inflow, float* outflow,  float* water_elevation) nogil

cdef extern from "levelpool_structs.h":
  ctypedef struct _MC_Levelpool:
    int lake_number
    float dam_length, area, max_depth
    float orifice_area, orifice_coefficient, orifice_elevation
    float weir_coefficient, weir_elevation, weir_length
    float initial_fractional_depth, water_elevation
  ctypedef struct _Reach:
    pass

cdef class MC_Levelpool(Reach):
  """

  """
  cpdef (float,float) run(self, float inflow, float lateral_inflow, float routing_period)
