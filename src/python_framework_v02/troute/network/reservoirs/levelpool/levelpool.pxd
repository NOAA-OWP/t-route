cimport numpy as np
"""
FIXME add some significant inline documentation
"""
from troute.network.reach cimport MC_Reach_Base_Class

############ Other Reservoir Interface ############
ctypedef enum compute_type: RESERVOIR_LP

cdef class MC_Levelpool(MC_Reach_Base_Class):
  """

  """
  cdef int lake_number
  cdef float dam_length, area, max_depth
  cdef float orifice_area, orifice_coefficient, orifice_elevation
  cdef float weir_coefficient, weir_elevation, weir_length
  cdef float initial_fractional_depth, water_elevation
  #Initialize level pool reservoir object
  cdef void* lp_handle;
  cpdef (float,float) run(self, float inflow, float lateral_inflow, float routing_period)
  cpdef float get_water_elevation(self)
