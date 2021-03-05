cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

from troute.network.reach cimport compute_type, _Reach
"""
Externally defined symbols
"""

cdef extern from "levelpool_structs.c":
  void init_levelpool_reach(_Reach* reach, int lake_number,
                            float dam_length, float area, float max_depth,
                            float orifice_area, float orifice_coefficient, float orifice_elevation,
                            float weir_coefficient, float weir_elevation, float weir_length,
                            float initial_fractional_depth, float water_elevation
  )
  void free_levelpool_reach(_Reach* reach)
  void route(_Reach* reach, float routing_period, float inflow, float lateral_inflow, float* outflow,  float* water_elevation) nogil

cdef void run(_Reach* reach, float inflow, float lateral_inflow, float routing_period, float* outflow,  float* water_elevation) nogil:
    route(reach, inflow, lateral_inflow, routing_period, outflow, water_elevation)

cdef class MC_Levelpool(Reach):
  """
    MC_Reservoir is a subclass of MC_Reach_Base_Class
  """

  def __init__(self, long id, int lake_number, long[::1] upstream_ids, args):
    """
      Construct the kernel based on passed parameters,
      which only constructs the parent class
    """
    super().__init__(id, upstream_ids, compute_type.RESERVOIR_LP)
    #init the backing struct, pass a dam_length of 10.0 for now
    #pass a negative water elevation, which causes the init to use the wrf hydro equation
    #TODO put in __calloc__
    init_levelpool_reach(&self._reach, lake_number,
                         10.0, args[0], args[1],
                         args[2], args[3], args[4],
                         args[5], args[6], args[7],
                         args[8], args[10])
    """
    self.lake_number = lake_number
    #TODO: Need new Lake Parm file, which now has dam_length
    #dam_length = wbody_parameters[wbody_index,1]
    #Setting default dam_length to 10
    self.dam_length = 10.0
    self.area = args[0]
    self.max_depth = args[1]
    self.orifice_area = args[2]
    self.orifice_coefficient = args[3]
    self.orifice_elevation  =  args[4]
    self.weir_coefficient = args[5]
    self.weir_elevation = args[6]
    self.weir_length = args[7]
    self.initial_fractional_depth  = args[8]
    #TODO: Read Water Elevation from Restart. Use below equation if no restart.
    #Equation below is used in wrf-hydro
    self.water_elevation = self.orifice_elevation + ((self.max_depth - self.orifice_elevation) * self.initial_fractional_depth)
    """
    """
    #Initialize level pool reservoir object
    with nogil:
      self.lp_handle = get_lp_handle()
      init_lp(self.lp_handle, &self.water_elevation, &self.area,
                   &self.weir_elevation, &self.weir_coefficient, &self.weir_length,
                   &self.dam_length, &self.orifice_elevation, &self.orifice_coefficient,
                   &self.orifice_area, &self.max_depth, &self.lake_number)
    #print(<int>self.lp_handle)
    """

  def __dealloc__(self):
    """

    """
    free_levelpool_reach(&self._reach)
    #free_lp(self.lp_handle)

  cpdef (float,float) run(self, float inflow, float lateral_inflow, float routing_period):
    cdef float outflow = 0.0
    cdef float water_elevation = 0.0
    with nogil:
      route(&self._reach, inflow, lateral_inflow, routing_period, &outflow,  &water_elevation)
      #printf("outflow: %f\n", outflow)
      return outflow, water_elevation#, self.water_elevation
