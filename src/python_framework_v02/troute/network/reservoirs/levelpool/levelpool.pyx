cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

"""
Externally defined symbols
"""
############ Level Pool Reservoir Interface ############
cdef extern void* get_lp_handle() nogil;

cdef extern void init_lp(void* handle, float *water_elevation, float *lake_area, float *weir_elevation,
                    float *weir_coefficient, float *weir_length, float *dam_length, float *orifice_elevation,
                    float *orifice_coefficient, float *orifice_area, float *max_depth, int *lake_number) nogil;

cdef extern void run_lp(void* handle, float *inflow, float *lateral_inflow,
                    float *water_elevation, float *outflow, float *routing_period) nogil;

cdef extern void free_lp(void* handle);

cdef class MC_Levelpool(MC_Reach_Base_Class):
  """
    MC_Reservoir is a subclass of MC_Reach_Base_Class
  """

  def __init__(self, lake_number, long[::1] upstream_ids, args):
    """
      Construct the kernel based on passed parameters,
      which only constructs the parent class
    """
    super().__init__(upstream_ids)
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
    #Initialize level pool reservoir object
    with nogil:
      self.lp_handle = get_lp_handle()
      init_lp(self.lp_handle, &self.water_elevation, &self.area,
                   &self.weir_elevation, &self.weir_coefficient, &self.weir_length,
                   &self.dam_length, &self.orifice_elevation, &self.orifice_coefficient,
                   &self.orifice_area, &self.max_depth, &self.lake_number)
    #print(<int>self.lp_handle)

  def __dealloc__(self):
    """

    """
    free(self._reach._upstream_ids)
    free_lp(self.lp_handle)

  cpdef (float,float) run(self, float inflow, float lateral_inflow, float routing_period):
    cdef float outflow = 0.0
    with nogil:
      run_lp(self.lp_handle, &inflow, &lateral_inflow, &self.water_elevation, &outflow, &routing_period)
      #printf("outflow: %f\n", outflow)
      return outflow, self.water_elevation#, self.water_elevation

  cpdef float get_water_elevation(self):
    #cdef float water_elevation

    with nogil:
      #water_elevation = self.water_elevation
      return self.water_elevation
