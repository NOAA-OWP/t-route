cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

from troute.network.reach cimport compute_type, _Reach
"""
Externally defined symbols
"""

cdef extern from "rfc_structs.c":
  void init_rfc_reach(_Reach* reach, int lake_number,
                            float dam_length, float area, float max_depth,
                            float orifice_area, float orifice_coefficient, float orifice_elevation,
                            float weir_coefficient, float weir_elevation, float weir_length,
                            float initial_fractional_depth, float water_elevation,
                            int reservoir_type, char *reservoir_parameter_file, char *start_date,
                            char *time_series_path, int forecast_lookback_hours
  )
  void free_rfc_reach(_Reach* reach)
  void route(_Reach* reach, float routing_period, float inflow, float lateral_inflow, float* outflow,  float* water_elevation) nogil

cdef void run(_Reach* reach, float inflow, float lateral_inflow, float routing_period, float* outflow,  float* water_elevation) nogil:
    route(reach, inflow, lateral_inflow, routing_period, outflow, water_elevation)

cdef class MC_RFC(Reach):
  """
    MC_Reservoir is a subclass of MC_Reach_Base_Class
  """

  def __init__(self, long id, int lake_number, long[::1] upstream_ids, args):
    """
      Construct the kernel based on passed parameters,
      which only constructs the parent class

      Params:
        id: long
          unique identity of the reach this reservoir represents
        lake_number: int TODO (long?)
          WRF_Hydro lake number of this reservoir
        upstream_ids: array[long]
          buffer/array of upstream identifiers which contribute flow to this reservoir
        args: list
          the levelpool paramters ordered as follows:
            area = args[0]
            max_depth = args[1]
            orifice_area = args[2]
            orifice_coefficient = args[3]
            orifice_elevation  =  args[4]
            weir_coefficient = args[5]
            weir_elevation = args[6]
            weir_length = args[7]
            initial_fractional_depth  = args[8]
            reservoir_type =  args[11]
            reservoir_parameter_file = args[12]
            start_date = args[13]
            time_series_path = args[14]
            forecast_lookback_hours = args[15]
    """
    super().__init__(id, upstream_ids, compute_type.RESERVOIR_RFC)
    # Note Some issues with __calloc__:
    # The python type isn't guaranteed to be properly constructed, so cannot depend on super class being constructured.
    # Thus I don't think we can put these C init functions in __calloc__, at least not in all cases.
    # init the backing struct, pass a dam_length of 10.0 for now
    #init_hybrid_reach(&self._reach, lake_number,
    #                     10.0, args[0], args[1],
    #                     args[2], args[3], args[4],
    #                     args[5], args[6], args[7],
    #                     args[8], args[10], args[11],
    #                     args[12], args[13], args[14],
    #                     args[15], args[16], args[17])

    #Setting default dam_length to 10
    dam_length = 10.0
    area = args[0]
    max_depth = args[1]
    orifice_area = args[2]
    orifice_coefficient = args[3]
    orifice_elevation  =  args[4]
    weir_coefficient = args[5]
    weir_elevation = args[6]
    weir_length = args[7]
    initial_fractional_depth  = args[8]
    reservoir_type =  args[11]
    reservoir_parameter_file = args[12]
    start_date = args[13]
    time_series_path = args[14]
    forecast_lookback_hours = args[15]

    #Check lengths of input strings to ensure that they do not exceed buffer size
    if (reservoir_parameter_file > 256):
       raise ValueError("reservoir_parameter_file path is too large. Length must be less than or equal to 256 characters.")

    # Note Some issues with __calloc__:
    # The python type isn't guaranteed to be properly constructed, so cannot depend on super class being constructured.
    # Thus I don't think we can put these C init functions in __calloc__, at least not in all cases.
    # init the backing struct, pass a dam_length of 10.0 for now
    init_rfc_reach(&self._reach, lake_number,
                         10.0, args[0], args[1],
                         args[2], args[3], args[4],
                         args[5], args[6], args[7],
                         args[8], args[10], args[11],
                         args[12], args[13], args[14],
                         args[15])

  def __dealloc__(self):
    """
    """
    free_rfc_reach(&self._reach)

  cpdef (float,float) run(self, float inflow, float lateral_inflow, float routing_period):
    """
      Run the rfc routing function
      Params:
        inflow: float
          inflow into the reservoir
        lateral_inflow: float
          lateral flows into the reservoir
        routing_period: float
          amount of time to simulatie reservoir operation for, outflow if valid until this time
      Return:
        outflow: float
          flow rate out of the reservoir valid for routing_period seconds
        water_elevation:
          reservoir water surface elevation after routing_period seconds
    """
    cdef float outflow = 0.0
    cdef float water_elevation = 0.0
    with nogil:
      route(&self._reach, inflow, lateral_inflow, routing_period, &outflow,  &water_elevation)
      #printf("outflow: %f\n", outflow)
      return outflow, water_elevation#, self.water_elevation

  @property
  def water_elevation(self):
    return self._reach.reach.rfc.water_elevation

  @property
  def lake_area(self):
    return self._reach.reach.rfc.area

  @property
  def weir_elevation(self):
    return self._reach.reach.rfc.weir_elevation

  @property
  def weir_coefficient(self):
    return self._reach.reach.rfc.weir_coefficient

  @property
  def weir_length(self):
    return self._reach.reach.rfc.weir_length

  @property
  def dam_length(self):
    return self._reach.reach.rfc.dam_length

  @property
  def orifice_elevation(self):
    return self._reach.reach.rfc.orifice_elevation

  @property
  def orifice_area(self):
    return self._reach.reach.rfc.orifice_area

  @property
  def max_depth(self):
    return self._reach.reach.rfc.max_depth

  @property
  def lake_number(self):
    return self._reach.reach.rfc.lake_number

  @property
  def initial_fractional_depth(self):
    return self._reach.reach.rfc.initial_fractional_depth

  @property
  def reservoir_type(self):
    return self._reach.reach.rfc.reservoir_type

  @property
  def reservoir_parameter_file(self):
    return self._reach.reach.rfc.reservoir_parameter_file

  @property
  def start_date(self):
    return self._reach.reach.rfc.start_date

  @property
  def time_series_path(self):
    return self._reach.reach.rfc.time_series_path

  @property
  def forecast_lookback_hours(self):
    return self._reach.reach.rfc.forecast_lookback_hours
