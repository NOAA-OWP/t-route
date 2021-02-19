from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.stdio cimport printf
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


############ Other Reservoir Interface ############
ctypedef enum compute_type: RESERVOIR_LP

cdef class compute_kernel_lp:
  """
    An encapsulation of a compute kernel that thread will execute
  """

  cdef int num_inputs;
  cdef int num_outputs;
  #cdef double* inputs;
  #cdef double* outputs;
  cdef float* inputs;
  cdef float* outputs;

  def __init__(self, num_inputs, num_outputs):
    """
      Allocate the input and output arrays
    """
    #self.inputs = <double*> PyMem_Malloc(num_inputs * sizeof(double))
    #self.outputs = <double*> PyMem_Malloc(num_outputs * sizeof(double))
    self.inputs = <float*> PyMem_Malloc(num_inputs * sizeof(float))
    self.outputs = <float*> PyMem_Malloc(num_outputs * sizeof(float))

  def __dealloc__(self):
    PyMem_Free(self.inputs)
    PyMem_Free(self.outputs)

cdef class lp_kernel(compute_kernel_lp):
  """
    Subclass for computing LevelPool reservoir
  """
  #Non python accessible attribute
  cdef void* lp_handle;

  #Reservoir attributes
  cdef float water_elevation
  cdef float lake_area
  cdef float weir_elevation
  cdef float weir_coefficient
  cdef float weir_length
  cdef float dam_length
  cdef float orifice_elevation
  cdef float orifice_coefficient
  cdef float orifice_area
  cdef float max_depth
  cdef int lake_number

  def __init__(self, num_inputs, num_outputs,
              water_elevation, lake_area, weir_elevation, weir_coefficient, weir_length,
              dam_length, orifice_elevation, orifice_coefficient, orifice_area,
              max_depth, lake_number):
    super().__init__(num_inputs, num_outputs)
    self.water_elevation = water_elevation
    self.lake_area = lake_area
    self.weir_elevation = weir_elevation
    self.weir_coefficient = weir_coefficient
    self.weir_length = weir_length
    self.dam_length = dam_length
    self.orifice_elevation = orifice_elevation
    self.orifice_coefficient = orifice_coefficient
    self.orifice_area = orifice_area
    self.max_depth = max_depth
    self.lake_number = lake_number

    with nogil:
      self.lp_handle = get_lp_handle()
      init_lp(self.lp_handle, &self.water_elevation, &self.lake_area,
                   &self.weir_elevation, &self.weir_coefficient, &self.weir_length,
                   &self.dam_length, &self.orifice_elevation, &self.orifice_coefficient,
                   &self.orifice_area, &self.max_depth, &self.lake_number)

  def __dealloc__(self):
    free_lp(self.lp_handle)

  cpdef float run(self, float inflow, float lateral_inflow, float routing_period):
    cdef float outflow = 0.0

    with nogil:
      run_lp(self.lp_handle, &inflow, &lateral_inflow, &self.water_elevation, &outflow, &routing_period)
      #printf("outflow: %f\n", outflow)
      return outflow

  cpdef float get_water_elevation(self):
    cdef float water_elevation

    with nogil:
      water_elevation = self.water_elevation
      return water_elevation

  cpdef int get_lake_number(self):
    cdef int lake_number

    with nogil:
      lake_number = self.lake_number
      return lake_number
