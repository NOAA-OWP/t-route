cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

cdef class Reach:
  """
    Encapsulation of basic reach properties and functionality
  """
  cdef readonly int id;
  #Keep a python list of ids only for pre/post processing and diagnostics
  cdef readonly long[::1] to_ids

  def __init__(self, int id, long[::1] to_ids):
    """
      Initialize the base properties
    """
    self.id = id
    self.to_ids = to_ids

cdef class Segment():
  """
    A Single routing segment
  """
  #cdef long id
  #cdef long upstream_id
  #cdef long downstream_id

  def __init__(self, id, upstream, downstream):
    """

    """
    self.id = id
    self.upstream_id = upstream
    self.downstream_id = downstream

cdef class MC_Segment(Segment):
  """
    A single segment that compes
  """

  def __init__(self, id, dt, dx, bw, tw, twcc, n, ncc, cs, s0, qdp, velp, depthp):
    """
      construct the kernel based on passed parameters
      TODO record segement topology better (i.e. up/down stream id)
    """
    super().__init__(id, -1, -1)
    self.dt = dt
    self.dx = dx
    self.bw = bw
    self.tw = tw
    self.twcc = twcc
    self.n = n
    self.ncc = ncc
    self.cs = cs
    self.s0 = s0
    self.qdp = qdp
    self.velp = velp
    self.depthp = depthp

cdef class MC_Reach():
  """
    A muskingcung reach
    TODO keep track of topology, subclass Reach
  """
  #TODO document these attributes. defined in pxd, commented here for reference
  #cdef readonly list segments
  #upstream_ids are stored as an actual ndarray so they can be used as an indexer
  #in other ndarrays.  This may need to change
  #cdef readonly np.ndarray upstream_ids

  def __init__(self, segments, long[::1] upstream_ids):
    self._num_segments = len(segments)
    self._num_upstream_ids = len(upstream_ids)

    self._segments = <_MC_Segment*> malloc(sizeof(_MC_Segment)*(self._num_segments))
    for i, segment in enumerate(segments):
      self._segments[i].id = segment.id
      self._segments[i].dt = segment.dt
      self._segments[i].dx = segment.dx
      self._segments[i].bw = segment.bw
      self._segments[i].tw = segment.tw
      self._segments[i].twcc = segment.twcc
      self._segments[i].n = segment.n
      self._segments[i].ncc = segment.ncc
      self._segments[i].cs = segment.cs
      self._segments[i].s0 = segment.s0
      self._segments[i].qdp = segment.qdp
      self._segments[i].velp = segment.velp
      self._segments[i].depthp = segment.depthp

    self.segments = segments
    self.upstream_ids = np.asarray(upstream_ids)
    #init the C struct
    self._reach._segments = self._segments
    self._reach._num_segments = self._num_segments
    self.num_segments = self._num_segments
    self._reach._upstream_ids = <long*>malloc(sizeof(long)*self._num_upstream_ids)
    self._reach._num_upstream_ids = self._num_upstream_ids
    for i in range(self._num_upstream_ids):
      self._reach._upstream_ids[i] = upstream_ids[i]

  def __dealloc__(self):
    """

    """
    free(self._segments)
    free(self._reach._upstream_ids)

  cdef route(self):
    """
      TODO implement?
    """
    #do a few things with gil
    #with nogil:
    #  pass
    pass
