cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free
from troute.network.reach cimport compute_type, _Reach

cdef extern from "mc_reach_structs.c":
  void init_mc_reach(_Reach* reach, int num_segments)
  void free_mc_reach(_Reach* reach)
  void set_mc_segment(_Reach* reach, int index, long id,
      float dt, float dx, float bw, float tw, float twcc,
      float n, float ncc, float cs, float s0,
      float qdp, float velp, float depthp)
  _MC_Segment get_mc_segment(_Reach* reach, int index) nogil

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

cdef class MC_Reach(Reach):
  """
    A muskingcung reach and subclass of MC_Reach_Base_Class
    TODO keep track of topology, subclass Reach
  """
  #TODO document these attributes. defined in pxd, commented here for reference
  #cdef readonly list segments
  #upstream_ids are stored as an actual ndarray so they can be used as an indexer
  #in other ndarrays.  This may need to change
  #cdef readonly np.ndarray upstream_ids

  def __init__(self, segments, long[::1] upstream_ids):
    """
      segments: ORDERED list of segments that make up the reach
    """
    #for now mc reaches aren't explicity identified, their segments are pass -1 for id
    super().__init__(-1, upstream_ids, compute_type.MC_REACH)

    self._num_segments = len(segments)
    init_mc_reach(&self._reach, self._num_segments)
    for i, s in enumerate(segments):
      #TODO  what about segment id???
      set_mc_segment(&self._reach, i, s.id,
                     s.dt, s.dx, s.bw, s.tw,
                     s.twcc, s.n, s.ncc, s.cs,
                     s.s0, s.qdp, s.velp, s.depthp)

    #self._segments = segments
    self.num_segments = self._num_segments

  #TODO implement __cinit__ for init_mc_reach?
  def __dealloc__(self):
    """

    """
    free_mc_reach(&self._reach)

  cdef route(self):
    """
      TODO implement?
    """
    #do a few things with gil
    #with nogil:
    #  pass
    pass

  def __getitem__(self, index):
    #TODO implement slicing, better errors
    if(index < 0):
      index = index + self._num_segments
    if(index > -1 and index <self._num_segments):
      return get_mc_segment(&self._reach, index)
    else:
      raise(IndexError)

  def __len__(self):
    return self._num_segments
