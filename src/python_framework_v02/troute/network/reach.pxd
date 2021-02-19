cimport numpy as np
"""
FIXME add some significant inline documentation
"""

cdef class Segment():
  """
    A Single routing segment
  """
  cdef readonly long id
  cdef long upstream_id
  cdef long downstream_id

cdef class MC_Segment(Segment):
  """
    A muskingcung segment
  """
  #TODO document these attributes
  cdef readonly float dt, dx, bw, tw, twcc, n, ncc, cs, s0
  cdef readonly float qdp, velp, depthp

cdef struct _MC_Segment:

  long id
  float dt, dx, bw, tw, twcc, n, ncc, cs, s0
  float qdp, velp, depthp

cdef struct _MC_Reach:
  _MC_Segment* _segments
  int _num_segments
  long* _upstream_ids
  int _num_upstream_ids

cdef class MC_Reach_Base_Class():
  cdef int _num_upstream_ids
  cdef _MC_Reach _reach
  cdef readonly np.ndarray upstream_ids

cdef class MC_Reservoir(MC_Reach_Base_Class):
  pass

cdef class MC_Reach(MC_Reach_Base_Class):
  """
    A muskingcung reach -> collection of ordered MC_Segments
  """
  cdef _MC_Segment* _segments
  cdef int _num_segments # C only accessible
  cdef readonly num_segments #Python accessible, readonly outside class
  cdef readonly list segments
  cdef route(self)
