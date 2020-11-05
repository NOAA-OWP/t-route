cimport numpy as np

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
  cdef float qdp, velp, depthp

cdef class MC_Reach():
  """
    A muskingcung reach -> collection of ordered MC_Segments
  """
  cdef readonly list segments
  cdef readonly np.ndarray upstream_ids
  cdef route(self)
