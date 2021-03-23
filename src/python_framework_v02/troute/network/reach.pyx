cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "reach_structs.c":
  void init_reach(_Reach* reach, long* upstream_ids, int num_upstream_ids, int reach_type)
  void free_reach(_Reach* reach)
  long* get_upstream_ids(_Reach* reach)

cdef class Reach:
  """
    Encapsulation of basic reach properties and functionality
  """

  def __init__(self, long id, long[::1] upstream_ids, int reach_type):
    """
      Construct the kernel based on passed parameters
    """
    self.to_ids = np.ones(1, dtype=np.int64 )*-1
    self._num_upstream_ids = len(upstream_ids)

    if( len(upstream_ids) > 0 ):
      init_reach(&self._reach, &upstream_ids[0], len(upstream_ids), reach_type)
    else:
      init_reach(&self._reach, NULL, 0, reach_type)
    #(<_Reach>self._reach).id = id
    self._reach.id = id

  def __dealloc__(self):
    """
      deallocate memory used for reach
    """
    free_reach(&self._reach)

  @property
  def id(self):
    """
      A unique integer identifier of the reach, implies order (ids < self.id) are upstream reaches,
      (ids > self.id) are downstream reaches
    """
    return self._reach.id

  @property
  def upstream_ids(self):
    """
      List of upstream segment identifiers which contribute flow to this reach
    """
    cdef long[::1] ids
    if(self._num_upstream_ids == 0):
      return []
    else:
      ids = <long[:self._num_upstream_ids]> get_upstream_ids(&self._reach)
      return ids

cdef class Segment():
  """
    A Single routing segment
  """

  def __init__(self, id, upstream, downstream):
    """

    """
    self.id = id
    self.upstream_id = upstream
    self.downstream_id = downstream
