from abc import ABC, abstractmethod



# -----------------------------------------------------------------------------
# Abstract Compute Class:
# Define all slots and pass function definitions to child classes
# -----------------------------------------------------------------------------
class AbstractCompute(ABC):
    """
    This just defines all of the slots that are be used by child classes.
    These need to be defined in a parent class so each child class can be
    combined into a single DataAssimilation object without getting a 
    'multiple-inheritance' error.
    """
    __slots__ = ["_reaches_ordered_bytw","_results",]
    
    def __init__(self,):
        """
        Run subnetworking pre-processing, then computing.
        """
        self._subset_domain()
        
        self._route()
        
    @property
    def get_output(self,):
        return self._results
    
    @abstractmethod
    def _subset_domain(self,):
        pass
    
    @abstractmethod
    def _route(self,):
        pass


# -----------------------------------------------------------------------------
# Compute class definitions:
#   1. serial
#   2. by_network
#   3. by_subnetwork_jit
#   4. by_subnetwork_jit_clustered: 
# -----------------------------------------------------------------------------
class serial(AbstractCompute):
    def __init__(self):
        """
        Serial compute class.
        """
        super().__init__()
    
    def _subset_domain(self,):
        #TODO Define subsetting method for serial
        self._reaches_ordered_bytw = {}
    
    def _route(self,):
        #TODO Define routing compute method for serial
        self._results = []
    
    
class by_network(serial):
    def __init__(self):
        """
        By Network compute class.
        
        #NOTE I think this can be a subclass of serial. It
        # just needs to route networks in parallel. -shorvath.
        """
        super().__init__()
    
    def _subset_domain(self,):
        #TODO Define subsetting method for by-network
        self._reaches_ordered_bytw = {}
    
    def _route(self,):
        #TODO Define routing compute method for by-network
        self._results = []
    
    
class by_subnetwork_jit(by_network):
    def __init__(self):
        """
        By Network JIT compute class.
        
        #NOTE I think this can be a subclass of by_network. It
        # just needs a couple extra steps to handle 'order',
        # e.g., 'flowveldepth_interorder'. -shorvath.
        """
        super().__init__()
    
    def _subset_domain(self,):
        #TODO Define subsetting method for by-network-jit
        self._reaches_ordered_bytw = {}
    
    def _route(self,):
        #TODO Define routing compute method for by-network-jit
        self._results = []
    
    
class by_subnetwork_jit_clustered(by_subnetwork_jit):
    def __init__(self):
        """
        By Network JIT Clustered compute class.
        
        #NOTE I think this can be a subclass of by_subnetwork_jit. It
        # just needs a couple extra steps to cluster subnetworks. -shorvath.
        """
        super().__init__()
    
    def _subset_domain(self,):
        #TODO Define subsetting method for by-network-jit-clustered
        self._reaches_ordered_bytw = {}
    
    def _route(self,):
        #TODO Define routing compute method for by-network-jit-clustered
        self._results = []