"""
See simple test example [here](https://colab.research.google.com/drive/1wrCl6bTQWNVMQj0U54vVdXMOx4mDf5_d#scrollTo=cOZV2reX1zDJ)<br><br>
NHD analysis from GIS (based on channels_nwm_v12_routeLink)
Total number of channel segments: 2,729,077
Order 10 terminal segments:   1  (total: 441)             1 Running total of terminal segments
Order 9 terminal segments:    2  (total: 3,334)           3
Order 8 terminal segments:    4  (total: 7,495)           7
Order 7 terminal segments:   31  (total: 20,956)         38
Order 6 terminal segments:   87  (total: 48,275)        125
Order 5 terminal segments:  230  (total: 91,386)        355
Order 4 terminal segments:  987  (total: 170,417)     1,342
Order 3 terminal segments: 1972  (total: 315,289)     3,314
Order 2 terminal segments: 5788  (total: 598,992)     9,102
Order 1 terminal segments: 5611  (total: 1,472,492)  14,713

##Define functions
"""
import time
class SuperNetwork():
    '''
    Any collection of networks, above the level of strict hydraulic connectivity,
    begins to be a human artifact and potentially arbitrary.
    Nonetheless, there may be cases where this is sensible, such as in the case
    of the Missississippi River SuperNetwork (with at least Missouri, Red,
    Arkansas, Illinois, and Ohio as contained networks.) Likewise, one could group
    all continental rivers, or rivers on either side of the continental divide,
    etc. In any case, this designation will almost always be manually applied,
    and in the test case, the SuperNetwork is simply the collection of all networks
    in the test dataset.
    '''
    def __init__(self):
        self.networkCollection = []

class Network():
    def __init__(self, networkID = None):
        self.networkID = networkID
        self.reachCollection = []
        self.nexusCollection = []
        head_segment = []
        trunk_segment = None
        head_reach = []
        trunk_reach = None
        head_nexus = []
        trunk_segement = None


