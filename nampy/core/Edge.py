# Note hyperedges aren't supported yet.
from .Object import Object
from .Node import Node
from .DictList import DictList

    
class Edge(Object):
    """

    """
    
    def __init__(self, node_pair_list, the_weight = 1):
        """

        """
        if ((type(node_pair_list[0]) == Node) & (type(node_pair_list[1]) == Node)):
            id = node_pair_list[0].id + '_' + node_pair_list[1].id
            Object.__init__(self, id)
            self._nodes = DictList(node_pair_list)
            self.weight = the_weight
            self._model = None
     
