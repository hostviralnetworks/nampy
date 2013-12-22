from copy import deepcopy
import re
from .Object import Object
from .DictList import DictList

class Node(Object):
    """

    
    """

    def __init__(self, id):
        """

        """
        Object.__init__(self, id)
        # Object includes notes and annotation fields
        self._edges = DictList()
        # Associate the nodes with the network
        # This allows us to convert to monopartite
        # and can keep the nodetype as a string
        # rather than a reference
        self._network = None
        self._nodetype = None
        self.source = 0
        self.notes = {}

    def add_edges(self, the_edge_list):
        self._edges.extend(the_edge_list)

    def remove_edges(self, the_edge_list):
        for the_edge in the_edge_list:
            self._edges.remove(the_edge)

    def set_nodetype(self, the_new_nodetype):
        setattr(self, '_nodetype', the_new_nodetype)


