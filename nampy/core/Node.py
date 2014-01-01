from copy import deepcopy
import re
from .Object import Object
from .DictList import DictList
#from .Edge import Edge


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
        self._edges.remove_subset(the_edge_list)

    def set_nodetype(self, the_new_nodetype):
        setattr(self, '_nodetype', the_new_nodetype)

    def get_nodetype(self):
        return self._nodetype

    def get_edges(self):
        the_edge_list = [x for x in self._edges]
        return the_edge_list
            
    def get_network(self):
        return self._network       
