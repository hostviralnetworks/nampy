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

    def get_connected_nodes(self):
        the_edge_list = self.get_edges()
        the_connected_nodes = set([])
        for the_edge in the_edge_list:
            the_node_pair = the_edge.get_node_pair()
            the_index = the_node_pair.index(self)
            if the_index == 1:
                the_connected_nodes.add(the_node_pair[0])
            else:
                the_connected_nodes.add(the_node_pair[1])
        return the_connected_nodes
            
    def get_network(self):
        return self._network       

    def copy(self):
        # Copy all of the information
        # in the node, but we won't copy over
        # the associations such as the
        # edges and network.
        the_copy = Node(self.id)
        # Nodetype is a string rather than
        # real association so we copy it over
        the_copy.set_nodetype(self.get_nodetype())
        the_copy.source = self.source
        the_copy.notes = deepcopy(self.notes)
        return the_copy
