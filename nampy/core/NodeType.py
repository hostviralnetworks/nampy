from copy import deepcopy
import re
from .Object import Object
from .DictList import DictList
from .Node import Node

class NodeType(Object):
    """

    """
    def __init__(self, id, **kwargs):
        # Will allow for a 'type' argument
        # to facilitate making multipartite
        # graphs
        Object.__init__(self, id, **kwargs)
        # Object includes notes and annotation fields
        self.nodes = DictList()
        self._network = None


    def add_nodes(self, the_node_list):
        """ Add nodes to the NodeType.
        
        Arguments:
         the_node_list: a list of nodes id's (strings) or nodes
          note NODE IDs MUST BE UNIQUE

        """
        # TODO: add a check for other nodetypes not using a node with the same id
        if sum([type(x) == str for x in the_node_list]) == len(the_node_list):
            the_node_list_by_id = the_node_list
        else:
            the_node_list_by_id = [x.id for x in self.nodes]
        existing_nodes_by_id = [x.id for x in self.nodes]
        # Safety check: filter out to make sure redundant nodes are not added
        the_node_list_by_id = [x for x in the_node_list_by_id if x not in existing_nodes_by_id]
        the_node_list = [Node(the_id) for the_id in the_node_list_by_id]
        [setattr(x, '_network', self._network) for x in the_node_list]
        if self.id != 'monopartite':
            [setattr(x, '_nodetype', self.id) for x in the_node_list]
        self.nodes.extend(the_node_list)


    def remove_nodes(self, the_node_list):
        if sum([type(the_node) == str for the_node in the_node_list]) == len(the_node_list):
            the_node_list = [self.nodes.get_by_id(the_node) for the_node in the_node_list]
        edges_to_remove = set([])
        for the_node in the_node_list:
            for the_edge in the_node._edges:
                edges_to_remove.add(the_edge)
            setattr(the_node, '_network', None)
            setattr(the_node, '_nodetype', None)
            self.nodes.remove(the_node)
        self._network.remove_edges(edges_to_remove)
            

