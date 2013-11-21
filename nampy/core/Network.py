from numpy import zeros
from scipy.sparse import dok_matrix
from copy import deepcopy
from .Object import Object
from .Edge import Edge
from .Node import Node
from .DictList import DictList
#


class Network(Object):
    """

    """

    def __init__(self, id):
        """

        """
        Object.__init__(self, id)
        self.nodes = DictList()
        self.edges = DictList()
        self.update_matrix()


    def add_nodes(self, the_node_list):
        """ Add notes to the model.
        
        Arguments:
         the_node_list: a list of nodes id's (strings) or nodes
          note NODE IDs MUST BE UNIQUE

        """
        if sum([type(x) == str for x in the_node_list]) == len(the_node_list):
            the_node_list_by_id = the_node_list
        else:
            the_node_list_by_id = [x.id for x in self.nodes]
        existing_nodes_by_id = [x.id for x in self.nodes]
        # Safety check: filter out to make sure redundant nodes are not added
        the_node_list_by_id = [x for x in the_node_list_by_id if x not in existing_nodes_by_id]
        the_node_list = [Node(the_id) for the_id in the_node_list_by_id]
        [setattr(x, '_model', self) for x in the_node_list]
        self.nodes.extend(the_node_list)


    def connect_node_pair(self, the_node_pair, the_weight = 1):
        """ Connect two nodes and create an edge.

        Arguments:
         the_node_list: a list of nodes id's (strings) or nodes
          note NODE IDs MUST BE UNIQUE

        Returns:
         the_edge
        
        """
        # Could add more proofing of the input
        # but this needs to be fast.
        the_edge = Edge(the_node_pair, the_weight)
        if type(the_edge) == Edge:
            self.edges.append(the_edge)
            for the_node in the_node_pair:
                the_node.add_edges([the_edge])
            setattr(the_edge, '_model', self)
        return the_edge
            


    def remove_edges(self, the_edge_list):
        if sum([type(the_edge) == str for the_edge in the_edge_list]) == len(the_edge_list):
            the_edge_list = [self.edges.get_by_id(the_edge) for the_edge in the_edge_list]
        for the_edge in the_edge_list:
            the_node_pair = the_edge._nodes
            the_node_pair[0].remove_edges([the_edge])
            the_node_pair[1].remove_edges([the_edge])
            setattr(the_edge, '_model', None)
            #the_index = self.edges.index(the_edge)
            #self.edges.pop(the_index)
            self.edges.remove(the_edge) 


    def remove_nodes(self, the_node_list):
        if sum([type(the_node) == str for the_node in the_node_list]) == len(the_node_list):
            the_node_list = [self.nodes.get_by_id(the_node) for the_node in the_node_list]
        edges_to_remove = []
        for the_node in the_node_list:
            for the_edge in the_node._edges:
                edges_to_remove.append(the_node)
            self.remove_edges(edges_to_remove)
            #the_index = self.nodes.index(the_node)
            self.nodes.remove(the_node)
            setattr(the_node, '_model', None)


    def update_matrix(self):
        if len(self.nodes) > 0:
            self.matrix = dok_matrix((len(self.nodes),len(self.nodes)), dtype = float)
            #self.matrix = dok_matrix((len(self.nodes),len(self.nodes)), dtype=float32)
            for the_node in self.nodes:
                for the_edge in the_node._edges:
                    i = self.nodes.index(the_edge._nodes[0])
                    j = self.nodes.index(the_edge._nodes[1])
                    self.matrix[i, j] = the_edge.weight
                    self.matrix[j, i] = the_edge.weight
        else:
            self.matrix = zeros((0,0), dtype = float)


    def update(self):
        self.update_matrix()
        for the_node in self.nodes:
            if len(the_node._edges) == 0:
                print "Warning, orphan node: %s" %(the_node.id)

