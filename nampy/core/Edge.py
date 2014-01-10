# Note hyperedges aren't supported yet.
from .Object import Object
from .Node import Node
from .DictList import DictList
from .parameters import edge_id_separator
    
class Edge(Object):
    """

    """
    
    def __init__(self, node_pair_list, the_weight = 1):
        """

        """
        if ((type(node_pair_list[0]) == Node) & (type(node_pair_list[1]) == Node)):
            id = node_pair_list[0].id + edge_id_separator + node_pair_list[1].id
            Object.__init__(self, id)
            self._nodes = DictList(node_pair_list)
            self.weight = the_weight
            self._network = None


    def update_id(self):
        node_pair_list = self._nodes
        self.id = node_pair_list[0].id + edge_id_separator + node_pair_list[1].id
        self._network.edges._generate_index()

    def get_node_pair(self):
        return self._nodes

    # This method might be a useful addition in the future?
    # For now, try to keep more detailed multi-class
    # operations in the manipulation module or as
    # components of the Network class
    # def replace_node(self, the_old_node, the_new_node):
    #    the_node_pair = self.get_node_pair()
    #    the_index = the_node_pair.index(the_old_node)
    #    the_edge._nodes[the_index] = the_new_node
    #    self.update_id()
