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

    # copy() methods would not be valid for edges.
    # No point in duplicating an existing edges,
    # though you might want to link nodes
    # in a different network in a similar way.
