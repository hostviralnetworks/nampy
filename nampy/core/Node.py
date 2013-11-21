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
        self._network = None
        self.source = 0

    def add_edges(self, the_edge_list):
        self._edges.extend(the_edge_list)

    def remove_edges(self, the_edge_list):
        self._edges.extend(the_edge_list)
        #for the_edge in the_edge_list:
        #    the_index = self._edges.index(the_edge)
        #    self._edges.pop(the_index)
