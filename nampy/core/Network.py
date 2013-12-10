from numpy import zeros
from scipy.sparse import dok_matrix
from copy import deepcopy
from .Object import Object
from .Edge import Edge
from .NodeType import NodeType
from .DictList import DictList
#

class Network(Object):
    """

    """

    def __init__(self, id):
        """

        """
        Object.__init__(self, id)
        self.nodetypes = DictList()
        self.edges = DictList()
        self.update_matrix()
        # monopartite is a reserved
        # string to indicate we won't enforce
        # that distinct node types be broken into
        # separate lists
        self.add_node_type('monopartite')

        
    def add_node_type(self, the_node_type_id):
        """ Add a node type to the model.
        
        Arguments:
         the_node_type_id: a string indicating the node type
         NOTE: MUST BE UNIQUE

        """
        if len(self.nodetypes) > 0:
            if the_node_type_id in [x.id for x in self.nodetypes]:
                return self.nodetypes.get_by_id(the_node_type)
        the_node_type = NodeType(the_node_type_id)
        setattr(the_node_type, '_network', self)
        self.nodetypes.append(the_node_type)
        return the_node_type
            


    def update_matrix(self):
        # TODO: add a check for bipartite matrices here
        
        if len([x.id for x in self.nodetypes]) == 1:
            if self.nodetypes[0].id == 'monopartite':
                the_node_dictlist = self.nodetypes[0].nodes
                if len(the_node_dictlist) > 0:
                    self.matrix = dok_matrix((len(the_node_dictlist),len(the_node_dictlist)), dtype = float)
                    for the_node in the_node_dictlist:
                        for the_edge in the_node._edges:
                            i = self.nodetypes[0].nodes.index(the_edge._nodes[0])
                            j = self.nodetypes[0].nodes.index(the_edge._nodes[1])
                            self.matrix[i, j] = the_edge.weight
                            self.matrix[j, i] = the_edge.weight
                    return self.matrix
        self.matrix = zeros((0,0), dtype = float)
        


    def update(self):
        self.update_matrix()
        all_node_ids = []
        for the_nodetype in self.nodetypes:
            for the_node in the_nodetype.nodes:
                if len(the_node._edges) == 0:
                    print "Warning, orphan node: %s" %(the_node.id)
                if the_node.id in all_node_ids:
                    print "Warning, repeated node ids. This will break things, please fix %s." %(the_node.id)
                all_node_ids.append(the_node.id)
        # TODO: add a check for multipartite networks so edges don't
        # connect nodes of the same nodetype


    def convert_to_multipartite(self):
        #
        from collections import OrderedDict

        if (len(self.nodetypes) == 1) & (self.nodetypes[0].id == 'monopartite'):
            nodetype_to_node_dict = OrderedDict()
            for the_nodetype in self.nodetypes:
                for the_node in the_nodetype.nodes:
                    the_nodetype_id = the_node._nodetype
                    if the_nodetype_id not in nodetype_to_node_dict.keys():
                        nodetype_to_node_dict[the_nodetype_id] = []
                    nodetype_to_node_dict[the_nodetype_id].append(the_node)
                    # Make sure all node_ids are unique to avoid problems here
                    # Warns in update_self
                
            existing_nodetype_ids = [x.id for x in self.nodetypes]

            if len(nodetype_to_node_dict.keys()) > 1:
                if 'monopartite' in nodetype_dict.keys():
                    print "Error, cannot convert to multipartite if nodes have monopartite nodetype designation."
                else:
                    for the_nodetype_id in nodetype_to_node_dict.keys():
                        the_nodetype = self.add_nodetype(the_nodetype_id)
                        the_nodes = nodetype_to_node_dict[the_nodetype_id]
                        the_nodetype.extend(the_nodes)
                old_nodetype = self.nodetype.get_by_id("monopartite")
                self.nodetype.pop(old_nodetype)  
            else:
                print "Need more than one nodetype designation to convert to multipartite."


    def convert_to_monopartite(self):
        #
        if (len(self.nodetypes) == 1) & (self.nodetypes[0].id == 'monopartite'):
            print "You already have a monopartite network. Exiting..."
        else:
            if 'monopartite' not in [x.id for x in nodetypes]:
                the_node_list = []
                for the_nodetype in self.nodetypes:
                    the_node_list += the_nodetype.nodes
                    self.pop(the_nodetype)
                the_nodetype = self.add_node_type("monopartite")
                the_nodetype.nodes.extend(the_node_list)
                    

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
            setattr(the_edge, '_network', self)
        return the_edge
