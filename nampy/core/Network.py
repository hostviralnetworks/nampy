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
        self.add_nodetype('monopartite')

        
    def add_nodetype(self, the_node_type_id):
        """ Add a node type to the model.
        
        Arguments:
         the_node_type_id: a string indicating the node type
         NOTE: MUST BE UNIQUE

        """
        if len(self.nodetypes) > 0:
            if the_node_type_id in [x.id for x in self.nodetypes]:
                return self.nodetypes.get_by_id(the_node_type_id)
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
        return self.matrix
        


    def update(self):
        self.nodetypes._generate_index()
        for the_nodetype in self.nodetypes:
            the_nodetype.nodes._generate_index()
        for the_edge in self.edges:
            node_pair_list = the_edge._nodes
            the_edge.id = node_pair_list[0].id + '___' + node_pair_list[1].id
        self.edges._generate_index()
        all_node_ids = []
        for the_nodetype in self.nodetypes:
            for the_node in the_nodetype.nodes:
                if len(the_node._edges) == 0:
                    print "Warning, orphan node: %s" %(the_node.id)
                if the_node.id in all_node_ids:
                    print "Warning, repeated node ids. This will break things, please fix %s." %(the_node.id)
                all_node_ids.append(the_node.id)
            the_nodetype.nodes._generate_index()
        self.update_matrix()
        # TODO: add a check for multipartite networks so edges don't
        # connect nodes of the same nodetype
        # no, this is not right, e.g. want to allow
        # host-host and host-virus interactions


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
                # might want to enforce this in the future
                # if 'monopartite' in nodetype_to_node_dict.keys():
                #     print "Error, cannot convert to multipartite if all nodes have monopartite nodetype designation."
                # else:
                old_nodetype = self.nodetypes.get_by_id("monopartite")
                for the_nodetype_id in nodetype_to_node_dict.keys():
                    the_nodetype = self.add_nodetype(the_nodetype_id)
                    the_nodes = nodetype_to_node_dict[the_nodetype_id]
                    old_nodetype.nodes.remove_subset(the_nodes)
                    the_nodetype.add_nodes(the_nodes)
                if len(old_nodetype.nodes) == 0:
                    setattr(old_nodetype, '_network', None)
                    self.nodetypes.remove(old_nodetype)  
            else:
                print "Need more than one nodetype designation to convert to multipartite."


    def convert_to_monopartite(self):
        #
        if (len(self.nodetypes) == 1) & (self.nodetypes[0].id == 'monopartite'):
            print "You already have a monopartite network. Exiting..."
        else:
            the_old_nodetype_ids = [x.id for x in self.nodetypes]
            if 'monopartite' not in [x.id for x in self.nodetypes]:
                the_nodetype = self.add_nodetype("monopartite")
            else:
                the_nodetype = self.nodetypes.get_by_id("monopartite")
            the_node_list = []
            for the_old_nodetype_id in the_old_nodetype_ids:
                if the_old_nodetype_id != "monopartite":
                    the_old_nodetype = self.nodetypes.get_by_id(the_old_nodetype_id)
                    the_node_list += the_old_nodetype.nodes
                    self.nodetypes.remove(the_old_nodetype)
                    setattr(the_old_nodetype, '_network', None)
            the_nodetype.nodes.extend(the_node_list)
                    

    def connect_node_pair(self, the_node_pair, the_weight = 1):
        """ Connect two nodes and creates an edge.

        Arguments:
         the_node_list: a list of 2 node objects
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


    def connect_node_pair_set(self, the_node_pair_list, **kwargs):
        """ Connect a set of node pairs to create edges.

        Arguments:
         the_node_pair_list: a list of [node1, node2]
          note NODE IDs MUST BE UNIQUE

        kwargs:
         the_weight_list:  alist of values to use for weights

        Returns:
         the_edge_list
        
        """
        # Could add more proofing of the input
        # but this needs to be fast.
        # We delay adding object references
        # because regenerating the indexes is slow
        from collections import defaultdict
        node_to_new_edge_dict = defaultdict(list)
        
        if 'the_weight_list' in kwargs:
            if len(kwargs['the_weight_list']) == len(the_node_pair_list):
                the_weight_list = kwargs['the_weight_list']
            else:
                the_weight_list = [1 for x in the_node_pair_list]
        else:
            the_weight_list = [1 for x in the_node_pair_list]

        the_edge_list = []
        for i, the_node_pair in enumerate(the_node_pair_list):
            the_weight = the_weight_list[i]
            the_edge = Edge(the_node_pair, the_weight)
            setattr(the_edge, '_network', self)
            the_edge_list.append(the_edge)
            for the_node in the_node_pair:
                node_to_new_edge_dict[the_node].append(the_edge)

        for the_node in node_to_new_edge_dict.keys():
            the_node.add_edges(node_to_new_edge_dict[the_node])

        self.edges.extend(the_edge_list)
        
        return the_edge_list
            

    def remove_edges(self, the_edge_list):
        if sum([type(the_edge) == str for the_edge in the_edge_list]) == len(the_edge_list):
            the_edge_list = [self.edges.get_by_id(the_edge) for the_edge in the_edge_list]
        else:
            # refilter these to make sure references are OK
            the_edge_list = [self.edges.get_by_id(x.id) for x in the_edge_list]

        from collections import defaultdict
        node_edge_remove_dict = defaultdict(set)

        for the_edge in the_edge_list:
            setattr(the_edge, '_network', None)
            the_node_set = the_edge.get_node_pair()
            for the_node in the_node_set:
                node_edge_remove_dict[the_node].add(the_edge)

        self.edges.remove_subset(the_edge_list)
        for the_node in node_edge_remove_dict.keys():
            the_node.remove_edges(node_edge_remove_dict[the_node])
            

    def get_node_locations(self):
        if (len(self.nodetypes) == 1) & (self.nodetypes[0] == 'monopatite'):
            the_node_locations = {x: 'monopartite' for x in self.nodetypes[0].nodes}
        else:
            the_node_locations = {}
            for the_nodetype in self.nodetypes:
                the_node_locations.update({x.id: the_nodetype.id for x in the_nodetype.nodes})
        return the_node_locations
