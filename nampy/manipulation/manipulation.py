
def add_source(the_network, source_dict, **kwargs):
    """ Script to add sources to a network for
    network analysis from a dictionary

    Arguments:
     the_network: a network object, modified in place
     source_dict: dict of nodes to add sources/sinks to
      must include 'value' entry

    kwargs:
     match_key_type: the key in the source_dict to use to make
      matches to the model.  If 'default', ids are used.

    Returns a tuple:
     the_network
     unmatched_ids_dict
                  
    """
    from copy import deepcopy

    continue_flag = True
    source_dict = deepcopy(source_dict)

    if 'match_key_type' in kwargs:
        match_key_type = kwargs['match_key_type'] 
    else:
        match_key_type = 'default'

    valid_match_key_types = ['default']
    for the_nodetype in the_network.nodetypes:
        for the_node in the_nodetype.nodes:
            valid_match_key_types += the_node.notes.keys()
    valid_match_key_types = list(set(valid_match_key_types))

    if match_key_type not in valid_match_key_types:
        print('Warning, %s not a valid key to match to for model nodes.  You must pick one from: %s' % (match_key_type, str(valid_match_key_types)))
        continue_flag = False            
        
    if continue_flag:
        dict_for_pairing = {}
        # might want to make this node type specific in the future
        for the_nodetype in the_network.nodetypes:
            for the_node in the_nodetype.nodes:
                if match_key_type == 'default':
                    pairing_key_list = [the_node.id]
                else:
                    pairing_key_list = the_node.notes[match_key_type]
                for pairing_key in pairing_key_list:
                    if pairing_key not in dict_for_pairing.keys():
                        dict_for_pairing[pairing_key] = {}
                        dict_for_pairing[pairing_key]['network_ids'] = []
                        dict_for_pairing[pairing_key]['source_dict_ids'] = []
                        dict_for_pairing[pairing_key]['network_ids'].append(the_node.id)
            for the_source_key in source_dict.keys():
                if match_key_type == 'default':
                    pairing_key_list = [the_source_key]
                else:
                    pairing_key_list = source_dict[the_source_key][match_key_type]
                for pairing_key in pairing_key_list:
                    if pairing_key not in dict_for_pairing.keys():
                        dict_for_pairing[pairing_key] = {}
                        dict_for_pairing[pairing_key]['network_ids'] = []
                        dict_for_pairing[pairing_key]['source_dict_ids'] = []
                    dict_for_pairing[pairing_key]['source_dict_ids'].append(the_source_key)

        unpaired_source_dict_keys = []
        pairing_key_model_one_to_source_dict_one = []
        unmatched_ids_dict = {}
        matched_ids_dict = {}
        unmatched_ids_dict['source_dict_one_to_network_many'] = {}
        unmatched_ids_dict['network_one_to_source_dict_many'] = {}
        unmatched_ids_dict['source_dict_many_to_network_many'] = {}
        unmatched_ids_dict['source_dict_none'] = {}
        unmatched_ids_dict['network_none'] = {}
        for the_pairing_key in dict_for_pairing.keys():
            the_mapping_type = ''
            network_ids = dict_for_pairing[the_pairing_key]['network_ids']
            source_dict_ids = dict_for_pairing[the_pairing_key]['source_dict_ids']
            if ((len(network_ids) == 1) & (len(source_dict_ids) == 1)):
                pairing_key_model_one_to_source_dict_one.append(the_pairing_key)
            elif ((len(network_ids) > 1) & (len(source_dict_ids) == 1)):
                the_mapping_type = 'source_dict_one_to_network_many'
                unpaired_source_dict_keys += source_dict_ids
            elif ((len(network_ids) == 1) & (len(source_dict_ids) > 1)):
                the_mapping_type = 'network_one_to_source_dict_many'
                unpaired_source_dict_keys += source_dict_ids
            elif ((len(network_ids) > 1) & (len(source_dict_ids) > 1)):
                the_mapping_type = 'network_many_to_source_dict_many'
                unpaired_source_dict_keys += source_dict_ids
            elif (len(network_ids) == 0):
                the_mapping_type = 'network_none'
                unpaired_source_dict_keys += source_dict_ids
            elif (len(source_dict_ids) == 0):
                the_mapping_type = 'source_dict_none'

            if len(the_mapping_type) >0:
                unmatched_ids_dict[the_mapping_type][the_pairing_key] = {}
                unmatched_ids_dict[the_mapping_type][the_pairing_key]['network_ids'] = network_ids
                unmatched_ids_dict[the_mapping_type][the_pairing_key]['source_dict_ids'] = source_dict_ids
            else:
                matched_ids_dict[the_pairing_key] = {}
                matched_ids_dict[the_pairing_key]['network_ids'] = network_ids
                matched_ids_dict[the_pairing_key]['source_dict_ids'] = source_dict_ids

        unpaired_source_dict_keys = list(set(unpaired_source_dict_keys))

        for the_nodettype in the_network.nodetypes:
            the_node_ids = [x.id for x in the_nodettype.nodes]
            for the_pairing_key in pairing_key_model_one_to_source_dict_one:
                network_id = dict_for_pairing[the_pairing_key]['network_ids'][0]
                source_dict_id = dict_for_pairing[the_pairing_key]['source_dict_ids'][0]
                if network_id in the_node_ids:
                    the_node = the_nodettype.nodes.get_by_id(network_id)
                    the_node.source = source_dict[source_dict_id]['value']
            
        print("Completed adding sources with %i pairings of %i sources in the source_dict and %i nodes in the model." % (len(pairing_key_model_one_to_source_dict_one), len(source_dict.keys()), len(the_network.nodetypes[0].nodes)))

        summary_dict = {}
        summary_dict['matched_ids'] = matched_ids_dict
        summary_dict['unmatched_ids'] = unmatched_ids_dict
            
        return the_network, summary_dict
    else:
        return None, None



def duplicate_node(the_node, new_id):
    """ Duplicates a node with associated data
    and adds it to the network as a duplicate
    node with a new id.

    Arguments:
     the_node: a node object to be copied
     new_id: string

    """
    
    from copy import deepcopy
    from ..core.Node import Node
    from ..core.Edge import Edge

    node_locations = the_node._network.get_node_locations()
    existing_ids = node_locations.keys()

    if new_id not in existing_ids:
        the_location = node_locations[the_node.id]
        the_new_node = Node(new_id)
        the_new_node._network = the_node._network
        the_new_node._nodetype = the_node._nodetype
        the_new_node.source = the_node.source
        the_new_node.notes = deepcopy(the_node.notes)
        the_edges_to_add = []
        for the_edge in the_node._edges:
            the_index = the_edge._nodes.index(the_node)
            if the_index == 0:
                the_nodes = [the_new_node, the_edge._nodes[1]]
            if the_index == 1:
                the_nodes = [the_edge._nodes[0], the_new_node]
            the_edge = Edge(the_nodes, the_edge.weight)
            the_edge._network = the_node._network
            the_edges_to_add.append(the_edge)
        the_new_node.add_edges(the_edges_to_add)
        the_nodetype_id = node_locations[the_node.id]
        the_node._network.nodetypes.get_by_id(the_nodetype_id).nodes.append(the_new_node)
        the_node._network.edges.extend(the_edges_to_add)
        # Do we really need to update here?
        the_node._network.update()
    return the_new_node


def transfer_edges(the_new_node, the_old_node):
    """ replace edge associations with a node.

    Arguments:
     the_new_node: a node object to be copied
     the_edges: a list of edges

    """
    
    from copy import deepcopy
    from ..core.Node import Node
    from ..core.Edge import Edge

    the_edges = the_old_node.get_edges()

    for the_edge in the_edges:
        the_index = the_edge._nodes.index(the_old_node)
        the_edge._nodes[the_index] = the_new_node
        the_edge.update_id()
        if the_edge.id not in [x.id for x in the_new_node._edges]:
            the_new_node._edges.append(the_edge)
            setattr(the_edge, '_network', the_new_node._network)
            if the_edge not in the_edge._network.edges:
                the_edge._network.edges.append(the_edge)
        else:
            # then this is a duplicate
            the_edge._network = None
            if the_edge in the_old_node._network.edges:
                the_old_node._network.edges.remove(the_edge)
                
            
        # there's a faster way to do this but leave it for now

    the_old_node.remove_edges(the_edges)

        
def merge_networks_by_node(the_first_network, the_second_network, new_network_id, **kwargs):
    """ Merges two networks based on node id's.

    Arguments:
     the_first_network
     the_second_network
     new_network_id

    kwargs:
    # merge_type: asymmetric or symmetric
     verbose

    """
    from ..core.Node import Node
    from ..core.Edge import Edge
    from ..core.Network import Network
    from ..core.NodeType import NodeType
    from copy import deepcopy

    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
    else:
        verbose = False

    the_network = Network(new_network_id)
    the_network_1_node_locations = the_first_network.get_node_locations()
    the_network_2_node_locations = the_second_network.get_node_locations()

    the_nodes_to_add = list(set(the_network_1_node_locations.keys() + the_network_2_node_locations.keys()))

    if verbose:
        print "Creating the nodes..."

    the_nodetype = the_network.nodetypes[0]
    the_nodes_to_add.sort()
    the_nodetype.add_nodes(the_nodes_to_add)

    for the_node in the_nodetype.nodes:
        if the_node.id in the_network_1_node_locations.keys():
            the_old_node = the_first_network.nodetypes.get_by_id(the_network_1_node_locations[the_node.id]).nodes.get_by_id(the_node.id)
            the_node._nodetype = the_old_node._nodetype
            the_node.source = the_old_node.source
            the_node.notes = deepcopy(the_old_node.notes)
            # if merge_type == 'symmetric':
        if the_node.id in the_network_2_node_locations.keys():
            # Give priority to data from the first network
            the_old_node = the_second_network.nodetypes.get_by_id(the_network_2_node_locations[the_node.id]).nodes.get_by_id(the_node.id)
            if the_node.id not in the_network_1_node_locations.keys():
                the_node._nodetype = the_old_node._nodetype
                the_node.source = the_old_node.source
                the_node.notes = deepcopy(the_old_node.notes)  
            else:
                for the_key in the_old_node.notes.keys():
                    if the_key not in the_node.notes:
                        the_node.notes[the_key] = deepcopy(the_old_node.notes[the_key])

    if verbose:
        print "     ... the nodes are created."
        print "Linking the nodes, this may take a while..."


    # here we prioritize edges in the first network for the purpose of weights and notes
    # is is assumed there is no redundancy in the edge ids
    the_edges_to_add = [x for x in the_first_network.edges]
    the_first_id = [x._nodes[0].id + '_' + x._nodes[1].id for x in the_first_network.edges]
    the_second_id = [x._nodes[1].id + '_' + x._nodes[0].id for x in the_first_network.edges]
    for the_edge in the_second_network.edges:
        if the_edge.id not in the_first_id:
            if the_edge.id not in the_second_id:
                the_edges_to_add.append(the_edge)

    the_first_nodes_to_link = [the_nodetype.nodes.get_by_id(x._nodes[0].id) for x in the_edges_to_add]
    the_second_nodes_to_link = [the_nodetype.nodes.get_by_id(x._nodes[1].id) for x in the_edges_to_add]
    the_node_pair_list = [[the_first_nodes_to_link[i], the_second_nodes_to_link[i]] for i in range(0, len(the_edges_to_add))]
    the_weights = [x.weight for x in the_edges_to_add]
    the_notes = [deepcopy(x.notes) for x in the_edges_to_add]
    
    the_edges = the_network.connect_node_pair_set(the_node_pair_list, the_weight_list = the_weights)
    
    for i, the_edge in enumerate(the_edges):
        the_edge.notes = the_notes[i]

    if verbose:
        print "Completed testing %s node pairs" %(str(len(the_edges)))      

    return the_network
