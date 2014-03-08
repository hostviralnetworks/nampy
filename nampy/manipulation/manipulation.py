from ..core.shared_functions import test_kwarg
from ..core.parameters import edge_id_separator
    
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


def transfer_edges(the_new_node, the_old_node, **kwargs):
    """ Move edge associations to another node.

    Arguments:
     the_new_node: a node object to be copied
     the_edges: the list of edges

    kwargs:
     verbose: inform if there are problems during edge transfer
     duplicate_target_id_behavior: if a duplicate target edge id is found in
                                   the target node 'delete' or 'keep' the old
                                   edge. 

    """
    
    from copy import deepcopy
    from ..core.Node import Node
    from ..core.Edge import Edge

    verbose = test_kwarg('verbose', kwargs, [True, False])
    # By default we will remove all old edge associations
    duplicate_target_id_behavior = test_kwarg('duplicate_target_id_behavior', kwargs, ['delete', 'keep'])
    
    the_edges = the_old_node.get_edges()
        
    the_edges_to_remove_from_old_node = []

    #    if the_new_node.id == '3107':
    #    import pdb
    #    pdb.set_trace()
    
    for the_edge in the_edges:
        # Could make these steps
        # into an Edge method, but this does
        # not make sense to stand alone
        the_node_pair = the_edge.get_node_pair()
        the_index = the_node_pair.index(the_old_node)
        the_edge._nodes[the_index] = the_new_node
        the_edge.update_id()
        # Now we add this edge to the new node
        # First we make sure an edge with the same id
        # doesn't already exist
        the_node_pair = the_edge.get_node_pair()
        the_first_id = the_node_pair[0].id + edge_id_separator + the_node_pair[1].id
        the_second_id = the_node_pair[1].id + edge_id_separator + the_node_pair[0].id
        the_new_node_edges = the_new_node.get_edges()
        the_new_node_edges_ids = [x.id for x in the_new_node_edges]
        if ((the_first_id not in the_new_node_edges_ids) & (the_second_id not in the_new_node_edges_ids)):
            the_new_node.add_edges([the_edge])
            # Also update network references
            setattr(the_edge, '_network', the_new_node.get_network())
            if the_edge not in the_edge._network.edges:
                the_edge._network.edges.append(the_edge)
            the_edges_to_remove_from_old_node.append(the_edge)
        else:
            # Then this could be a duplicate id that we 
            # should delete, but before deleting make sure it is not
            # the same edge already present in the_new_node
            if the_edge not in the_new_node_edges:
                if the_first_id in the_new_node_edges_ids:
                    dupe_id = the_first_id
                elif the_second_id in the_new_node_edges_ids:
                    dupe_id = the_second_id  
                if verbose:
                    print 'Warning, edge transfer failed. Edge id %s already exists for target node.' % dupe_id
                # If this is the case, perform the selected fail behavior
                if duplicate_target_id_behavior == 'delete':
                    setattr(the_edge, '_network', None)
                    # Remove the duplicate from the network
                    if the_edge in the_old_node._network.edges:
                        the_old_node._network.edges.remove(the_edge)
                    the_edges_to_remove_from_old_node.append(the_edge)
                else:
                    # 'keep' is the other option,
                    # we should revert the edge
                    # in this case
                    the_node_pair = the_edge.get_node_pair()
                    # This should still be in memory
                    # the_index = the_node_pair.index(the_old_node)
                    the_edge._nodes[the_index] = the_old_node
                    the_edge.update_id()
            # If the edge is already in the new node, then
            # go ahead and perform the default. 
            else:
                if verbose:
                    print 'Warning, edge transfer failed. Edge id %s already exists for target node.' % the_edge.id
                # If this is the case, perform the selected fail behavior
                if duplicate_target_id_behavior == 'delete':
                    setattr(the_edge, '_network', None)
                    # Remove the duplicate from the network
                    if the_edge in the_old_node._network.edges:
                        the_old_node._network.edges.remove(the_edge)
                    the_edges_to_remove_from_old_node.append(the_edge)
                else:
                    # 'keep' is the other option,
                    # we should revert the edge
                    # in this case
                    the_node_pair = the_edge.get_node_pair()
                    # This should still be in memory
                    # the_index = the_node_pair.index(the_old_node)
                    the_edge._nodes[the_index] = the_old_node
                    the_edge.update_id()                
        
            
    # remove the edges from the old node
    the_old_node.remove_edges(the_edges_to_remove_from_old_node)


def merge_nodes(the_node_1, the_node_2):
    """ Merge two nodes and remove the second
    from the network

    Arguments:
     the_first_node: 'dominant' node
     the_second_node: node to be eliminated
    

    """
    continue_flag = True

    if the_node_1.get_network() != the_node_1.get_network():
        confinue_flag = False
        print 'Both nodes must belong to the same network, exiting...'

    if continue_flag:
        transfer_edges(the_node_1, the_node_2)
        the_node_1.notes = merge_notes(the_node_1, the_node_2)
        the_network = the_node_2.get_network()
        the_node_location_dict = the_network.get_node_locations()
        the_nodetype_id = the_node_location_dict[the_node_2.id]
        the_network.nodetypes.get_by_id(the_nodetype_id).remove_nodes([the_node_2])

    return the_node_1


def merge_notes(the_object_1, the_object_2):
    """
    """
    from copy import deepcopy
    continue_flag = True
    merge_notes = {}
    if ('notes' not in dir(the_object_1)) | ('notes' not in dir(the_object_2)):
        print 'Both objects must have a notes attribute, exiting...'
        continue_flag = False
    if continue_flag:
        notes_1 = the_object_1.notes
        notes_2 = the_object_2.notes
        all_keys = list(set(notes_1.keys() + notes_2.keys()))
        for the_key in all_keys:
            combined_data = []
            if the_key in notes_1.keys():
                if type(notes_1[the_key]) == list:
                    combined_data += deepcopy(notes_1[the_key])
                else:
                    combined_data += [deepcopy(notes_1[the_key])]
            if the_key in notes_2.keys():
                if type(notes_2[the_key]) == list:
                    combined_data += deepcopy(notes_2[the_key])
                else:
                    combined_data += [deepcopy(notes_2[the_key])]
            merge_notes[the_key] = list(set(combined_data))

    return merge_notes


def merge_networks_by_node(the_first_network, the_second_network, new_network_id, **kwargs):
    """ Merges two networks based on node id's.

    Arguments:
     the_first_network
     the_second_network
     new_network_id

    kwargs:
     verbose

    """
    from ..core.Node import Node
    from ..core.Edge import Edge
    from ..core.Network import Network
    from ..core.NodeType import NodeType
    from copy import deepcopy
        
    verbose = test_kwarg('verbose', kwargs, [False, True])
    
    the_network = Network(new_network_id)
    the_network_1_node_locations = the_first_network.get_node_locations()
    the_network_2_node_locations = the_second_network.get_node_locations()
    the_first_network_nodetype_ids = list(set(the_network_1_node_locations.values()))
    the_first_network_nodetype_locations = {}
    for the_nodetype_id in the_first_network_nodetype_ids:
        the_first_network_nodetype_locations[the_nodetype_id] = []
    for the_node_id in the_network_1_node_locations.keys():
        the_first_network_nodetype_locations[the_network_1_node_locations[the_node_id]].append(the_node_id)
    the_second_network_nodetype_ids = list(set(the_network_2_node_locations.values()))
    the_second_network_nodetype_locations = {}
    for the_nodetype_id in the_second_network_nodetype_ids:
        the_second_network_nodetype_locations[the_nodetype_id] = []
    for the_node_id in the_network_2_node_locations.keys():
        the_second_network_nodetype_locations[the_network_2_node_locations[the_node_id]].append(the_node_id)

    if verbose:
        print "Creating the nodes..."  

    the_nodes_to_add = []
    the_node_ids_to_add = []
    for the_nodetype_id in the_first_network_nodetype_locations.keys():
        the_nodetype = the_first_network.nodetypes.get_by_id(the_nodetype_id)
        for the_node_id in the_first_network_nodetype_locations[the_nodetype_id]:
            the_nodes_to_add.append(the_nodetype.nodes.get_by_id(the_node_id).copy())
            the_node_ids_to_add.append(the_node_id)
    for the_nodetype_id in the_second_network_nodetype_locations.keys():
        the_nodetype = the_second_network.nodetypes.get_by_id(the_nodetype_id)
        for the_node_id in the_second_network_nodetype_locations[the_nodetype_id]:
            if the_node_id not in the_node_ids_to_add:
                the_nodes_to_add.append(the_nodetype.nodes.get_by_id(the_node_id).copy())
                the_node_ids_to_add.append(the_node_id)
                # Otherwise we let information in the existing node override
 
    the_nodetype = the_network.nodetypes[0]
    the_nodetype.add_nodes(the_nodes_to_add)

    if verbose:
        print "     ... the nodes are created."
        print "Linking the nodes, this may take a while..."
        
    # here we prioritize edges in the first network for the purpose of weights and notes
    the_edges_to_add = [x for x in the_first_network.edges]
    the_first_id = [x._nodes[0].id + edge_id_separator + x._nodes[1].id for x in the_first_network.edges]
    the_second_id = [x._nodes[1].id + edge_id_separator + x._nodes[0].id for x in the_first_network.edges]
    for the_edge in the_second_network.edges:
        if the_edge.id not in the_first_id:
            if the_edge.id not in the_second_id:
                the_edges_to_add.append(the_edge)

    the_first_nodes_to_link = [the_nodetype.nodes.get_by_id(x._nodes[0].id) for x in the_edges_to_add]
    the_second_nodes_to_link = [the_nodetype.nodes.get_by_id(x._nodes[1].id) for x in the_edges_to_add]
    the_node_pair_list = [[the_first_nodes_to_link[i], the_second_nodes_to_link[i]] for i in range(0, len(the_edges_to_add))]
    the_weights = [x.weight for x in the_edges_to_add]
    the_notes = [deepcopy(x.notes) for x in the_edges_to_add]
    
    the_edges = the_network.connect_node_pair_list(the_node_pair_list, the_weight_list = the_weights)
    
    for i, the_edge in enumerate(the_edges):
        the_edge.notes = the_notes[i]
        
    if verbose:
        print "Completed testing %s node pairs" %(str(len(the_edges)))      

    return the_network


def make_subnetwork(the_network, the_node_id_list, the_id):
    """ Make an independent network that
    is a subset of network components consisting of the
    node ids in the node_id_list and any edges
    in the network that happen to connect them.

    Arguments:
     the_network: a NAMpy network object
     the_node_id_list: a list of node ids (not node objects)
     the_id: an id string for the new "sub" network.

    Returns:
     the_subnetwork: a new network object that has the
      indicated nodes and selected edges

    """
    from ..core.Network import Network
    from ..core.Node import Node

    # A quick fix up to be nice for those that insist on using Node objects
    nodes_to_pop = []
    ids_to_add = []
    for the_index, the_node_id in enumerate(the_node_id_list):
        if type(the_node_id) == Node:
            nodes_to_pop.append(the_node_id)
            ids_to_add.append(the_node_id.id)
    for the_node in nodes_to_pop:
        the_index = the_node_id_list.index(the_node)
        the_node_id_list.pop(the_index)
    the_node_id_list += ids_to_add
            
    the_subnetwork = Network(the_id)
    the_node_locations = the_network.get_node_locations()
    the_node_id_list = [x for x in the_node_id_list if x in the_node_locations.keys()]
    node_ids_to_remove = [x for x in the_node_locations.keys() if x not in the_node_id_list]
    for the_node_id in node_ids_to_remove:
        the_node_locations.pop(the_node_id)
    the_nodetype_ids = list(set(the_node_locations.values()))
    the_nodetype_locations = {}
    for the_nodetype_id in the_nodetype_ids:
        the_nodetype_locations[the_nodetype_id] = []
    for the_node_id in the_node_locations.keys():
        the_nodetype_locations[the_node_locations[the_node_id]].append(the_node_id)
    the_nodes_to_add = []
    for the_nodetype_id in the_nodetype_ids:
        the_nodetype = the_network.nodetypes.get_by_id(the_nodetype_id)
        the_nodes_to_add += [the_node.copy() for the_node in the_nodetype.nodes if the_node.id in the_node_id_list]
    the_subnetwork.nodetypes[0].add_nodes(the_nodes_to_add)
    node_pairs_to_connect = []
    the_weight_list = []
    for the_edge in the_network.edges:
        edge_node_id_list = [x.id for x in the_edge.get_node_pair()]
        if (edge_node_id_list[0] in the_node_id_list) and (edge_node_id_list[1] in the_node_id_list): 
            node_pairs_to_connect.append([the_subnetwork.nodetypes[0].nodes.get_by_id(edge_node_id_list[0]), the_subnetwork.nodetypes[0].nodes.get_by_id(edge_node_id_list[1])])
            the_weight_list.append(the_edge.weight)
    the_subnetwork.connect_node_pair_list(node_pairs_to_connect, the_weight_list = the_weight_list)
    return the_subnetwork
        
