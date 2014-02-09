
def write_network_textfile(the_network, **kwargs):
    """ Write a simple tab-delimited textfile that can serve as
    table to import the network to cytosape

    Arguments:
     the_network: a nampy network object.  Note node_id_1
     and node_id_2 are extracted from the edges and
     can be used to define the "interaction" in Cytoscape.

    kwargs:
     properties_dict: a dicts of additional
      edge properties to write, with the
      property as the top level key.  These will 
      have each key corresponding to an edge ID.
     exclude_nodetypes: a list of nodetypes to avoid 
      including in the output.  Edges connecting
      nodes to these nodetypes will be ignored.
     

    """
    from .networkio import write_dict_to_textfile
    from ..core import NodeType
    continue_flag = True
    
    if 'properties_dict' in kwargs:
        properties_dict = kwargs['properties_dict']
    else:
        properties_dict = {}

    exclude_nodetypes = []
    if 'exclude_nodetypes' in kwargs:
        test_nodetypes = kwargs['exclude_nodetypes']
        for the_nodetype in test_nodetypes:
            if type(the_nodetype) == NodeType:
                exclude_nodetypes.append(the_nodetype.id)  
            else:
                exclude_nodetypes.append(the_nodetype) 
    
    the_output_dict = {}

    if continue_flag:
        the_edge_ids = []
        for the_edge in the_network.edges:
            the_node_pair = the_edge.get_node_pair()
            the_node_1 = the_node_pair[0]
            the_node_2 = the_node_pair[1]
            if ((the_node_1.get_nodetype() not in exclude_nodetypes) & (the_node_2.get_nodetype() not in exclude_nodetypes)):
                the_output_dict[the_edge.id] = {}
                the_output_dict[the_edge.id]['node_1_id'] = the_node_pair[0].id
                the_output_dict[the_edge.id]['node_2_id'] = the_node_pair[1].id
                the_output_dict[the_edge.id]['weight'] = the_edge.weight
                the_edge_ids.append(the_edge.id)

        # Try to write notes to file as individual
        # entries rather than lists when possible.
        list_note_keys = set([])
        one_element_note_keys = set([])        
        for the_edge_id in the_edge_ids:
            the_edge = the_network.edges.get_by_id(the_edge_id)
            for the_key in the_edge.notes.keys():
                if len(the_edge.notes[the_key]) == 1:
                    if the_key not in list_note_keys:
                        one_element_note_keys.add(the_key)
                else:
                    if the_key in one_element_note_keys:
                        one_element_note_keys.pop(the_key)
                    list_note_keys.add(the_key)
        
        for the_edge_id in the_edge_ids:
            the_edge = the_network.edges.get_by_id(the_edge_id)
            for the_key in the_edge.notes.keys():
                if the_key in list_note_keys:
                    the_output_dict[the_edge.id][the_key] = the_edge.notes[the_key]
                elif the_key in one_element_note_keys:
                    the_output_dict[the_edge.id][the_key] = the_edge.notes[the_key][0]

        for the_property in properties_dict.keys():
            for the_id in properties_dict[the_property].keys():
                if the_id in the_output_dict.keys():
                    if the_id in the_edge_ids:
                        the_output_dict[the_id][the_property] = properties_dict[the_property][the_id]

        write_dict_to_textfile(the_network.id + '_network_table.txt', the_output_dict, 'model_edge_id')
            
        
def write_node_attributes_to_textfile(the_network, **kwargs):
    """ Write a simple tab-delimited textfile that can serve as
    table to import the network to cytosape

    Arguments:
     the_network: a nampy network object.

    kwargs:
     properties_dict: a dicts of 
      additional node properties to write, with the
      property as the top level key.  These will have
      each key corresponding to a node ID.
     exclude_nodetypes: a list of nodetypes to avoid 
      including in the output.
     

    """
    from .networkio import write_dict_to_textfile
    from ..core import NodeType
    continue_flag = True

    
    if 'properties_dict' in kwargs:
        properties_dict = kwargs['properties_dict']
    else:
        properties_dict = {}

    exclude_nodetypes = []
    if 'exclude_nodetypes' in kwargs:
        test_nodetypes = kwargs['exclude_nodetypes']
        for the_nodetype in test_nodetypes:
            if type(the_nodetype) == NodeType:
                exclude_nodetypes.append(the_nodetype.id)
            else:
                exclude_nodetypes.append(the_nodetype) 
            
    the_output_dict = {}

    if continue_flag:
        the_node_ids = []
        for the_nodetype in the_network.nodetypes:
            if the_nodetype.id not in exclude_nodetypes:
                for the_node in the_nodetype.nodes:
                    checked_nodetype = the_node.get_nodetype()
                    if checked_nodetype not in exclude_nodetypes:
                        the_output_dict[the_node.id] = {}
                        the_output_dict[the_node.id]['nodetype'] = the_node.get_nodetype()
                        the_output_dict[the_node.id]['source'] = the_node.source
                        the_node_ids.append(the_node.id)

        # Try to write notes to file as individual
        # entries rather than lists when possible.
        for the_nodetype in the_network.nodetypes:
            list_note_keys = set([])
            one_element_note_keys = set([])
            
            for the_node in the_nodetype.nodes:
                if the_node.id in the_node_ids:
                    for the_key in the_node.notes.keys():
                        if len(the_node.notes[the_key]) == 1:
                            if the_key not in list_note_keys:
                                one_element_note_keys.add(the_key)
                        else:
                            if the_key in one_element_note_keys:
                                one_element_note_keys.pop(the_key)
                            list_note_keys.add(the_key)
                            
            for the_node in the_nodetype.nodes:
                if the_node.id in the_node_ids:
                    for the_key in the_node.notes.keys():
                        if the_key in list_note_keys:
                            the_output_dict[the_node.id][the_key] = the_node.notes[the_key]
                        elif the_key in one_element_note_keys:
                            the_output_dict[the_node.id][the_key] = the_node.notes[the_key][0]

        for the_property in properties_dict.keys():
            for the_id in properties_dict[the_property].keys():
                if the_id in the_output_dict.keys():
                    the_output_dict[the_id][the_property] = properties_dict[the_property][the_id]

        write_dict_to_textfile(the_network.id + '_node_attribute_table.txt', the_output_dict, 'node_id')
