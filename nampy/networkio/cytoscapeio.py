from ..core.shared_functions import test_kwarg
from ..core.parameters import edge_id_separator
from ..core.shared_functions import test_kwarg

def write_network_textfile(the_network, **kwargs):
    """ Write a simple tab-delimited textfile that can serve as
    table to import the network to cytoscape

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
     include_orphans: [False (default), True]
      whether to check for orphan nodes in the network
      include in the output.
     

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
                
    include_orphans = test_kwarg('include_orphans', kwargs, [False, True])
    the_output_dict = {}

    if continue_flag:
        the_edge_ids = []
        the_node_ids = set([])
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
                the_node_ids.update([the_node_pair[0].id, the_node_pair[1].id])

        # Try to write notes to file as individual
        # entries rather than lists when possible.
        list_note_keys = set([])
        one_element_note_keys = set([])        
        for the_edge_id in the_edge_ids:
            the_edge = the_network.edges.get_by_id(the_edge_id)
            for the_key in the_edge.notes.keys():
                if len(the_edge.notes[the_key]) < 2:
                    if the_key not in list_note_keys:
                        one_element_note_keys.add(the_key)
                else:
                    if the_key in one_element_note_keys:
                        one_element_note_keys.remove(the_key)
                    list_note_keys.add(the_key)
        
        for the_edge_id in the_edge_ids:
            the_edge = the_network.edges.get_by_id(the_edge_id)
            for the_key in the_edge.notes.keys():
                if the_key in list_note_keys:
                    the_output_dict[the_edge.id][the_key] = the_edge.notes[the_key]
                elif the_key in one_element_note_keys:
                    if len(the_edge.notes[the_key]) == 1:
                        the_output_dict[the_edge.id][the_key] = the_edge.notes[the_key][0]
                    else:
                        # Otherwise this is empty
                        the_output_dict[the_edge.id][the_key] = ''

        if include_orphans:
            the_node_location_dict = the_network.get_node_locations()
            for the_node_id in the_node_location_dict.keys():
                if the_node_id not in the_node_ids:
                    the_node = the_network.get_by_id(the_node)
                    if the_node.get_nodetype() not in exclude_nodetypes:
                        # Guess we effectively just assign the other to blank
                        the_edge_id = the_node + edge_id_separator
                        the_output_dict[the_edge_id] = {}
                        the_output_dict[the_edge_id]['node_1_id'] = the_node_id
                        the_output_dict[the_edge_id]['node_2_id'] = ''
                        # '' seems safer than 0, 
                        # ideally a NULL value would be available.
                        the_output_dict[the_edge_id]['weight'] = ''

        for the_property in properties_dict.keys():
            for the_id in properties_dict[the_property].keys():
                if the_id in the_output_dict.keys():
                    if the_id in the_edge_ids:
                        the_output_dict[the_id][the_property] = properties_dict[the_property][the_id]



        write_dict_to_textfile(the_network.id + '_network_table.txt', the_output_dict, top_key = 'model_edge_id', subfields_on_top = True)
            
        
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
     replace_lists: [True (default), False]
      if True, node notes that cannot be converted to
      1 element strings will be converted to strings
      with list items separated by '|'
     

    """
    from .networkio import write_dict_to_textfile
    from ..core import NodeType
    continue_flag = True

    replace_lists = test_kwarg('replace_lists', kwargs, [True, False])
    
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
                        if len(the_node.notes[the_key]) < 2:
                            # empty or len 1 is OK
                            if the_key not in list_note_keys:
                                one_element_note_keys.add(the_key)
                        else:
                            if the_key in one_element_note_keys:
                                one_element_note_keys.remove(the_key)
                            list_note_keys.add(the_key)
                            
            for the_node in the_nodetype.nodes:
                if the_node.id in the_node_ids:
                    for the_key in the_node.notes.keys():
                        if the_key in list_note_keys:
                            if not replace_lists:
                                the_output_dict[the_node.id][the_key] = the_node.notes[the_key]
                            else:
                                the_string_to_add = ''
                                for the_index, the_item in enumerate(the_node.notes[the_key]):
                                    if the_index > 0:
                                        the_string_to_add += '|'
                                    the_string_to_add += the_item
                                the_output_dict[the_node.id][the_key] = the_string_to_add
                                
                        elif the_key in one_element_note_keys:
                            if len(the_node.notes[the_key]) == 1:
                                the_output_dict[the_node.id][the_key] = the_node.notes[the_key][0]
                            else:
                                # Otherwise this is empty
                                the_output_dict[the_node.id][the_key] = ''

        for the_property in properties_dict.keys():
            for the_id in properties_dict[the_property].keys():
                if the_id in the_output_dict.keys():
                    the_output_dict[the_id][the_property] = properties_dict[the_property][the_id]

        write_dict_to_textfile(the_network.id + '_node_attribute_table.txt', the_output_dict, top_key = 'node_id', subfields_on_top = True)


