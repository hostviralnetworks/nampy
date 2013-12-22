
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


