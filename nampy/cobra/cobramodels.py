# Some parts of this module require a functional COBRApy

    
def translate_cobra_model(cobra_model):
    """ Script to build a NAMpy model
    from a COBRApy model. The COBRApy model
    should be in memory - e.g. you need a
    functional cobrapy for this.

    Arguments:
     cobra_model

    Returns:
     network_model

                  
    """
    from nampy import Network, Node, Edge
    from copy import deepcopy

    the_network = Network(cobra_model.id)

    the_nodes_to_add = []
    for the_reaction in cobra_model.reactions:
        the_nampy_reaction = Node(the_reaction.id)
        the_nampy_reaction.notes = deepcopy(the_reaction.notes)
        the_nampy_reaction.set_nodetype('reaction')
        the_nampy_reaction.notes['reaction'] = [deepcopy(the_reaction.reaction)]
        the_nampy_reaction.notes['name'] = [deepcopy(the_reaction.name)]
        the_nampy_reaction.notes['gene_reaction_rule'] = deepcopy(the_reaction.gene_reaction_rule)
        the_nampy_reaction.notes['bounds'] = [the_reaction.lower_bound, the_reaction.upper_bound]
        the_nodes_to_add.append(the_nampy_reaction)
    
    for the_metabolite in cobra_model.metabolites:
        the_nampy_metabolite = Node(the_metabolite.id)
        the_nampy_metabolite.notes = deepcopy(the_metabolite.notes)
        the_nampy_metabolite.notes['formula'] = [deepcopy(the_metabolite.formula)]
        the_nampy_metabolite.notes['name'] = [deepcopy(the_metabolite.name)]
        the_nampy_metabolite.notes['charge'] = [deepcopy(the_metabolite.charge)]
        the_nampy_metabolite.notes['compartment'] = [deepcopy(the_metabolite.compartment)]
        the_nampy_metabolite.set_nodetype('metabolite')
        the_nodes_to_add.append(the_nampy_metabolite)

    for the_gene in cobra_model.genes:
        the_nampy_gene = Node(the_gene.id)
        the_nampy_gene.notes = deepcopy(the_gene.notes)
        the_nampy_gene.notes['name'] = deepcopy(the_gene.name)
        the_nampy_gene.set_nodetype('gene')
        the_nodes_to_add.append(the_nampy_gene)
        
    the_network.nodetypes[0].add_nodes(the_nodes_to_add)

    the_edges_to_add = []
    the_nodetype = the_network.nodetypes[0]
    # Now we need to copy the edges
    for the_metabolite in cobra_model.metabolites:
        the_nampy_metabolite = the_nodetype.nodes.get_by_id(the_metabolite.id)
        the_reactions = the_metabolite.get_reaction()
        for the_reaction in the_reactions:
            the_nampy_reaction = the_nodetype.nodes.get_by_id(the_reaction.id)
            the_edges_to_add.append([the_nampy_reaction, the_nampy_metabolite])

    for the_reaction in cobra_model.reactions:
        the_nampy_reaction = the_nodetype.nodes.get_by_id(the_reaction.id)
        the_genes = the_reaction.get_gene()
        for the_gene in the_genes:
            the_nampy_gene = the_nodetype.nodes.get_by_id(the_gene.id)
            the_edges_to_add.append([the_nampy_reaction, the_nampy_gene])

    the_network.connect_node_pair_set(the_edges_to_add)
    # We also need gene reaction rules.
    return the_network
    
            
        
    
def convert_to_cobra_model(the_network):
    """ Take a generic NAMpy model and convert to
    a COBRA model.  You need a functional
    COBRApy for this.

    Arguments:
     the_network

    kwargs:
     flux_bounds


    """
    continue_flag = True
    try:
        from cobra import Model, Reaction, Metabolite
    except:
        print 'This function requires a functional COBRApy, exiting ...'

    if continue_flag:
        __default_objective_coefficient = 0
        
        if 'flux_bounds' in kwargs:
            flux_bounds = kwargs['flux_bounds'] 
        else:
            flux_bounds = len(the_nodes)

        metabolite_dict = {}
        for the_node in the_network.Nodes: 
            the_metabolite = Metabolite(the_node.id)
            metabolite_dict.update({the_node.id: the_metabolite})

        cobra_reaction_list = []
        for the_edge in the_network.edges:
            the_reaction = Reaction(the_edge.id)
            cobra_reaction_list.append(the_reaction)
            the_reaction.upper_bound = flux_bounds
            the_reaction.lower_bound = -1 * flux_bounds
            cobra_metabolites = {}
            the_metabolite_id_1 = the_edge._nodes[0].id
            the_metabolite_id_2 = the_edge._nodes[1].id
            cobra_metabolites[metabolite_dict[the_metabolite_id_1]] = 1.
            cobra_metabolites[metabolite_dict[the_metabolite_id_2]] = -1.
            reaction.add_metabolites(cobra_metabolites)
            reaction.objective_coefficient = __default_objective_coefficient

        cobra_model = Model(model_id)
        cobra_model.add_reactions(cobra_reaction_list)

        return cobra_model
    else:
        return None


def add_flux_source(cobra_model, source_dict, **kwargs):
    """ Script to add sources to a cobra_model for
    network analysis from a dictionary

    Arguments:
     cobra_model: a cobra_model object, modified in place
     source_dict: dict of nodes to add sources/sinks to
      must include 'value' entry

    kwargs:
     match_key_type: the key in the source_dict to use to make
      matches to the model.  If 'default', ids are used.
     bound_type: 'symmetric', 'source' or 'sink'.

    Returns a tuple:
     cobra_model
     unmatched_ids_dict
                  
    """
    continue_flag = True
    try:
        from cobra import Reaction
    except:
        print 'This function requires a COBRApy, exiting...'
        continue_flag = False
    from copy import deepcopy
    __default_objective_coefficient = 0

    source_dict = deepcopy(source_dict)


    if 'match_key_type' in kwargs:
        match_key_type = kwargs['match_key_type'] 
    else:
        match_key_type = 'default'

    if 'bound_type' in kwargs:
        bound_type = kwargs['bound_type']
        if bound_type not in ['source', 'sink']:
            print('Invalid bound type, must be "source" or "sink".')
            continue_flag = False            
    else:
        bound_type = 'source'

    valid_match_key_types = ['default']
    for the_metabolite in cobra_model.metabolites:
        valid_match_key_types += the_metabolite.notes.keys()
    valid_match_key_types = list(set(valid_match_key_types))

    if match_key_type not in valid_match_key_types:
        print('Warning, %s not a valid key to match to for model nodes.  You must pick one from: %s' % (match_key_type, str(valid_match_key_types)))
        continue_flag = False

    # Note positive numbers indicate a skink
    # for exchange reactions, this is the
    # standard usually employed in COBRA
    for source_dict_id in source_dict.keys():
        the_value = source_dict[source_dict_id]['value']
        if bound_type == 'source':
            source_dict[source_dict_id]['upper_bound'] = 0.
            source_dict[source_dict_id]['lower_bound'] = -1. * the_value
        elif bound_type == 'sink':
            source_dict[source_dict_id]['upper_bound'] = 1. * the_value
            source_dict[source_dict_id]['lower_bound'] = 0.              
        
    if continue_flag:
        dict_for_pairing = {}
        for the_metabolite in cobra_model.metabolites:
            if match_key_type == 'default':
                pairing_key_list = [the_metabolite.id]
            else:
                pairing_key_list = the_metabolite.notes[match_key_type]
            for pairing_key in pairing_key_list:
                if pairing_key not in dict_for_pairing.keys():
                    dict_for_pairing[pairing_key] = {}
                    dict_for_pairing[pairing_key]['cobra_model_ids'] = []
                    dict_for_pairing[pairing_key]['source_dict_ids'] = []
                dict_for_pairing[pairing_key]['cobra_model_ids'].append(the_metabolite.id)
        for the_source_key in source_dict.keys():
            if match_key_type == 'default':
                pairing_key_list = [the_source_key]
            else:
                pairing_key_list = source_dict[the_source_key][match_key_type]
            for pairing_key in pairing_key_list:
                if pairing_key not in dict_for_pairing.keys():
                    dict_for_pairing[pairing_key] = {}
                    dict_for_pairing[pairing_key]['cobra_model_ids'] = []
                    dict_for_pairing[pairing_key]['source_dict_ids'] = []
                dict_for_pairing[pairing_key]['source_dict_ids'].append(the_source_key)

        unpaired_source_dict_keys = []
        pairing_key_model_one_to_source_dict_one = []
        unmatched_ids_dict = {}
        unmatched_ids_dict['source_dict_one_to_cobra_model_many'] = {}
        unmatched_ids_dict['cobra_model_one_to_source_dict_many'] = {}
        unmatched_ids_dict['source_dict_many_to_cobra_model_many'] = {}
        unmatched_ids_dict['source_dict_none'] = {}
        unmatched_ids_dict['cobra_model_none'] = {}
        for the_pairing_key in dict_for_pairing.keys():
            the_mapping_type = ''
            cobra_model_ids = dict_for_pairing[the_pairing_key]['cobra_model_ids']
            source_dict_ids = dict_for_pairing[the_pairing_key]['source_dict_ids']
            if ((len(cobra_model_ids) == 1) & (len(source_dict_ids) == 1)):
                pairing_key_model_one_to_source_dict_one.append(the_pairing_key)
            elif ((len(cobra_model_ids) > 1) & (len(source_dict_ids) == 1)):
                the_mapping_type = 'source_dict_one_to_cobra_model_many'
                unpaired_source_dict_keys += source_dict_ids
            elif ((len(cobra_model_ids) == 1) & (len(source_dict_ids) > 1)):
                the_mapping_type = 'cobra_model_one_to_source_dict_many'
                unpaired_source_dict_keys += source_dict_ids
            elif ((len(cobra_model_ids) > 1) & (len(source_dict_ids) > 1)):
                the_mapping_type = 'cobra_model_many_to_source_dict_many'
                unpaired_source_dict_keys += source_dict_ids
            elif (len(cobra_model_ids) == 0):
                the_mapping_type = 'cobra_model_none'
                unpaired_source_dict_keys += source_dict_ids
            elif (len(source_dict_ids) == 0):
                the_mapping_type = 'source_dict_none'

            if len(the_mapping_type) >0:
                unmatched_ids_dict[the_mapping_type][the_pairing_key] = {}
                unmatched_ids_dict[the_mapping_type][the_pairing_key]['cobra_model_ids'] = cobra_model_ids
                unmatched_ids_dict[the_mapping_type][the_pairing_key]['source_dict_ids'] = source_dict_ids

        unpaired_source_dict_keys = list(set(unpaired_source_dict_keys))

        cobra_reaction_list = []
        for the_pairing_key in pairing_key_model_one_to_source_dict_one:
            cobra_model_id = dict_for_pairing[the_pairing_key]['cobra_model_ids'][0]
            source_dict_id = dict_for_pairing[the_pairing_key]['source_dict_ids'][0]
            the_reaction_id = 'EX_' + cobra_model_id
            reaction = Reaction(the_reaction_id)
            cobra_reaction_list.append(reaction)
            if 'upper_bound' in source_dict[source_dict_id].keys():
                reaction.upper_bound = source_dict[source_dict_id]['upper_bound']
            else:
                reaction.upper_bound = __default_upper_bound
            if 'lower_bound' in source_dict[source_dict_id].keys():
                reaction.lower_bound = source_dict[source_dict_id]['lower_bound']
            else:
                reaction.lower_bound = __default_lower_bound
            if 'objective_coefficient' in source_dict[source_dict_id].keys():
                reaction.objective_coefficient = source_dict[source_dict_id]['objective_coefficient']
            else:
                reaction.objective_coefficient = __default_objective_coefficient
                
            cobra_metabolite = cobra_model.metabolites.get_by_id(cobra_model_id)
            reaction.add_metabolites({cobra_metabolite: -1})
            
        print("Completed adding sources with %i pairings of %i in the source_dict and %i in the model." % (len(source_dict.keys()) - len(unpaired_source_dict_keys), len(source_dict.keys()), len(cobra_model.metabolites)))
            
        cobra_model.add_reactions(cobra_reaction_list)
        return cobra_model, unmatched_ids_dict
    else:
        return None, None

