# Some parts of this module require a functional COBRApy

    
def translate_cobra_model(cobra_model):
    """ Script to build a NAMpy model
    from a COBRApy model. The COBRApy model
    should be in memory - e.g. you need a
    functional COBRApy for this.

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
    a COBRA model.  The model is assumed to be monopartite.
    You need a functional COBRApy for this.

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
        for the_node in the_network.nodetypes[0].nodes: 
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


