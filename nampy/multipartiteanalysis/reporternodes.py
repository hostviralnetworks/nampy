

def evaluate_reaction_pvalues(the_network, expression_dict, **kwargs):
    """ pvalues for gene changes are combined
    to get reaction level changes

    this function inspired by panalty evaluation in gim3e:
    Schmidt, B. J., Ebrahim, A., Metz, T. O., Adkins, J. N., Palsson, B. O., & Hyduke, D. R. (2013). 
    GIM3E: condition-specific models of cellular metabolism developed from metabolomics 
    and expression data. Bioinformatics, 29(22), 2900–2908.

    Arguments:
     the_network: a nampy network from a COBRApy model, 
                  converted to multipartite
     expression_dict: a dict of pvalues

    Returns:
     penalites: a dict with reaction id's as keys
                and two-tailed test for change
                p-values values

    
    """
    import re
    from copy import deepcopy
    # # First we test for gene P/A calls.
    # # Declare P/A based on cutoff
    # gene_pa_dict = {}
    # for cur_gene in expression_dict.keys():
    #     if cur_gene in new_cobra_model.genes:
    #         if expression_dict[cur_gene] > threshold:
    #             gene_pa_dict[cur_gene] = 1
    #         else:
    #             gene_pa_dict[cur_gene] = 0

    allowed_settings = ['no_change', 'ignore_reaction']
    if 'missing_genes' in kwargs:
        if kwargs['missing_genes'] in allowed_settings:
            missing_genes = kwargs['missing_genes']
            # ignore_reaction not yet implemented
            missing_genes = allowed_settings[0]
        else:
            missing_genes = allowed_settings[0]
    else:
        missing_genes = allowed_settings[0]

    # settings for reactions with missing gene_reaction_rules?
    # seems only fair that we should ignore these reactions.

    reaction_pval_dict = {}
    # We also keep a list to track if a reaction is inactivated
    # due to GE data
    reaction_inactivated = []
    for test_reaction in the_network.nodetypes.get_by_id('reaction').nodes:
        # Ignore reactions without known gene reaction relations.
        # len(test_reaction.notes['gene_reaction_rule']) == 0:
        # penalties[test_reaction.id] = 0
        if (len(test_reaction.notes['gene_reaction_rule']) > 0) & (test_reaction.notes['gene_reaction_rule'] != 's0001'):
            # Otherwise, we will have to evaluate the GRR
            the_gene_reaction_rule = deepcopy(test_reaction.notes['gene_reaction_rule'])
            reaction_gene_dict = {}
            current_gene_list = []
            for the_edge in test_reaction.get_edges():
                the_node_pair = the_edge.get_node_pair()
                the_index = the_node_pair.index(test_reaction)
                if the_index == 1:
                    if the_node_pair[0].get_nodetype() == 'gene':
                        current_gene_list.append(the_node_pair[0])
                else:
                    if the_node_pair[1].get_nodetype() == 'gene':
                        current_gene_list.append(the_node_pair[1])                    
            
            for the_gene in current_gene_list:
                if the_gene.id in expression_dict:
                    reaction_gene_dict[the_gene.id] = expression_dict[the_gene.id]
                else:
                    if missing_genes == 'no_change':
                        reaction_gene_dict[the_gene.id] = 1.
                    else:
                        # not yet implemented
                        reaction_gene_dict[the_gene.id] = 1.
                        
            # Now evaluate the reaction
            for the_gene in current_gene_list:
                the_gene_re = re.compile('(^|(?<=( |\()))%s(?=( |\)|$))'%re.escape(the_gene.id))
                the_gene_reaction_rule = the_gene_re.sub(str(reaction_gene_dict[the_gene.id]), the_gene_reaction_rule)
            reaction_pval_dict[test_reaction.id] = evaluate_expression_string(the_gene_reaction_rule)
    return reaction_pval_dict

     
import re
number_finder = re.compile("[\d]+\.?[\d]*")

class penalty_number:
    """
    this class inspired by panalty evaluation in gim3e:
    Schmidt, B. J., Ebrahim, A., Metz, T. O., Adkins, J. N., Palsson, B. O., & Hyduke, D. R. (2013). 
    GIM3E: condition-specific models of cellular metabolism developed from metabolomics 
    and expression data. Bioinformatics, 29(22), 2900–2908.
    """
    def __init__(self, value):
        self.str_value = str(value)
        self.value = float(value)
    
    def __add__(self, other):
        # addition is like OR
        return penalty_number(min(self.value, other.value))
    
    def __mul__(self, other):
        # multiplication is like AND
        return penalty_number(max(self.value, other.value))


def evaluate_expression_string(penalty_string, **kwargs):
    """ penalty string will have:
        * 'or' statements which need to be converted to min
        * 'and' statements which need to be converted to max
    >>> evaluate_penalty("(1 and 2 and 3)")
    max(1, 2, 3)

    this function inspired by panalty evaluation in gim3e:
    Schmidt, B. J., Ebrahim, A., Metz, T. O., Adkins, J. N., Palsson, B. O., & Hyduke, D. R. (2013). 
    GIM3E: condition-specific models of cellular metabolism developed from metabolomics 
    and expression data. Bioinformatics, 29(22), 2900–2908.

    """
    # if there are no ands or ors, we are done
    penalty_string = penalty_string.lower()  # don't want to deal with cases
    
    if "and" not in penalty_string and "or" not in penalty_string:
        return eval(penalty_string)
    # we will replace AND and OR with addition and multiplication
    # equivalent to max/min
    penalty_string = penalty_string.replace("or", "*").replace("and", "+")
    # replace the numbers with the custom class which have overloaded AND/OR
    values = [penalty_number(i) for i in number_finder.findall(penalty_string)]
    values_strings = tuple("values[%i]" % i for i in range(len(values)))
    penalty_string = number_finder.sub("%s", penalty_string)
    penalty_string = penalty_string % values_strings
    return eval(penalty_string).value


def calculate_reporter_scores(the_network, hyperedge_score_dict, score_nodetype, hyperedge_nodetype, **kwargs):
    """ This algorithm is based on the reporter methods
    published by Jens Nielsen's group:

    Patil, K. R., & Nielsen, J. (2005). Uncovering transcriptional 
    regulation of metabolism by using metabolic network topology. PNAS, 102(8), 2685–9.

    Oliveira, A. P., Patil, K. R., & Nielsen, J. (2008). Architecture of transcriptional regulatory 
    circuits is knitted over the topology of bio-molecular interaction networks. BMC systems biology, 
    2, 17.

    And also inspired by functions in COBRApy
    Ebrahim, A., Lerman, J. A., Palsson, B. O., & Hyduke, D. R. (2013). COBRApy: 
    COnstraints-Based Reconstruction and Analysis for Python. BMC systems biology, 7(1), 74.

    Arguments:
     the_network: a NAMpy network with appropriately defined nodetypes
     hyperedge_score_dict: scores for the hyperedges.  here, the nodes
                           scored in the source hyperedge_score_dict are called
                           hyperedges
     score_nodetype: the type of node to score, e.g. "reaction," "metabolite," "gene"
     hyperedge_nodetype: the 'nodetype' the source scores are for, e.g. "reaction," "metabolite," "gene"

    kwargs:
     score_type: 'p' or 'z'
     number_of_randomizations
     background_correction
     degree: aggregate all n-degree neighbors into the reporter metabolites.  
             The default, 1, captures nearest-neighbors
     verbose: If true, output information on the background correction calculation,
              which can take some time when the number_of_randomizations is large or
              the degree is high.
    aggregation_type: 'size-independent' or 'mean'
                      'size-independent' is similar to Stouffer's method, z-scores divided by 1/(k**.5)
                      and was employed in Patil 2005.  This method is not supported for > degree 1
                      'mean' aggregates z-scores by the average and was employed in Oliveira 2008
                      note that if background_correct is on, these methods will yield identical,
                      ~size-independent results.

    
    """
    from scipy.stats import norm
    from numpy.random import randint
    from numpy import mean, std, isnan, array, where, zeros, inf
    from time import time
    from math import floor, ceil
    from ..core.Node import Node
    continue_flag = True

    allowed_values = ['p', 'z']
    if 'score_type' in kwargs:
        if kwargs['score_type'] in allowed_values:
            score_type = kwargs['score_type']
        else:
            score_type = allowed_values[0]
    else:
        score_type = allowed_values[0]

    if 'number_of_randomizations' in kwargs:
        number_of_randomizations = kwargs['number_of_randomizations']
    else:
        number_of_randomizations = 1000

    allowed_values = [True, False]
    if 'background_correction' in kwargs:
        if kwargs['background_correction'] in allowed_values:
            background_correction = kwargs['background_correction']
        else:
            background_correction = allowed_values[0]
    else:
        background_correction = allowed_values[0]

    if 'degree' in kwargs:
        degree = kwargs['degree']
    else:
        degree = 1

    allowed_values = [True, False]
    if 'verbose' in kwargs:
        if kwargs['verbose'] in allowed_values:
            verbose = kwargs['verbose']
        else:
            verbose = allowed_values[0]
    else:
        verbose = allowed_values[0]        

    allowed_values = ['size-independent', 'mean']
    if 'aggregation_type' in kwargs:
        if kwargs['aggregation_type'] in allowed_values:
            aggregation_type = kwargs['aggregation_type']
        else:
            aggregation_type = allowed_values[0]
    else:
        aggregation_type = allowed_values[0] 

    allowed_values = [x.id for x in the_network.nodetypes]
    if score_nodetype not in allowed_values:
        print 'A valid nodetype is not selected for score_nodetype, exiting...'
        confinue_flag = False

    allowed_values = [x.id for x in the_network.nodetypes]
    if score_nodetype not in allowed_values:
        print 'A valid nodetype is not selected for score_nodetype, exiting...'
        confinue_flag = False
    if hyperedge_nodetype not in allowed_values:
        print 'A valid nodetype is not selected for hyperedge_nodetype, exiting...'
        confinue_flag = False
    
    if ((aggregation_type == 'size-independent') & (degree > 1)):
        print("Warning, 'size-independent' aggregation option has not been investigated when degree > 1, the background correction results may be suboptimal.")

    if continue_flag:

        # Make sure these are id's, will use to grab objects from the current model
        from copy import deepcopy
        hyperedge_score_dict = deepcopy(hyperedge_score_dict)
        for the_key in hyperedge_score_dict.keys():
            if type(the_key) == Node:
                hyperedge_score_dict[the_key.id] = hyperedge_score_dict[the_key]
                hyperedge_score_dict.pop(the_key)
            
        model_node_ids = [x.id for x in the_network.nodetypes.get_by_id(score_nodetype).nodes]
        model_hyperedge_ids =  [x.id for x in the_network.nodetypes.get_by_id(hyperedge_nodetype).nodes]

        hyperedge_list = [the_network.nodetypes.get_by_id(hyperedge_nodetype).nodes.get_by_id(hyperedge_id) for hyperedge_id in hyperedge_score_dict.keys() if hyperedge_id in model_hyperedge_ids]
        the_scores = [hyperedge_score_dict[the_hyperedge.id] for the_hyperedge in hyperedge_list]

        # minimum and maximum p-values are used to prevent numerical problems.
        the_scores = array(the_scores)
        if score_type == 'p':
            minimum_p = min(the_scores[the_scores.nonzero()[0]])
            maximum_p = max(the_scores[where(the_scores < 1)[0]])
            the_scores[where(the_scores < minimum_p)] = minimum_p
            the_scores[where(the_scores > maximum_p)] = maximum_p
            the_scores = norm.ppf(the_scores)
        elif score_type == 'z':
            minimum_z = min(the_scores[where(the_scores>-inf)])
            maximum_z = max(the_scores[where(the_scores<inf)])
            the_scores[where(the_scores < minimum_z)] = minimum_z
            the_scores[where(the_scores > maximum_z)] = maximum_z
        # update the dictionary with the new scores
        hyperedge_score_dict = dict(zip(hyperedge_list, the_scores))

        # Get the connectivity for each node of interest,
        # where we have data for at least 1 connected hyperedge
        the_nodes = set([])
        for the_hyperedge in hyperedge_score_dict.keys():
            the_connected_nodes = the_hyperedge.get_connected_nodes()
            for the_node in the_connected_nodes:
                if the_node.get_nodetype() == score_nodetype:
                    the_nodes.add(the_node)
        # Observe that there may be some asymmetry in the node relations.
        # E.g. a metabolite can still turnover/change if at least
        # one reaction is present, but a reaction cannot proceed if
        # a single metabolite is missing/ cannot sustain turnover.
        # However, we simply score based on available connected p-values here

        # Calculate the score for each metabolite
        start_time = time()
        # To speed up calculations, get a list of all
        # neighbors for each node
        node_neighbor_dict = {}
        # We also need to calculate the uncorrected scores 
        # for the first degree
        # neighboring nodes in any case
        first_degree_node_scores = {}
        neighbor_hyperedge_dict = {}
        for the_node in the_nodes:
            # Here, we count a hyperedge as a neighbor only if
            # the connection linking it to the present node is scored
            neighbor_hyperedges = set([x for x in the_node.get_connected_nodes() if (x in hyperedge_list)])
            neighbor_hyperedge_dict[the_node.id] = [x.id for x in neighbor_hyperedges]
            neighbor_nodes = set([])
            for the_hyperedge in neighbor_hyperedges:
                # Only nodes connected to some associated data, already in 
                # the_nodes, will be useful for calculating z-scores
                neighbor_nodes = neighbor_nodes.union(set([x.id for x in the_hyperedge.get_connected_nodes() if x in the_nodes]))
            neighbor_nodes.remove(the_node.id)
            node_neighbor_dict[the_node.id] = neighbor_nodes
            number_of_connections = len(neighbor_nodes)
            tmp_score = 0.
            for the_hyperedge in neighbor_hyperedges:
                tmp_score += hyperedge_score_dict[the_hyperedge]
            if aggregation_type == "size-independent":
                first_degree_node_scores[the_node.id] = tmp_score / (number_of_connections ** 0.5)
            else:
                # only other option at this time is "mean"
                first_degree_node_scores[the_node.id] = tmp_score / number_of_connections            

        node_scores, node_neighbors = calculate_node_scores(node_neighbor_dict, first_degree_node_scores, degree, aggregation_type = aggregation_type)
    
        if verbose:
            print("Time to identify connections and calculate scores: " + str(time() - start_time))

        if background_correction:
            start_time = time()
            # safe to eliminate 0, these metabolite scores are 0.
            #test_size_list = list(set(metabolite_connections_n.values()).difference(set([0])))
            #test_size_list.sort()
            the_node_ids = node_scores.keys()
            first_degree_node_random_scores = zeros((len(the_node_ids), number_of_randomizations))
            for node_index, the_node_id in enumerate(the_node_ids):
                # The randomization for degree = 1 is pretty straightforward,
                # we just need a background based on the number of connections.
                # However, for higher degrees, the configuration
                # may also play a role in the score.  In this case we need to generate a
                # score distribution for the connections of each first degree node,
                # then we can deduce a score distribution for the
                # neighbor set of each metabolite.
                # Basically, what we're doing here to start is that for each i we select i
                # scores number_of_randomizations times
                i = len(node_neighbor_dict[the_node_id])
                the_random_indices = randint(0, len(the_scores), size=(number_of_randomizations, i))
                
                if aggregation_type == "size-independent":
                    random_score_distribution = (the_scores[the_random_indices]).sum(1) / (i ** 0.5)
                else:
                    # only other option at this time is "mean"
                    random_score_distribution = (the_scores[the_random_indices]).sum(1) / i
                # keep this as a matrix for now
                first_degree_node_random_scores[node_index,:] = random_score_distribution

            node_random_distribution_array = zeros((len(the_node_ids), number_of_randomizations))
            for randomization_index in range(0, number_of_randomizations):
                current_node_random_distribution_dict = {x: y for x, y in zip(the_node_ids, first_degree_node_random_scores[:, randomization_index])}
                # To save memory we delete this column after calculating the scores.
                random_node_scores, temp = calculate_node_scores(node_neighbor_dict, current_node_random_distribution_dict, degree, aggregation_type = aggregation_type)
                for index, the_node_id in enumerate(the_node_ids):
                    node_random_distribution_array[index, randomization_index] = random_node_scores[the_node_id]
                
                if verbose:
                    if randomization_index < (number_of_randomizations - 1):
                        if (randomization_index + 1) % 100 == 0:
                            scaling_factor = ((number_of_randomizations - randomization_index - 1) / (randomization_index + 1))
                            pass_s = time() - start_time
                            remaining_s = pass_s * scaling_factor
                            remaining_h = floor((remaining_s)/3600)
                            remaining_m = floor(((remaining_s)/3600 - remaining_h) * 60)
                            print("Calculated randomization "+ str(randomization_index + 1) +" of " + str(number_of_randomizations) + ".  Remaining time estimate: %0.0f hr %0.0f min." % (remaining_h, remaining_m))            
            for index, the_node_id in enumerate(the_node_ids):
                random_score_distribution = node_random_distribution_array[index, :]
                node_scores[the_node_id] = (node_scores[the_node_id] - mean(random_score_distribution)) / std(random_score_distribution)

        # For summary purposes also return information on metabolite and reaction neighbors
        node_neighbor_hyperedge_n = {}
        node_neighbor_n = {}
        node_pvalue_dict = {}
        for the_node_id in node_neighbors.keys():
            node_neighbor_n[the_node_id] = len(node_neighbors[the_node_id])
            # Note that large positive z-scores correspond to reporters that change
            # more than the mean.  We effectively want a one-tail test
            # here but we assign smaller p-values to larger z-scores.
            # Note the reporter z-scores are in fact normally discributed due to
            # the background correction.
            if background_correction:
                node_pvalue_dict[the_node_id] = norm.cdf(-1 * node_scores[the_node_id])
            else:
                node_pvalue_dict[the_node_id] = None
            
        return_dictionary = {'z_scores': node_scores,
                         'n_neighbors':  node_neighbor_n,
                         'directly_connected_hyperedges': neighbor_hyperedge_dict,
                         'p_values': node_pvalue_dict}

        return return_dictionary
    else:
        return None


def calculate_node_scores(node_neighbor_dict, first_degree_node_scores, degree, **kwargs):
    """ Score the nodes.
    
    Arguments:
     node_neighbor_dict: a dictionary of node_id : set(id_of_neighbor)
      where id_of_neighbor are nearest node neighbors (not hyperedges)
     first_degree_metabolite_scores: scores for nodes
      (no neighbors included) of the form node_id: score
     degree: degree of neighboring nodes to consider

    kwargs:
     aggregation_type

    Returns:
     node_scores: scores aggregated up to a specified degree neighbors
     node_neighbors: neighboring emtabolites of up to the specified degree
     

    """
    allowed_values = ['size-independent', 'mean']
    if 'aggregation_type' in kwargs:
        if kwargs['aggregation_type'] in allowed_values:
            aggregation_type = kwargs['aggregation_type']
        else:
            aggregation_type = allowed_values[0]
    else:
        aggregation_type = allowed_values[0] 

    
    the_nodes = node_neighbor_dict.keys()
    node_scores = {}
    node_neighbors = {}
    for key_node_id in the_nodes:
        for current_degree in range(1, degree + 1):
            if current_degree == 1:
                all_neighbor_nodes = set([key_node_id])
                neighbor_nodes_tested = set([])
            else:
                last_version = list(all_neighbor_nodes)
                [all_neighbor_nodes.update(node_neighbor_dict[x]) for x in last_version]
        all_neighbor_nodes = list(all_neighbor_nodes)
        neighbor_scores = []
        neighbor_scores = [first_degree_node_scores[x] for x in all_neighbor_nodes]
        if aggregation_type == "size-independent":
            node_scores[key_node_id] = sum(neighbor_scores) / (len(neighbor_scores) ** 0.5)
        else:
            # only other option at this time is "mean"
            node_scores[key_node_id] = sum(neighbor_scores) / len(neighbor_scores)
        node_neighbors[key_node_id] = all_neighbor_nodes
    return node_scores, node_neighbors

