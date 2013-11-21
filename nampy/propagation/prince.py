def prince(the_network, **kwargs):
    """
    Based on Vanunu, O., Magger, O., Ruppin, E., Shlomi, T., & Sharan, R. (2010). 
    Associating genes and protein complexes with disease via network propagation. 
    PLoS computational biology, 6(1), e1000641. doi:10.1371/journal.pcbi.1000641
    Credits for the nampy translation: Dorothea Emig-Agius, Greg Hannum, Brian Schmidt, 2013

	flow, pump information through network from source nodes  
	 
	formula:
	F(u)^t = alpha * W' * F(u)^t-1 + (1-alpha) * Y
	
	W': adjmatrix / weightmatrix with entries normalized by row sums (degree-normalized or weight-normalized)
	W'(i,j) = W(i,j) / sqrt(D(i,i), D(j,j)) 
	D(i,i) = row sum_j: wij
	
	Y = vector with prior knowledge, all source nodes set to 1, rest to 0
	
	alpha > 0.5, apparently does not make any difference, thus set it to 0.8
	iterate and when everything is in steady-state add a round with alpha = 1


    Arguments: 
     the_network

    kwargs:
     alpha
     verbose

    Returns:
     propagation_score: a dict of {node_id: scores}


    """
    from scipy import array
    from numpy import sqrt, zeros, dot

    the_dim = len(the_network.nodes)

    # Make sure the matrix is updated
    the_network.update()
    the_matrix = the_network.matrix

    if 'alpha' in kwargs: 
        alpha = kwargs['alpha']
    else:
        alpha = 0.8

    if 'verbose' in kwargs: 
        verbose = kwargs['verbose']
    else:
        verbose = False

    if verbose:
        print "Running PRINCE."

    diagonal = the_matrix.sum(1)

    the_norm = sqrt(diagonal * diagonal.transpose())

    # Better to perform computations on dense
    w_prime = the_matrix.todense() / the_norm
    if verbose:
        print "W' done."

    # now set Y
    # set all source nodes to 1, rest to 0
    ft1 = zeros((the_dim, 1))
    y = zeros((the_dim, 1))
    for i, the_node in enumerate(the_network.nodes):
        ft1[i] = the_node.source
        y[i] = ft1[i] * (1.0 - alpha)
    if verbose:
        print "Y done."

    # now compute the propagation...
    ft = zeros((the_dim, 1))
    continue_propagation = True 
    l1norm_cutoff = 1E-6
    
    while continue_propagation:
        ft = alpha * dot(w_prime,ft1) + y
        if (abs(ft - ft1)).sum() < l1norm_cutoff:
            continue_propagation = False
        else:
            ft1 = ft
    print "Transition done"

	# and now the final smoothing step with alpha = 1... (which means, no y added in the last step)
    ft = 1. * dot(w_prime, ft)

    result_dict = {}
    for i, the_node_id in enumerate([x.id for x in the_network.nodes]):
        result_dict[the_node_id] = ft[i,0]

    return result_dict
