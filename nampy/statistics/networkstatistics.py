# Full functionality of this module requires a functional rpy2 and R
from ..core.shared_functions import test_kwarg


def mtcorrect(p_value_dict, method = '"BY"'):
    """ Apply MT correction.  This is a wrapper for R's p.adjust function.

    Arguments:
     p_value_dict: a dict with keys = probe names, and values of p-values
     method: MT correction method, from R.
      "holm"
      "hochberg",
      "hommel"
      "bonferroni"
      "BH"
      "BY"
      "fdr"
      "none"
    Returns:
     adjusted_p
     
 
    """
    continue_flag = True
    try:
        from rpy2.robjects import r
        from rpy2 import robjects
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()
    except:
        print "This module requires a functional rpy2 and R, exiting..."
        continue_flag = False

    if continue_flag:
        row_names = [x for x in p_value_dict.keys()]
        p_values_list = [p_value_dict[id] for id in row_names]
        # need to create an r object first
        p_values_list_r = robjects.FloatVector(p_values_list)
        # need to assign the r object into the r namespace
        r.assign('p_values_list_r', p_values_list_r)
        r('corrected_data = p.adjust(p_values_list_r, method = ' + str(method) + ')')
        adjusted_p = robjects.numpy2ri.ri2numpy(r('corrected_data'))
        adjusted_p.tolist
        adjusted_p = {id: adjusted_p[i] for i, id in enumerate(row_names)}
        return adjusted_p
    else:
        return {}


    
def get_pvalue_from_scores(result_dict, permutation_dict, **kwargs):
    """ Uses network scores to calculate p-values.

    Arguments:
     result_dict: a dict of results from the network propagation.  e.g. {node id : score value}
     permutation_dict: a dict of permutation results, e.g. {node id: [list of score values]}

    kwargs:
     verbose

    returns:
     ttp_dict: dictionary of p-values from a two-tailed test of significant
     side_dict: a dict to indicate if the result lies in the upper or lower tail ('+' or '-')
 
    """
    continue_flag = True

    if 'verbose' in kwargs: 
        verbose = kwargs['verbose']
    else:
        verbose = False
    
    try:
        from rpy2.robjects import r
        from rpy2 import robjects
        from rpy2.robjects import numpy2ri
        numpy2ri.activate()
    except:
        print "This module requires a functional rpy2 and R, exiting..."
        continue_flag = False

    if continue_flag:
        from time import time
        start_time = time()
    
        ttp_dict = {}
        side_dict = {}
    
        # Also require the pareto extrapolation code
        the_file = 'ParetoExtrapolation.R'
        import nampy
        import platform 
        # This may need some fix, I've only tried
        # this on OSX so far.
        if platform.system() == 'Windows':
            pareto_file_path = nampy.__path__[0] + '\\rfiles\\'
        else:
            pareto_file_path = nampy.__path__[0] + '/rfiles/'
        r('source(%s)' %('"' + pareto_file_path + the_file + '"'))
        test_ids = result_dict.keys()
        n_tests = len(test_ids)
        start_time = time()
        for i, the_id in enumerate(test_ids):
            test_statistic = result_dict[the_id]
            null_distribution = permutation_dict[the_id]
            n_permutations = len(null_distribution)
            # This won't work by definition if we have less than 20 permutations
            if n_permutations < 20:
                if verbose:
                    print 'Warning, null distribution is too small.  Run more permutations.  Exiting...'
                return (None, None)
            M1 = sum(null_distribution > test_statistic)
            M2 = sum(null_distribution < test_statistic)
            # Want to know which end the test statistics lies towards
            if M2 < M1:
                side = '-'
                M = M2
            else:
                side = '+'
                M = M1
            if (M >= 10):
                estimated_p = 2. * M / n_permutations
            else:
                # need to create an r object first
                null_distribution_r = robjects.FloatVector(null_distribution)
                # need to assign the r object into the r namespace
                r.assign('null_distribution_r', null_distribution_r)
                r.assign('test_statistic', test_statistic)
                # now run our r scripts
                r('fit = getParetoFit(null_distribution_r, side="two-sided")') 
                r('distCDF = paretoExtrapolate(test_statistic, fit)')
                estimated_p = robjects.numpy2ri.ri2numpy(r('distCDF'))[0]
            ttp_dict[the_id] = estimated_p
            side_dict[the_id] = side
            if verbose:
                if (i+1) % 1000 == 0:
                    print 'Test %i of %i, el: %f hr' %((i + 1), n_tests, (time() - start_time)/3600.)
        return ttp_dict, side_dict
    else:
        return {}, {}


def ttp_to_otp(ttp_dict, side_dict):
    """ Take p-values from a two-tailed test and
    convert to values from a one-tailed test.

    Note the convention assumed here is that
    increases correspond to p-values near zero
    in the one-tailed test, so you this is 
    technically backwards.  This is just more relevant if
    we are looking for increases in a network score.

    Arguments:
     ttp_dict: a dictionary of {id: two-tail p-value}
     side_dict: a dictionary of {id: '+' or '-'}

    Returns:
     otp_dict
     
    """
    
    otp_dict = {}
    for the_key in ttp_dict.keys():
        the_dir = side_dict[the_key]
        if the_dir == '+':
            otp_dict[the_key] = ttp_dict[the_key] / 2.
        else:
            otp_dict[the_key] = 1. - ttp_dict[the_key] / 2.
            
    return otp_dict


def t_uneqvar(list_1, list_2, **kwargs):
    """ Performs a t-test without the equal variance
    assumption of Student's t.  For example, see: 
    Ruxton, G. D. (2006). 
    The unequal variance t-test is an underused alternative to 
    Student's t-test and the Mann-Whitney U test
    Behavioral Ecology, 17(4), 688–690.

    Arguments:
     list_1, list_2: list of values from the first and final condition, respectively

    Returns: a dict containing keys:
     'p': the p-value resulting from a two-tailed test for change
          note two-tailed p-values are preferred to avoid numerical issues
          associated with highly significant p-values
     'dir': the direction
     't': the t statistic
     'df': the calculated degrees of freedom

    """
    from scipy.stats import t 
    from numpy import std, mean
    from math import sqrt

    two_tailed = test_kwarg('two_tailed', kwargs, [True, False])
    the_return_dict = {}

    var_1 = (std(list_1, ddof = 1))**2
    var_2 = (std(list_2, ddof = 1))**2
    the_u = var_2 / var_1
    n_1 = len(list_1) * 1.
    n_2 = len(list_2) * 1.
    df = (1./n_1 + the_u/n_2)**2/(1/(n_1**2*(n_1-1)) + the_u**2/(n_2**2*(n_2-1)))
    # Use 1 - 2 here before calculating p so more positive changes corresponse to smaller p
    t_val = (mean(list_1) - mean(list_2)) / sqrt((var_1 / n_1) + (var_2 / n_2))
    # One-sided p
    the_p = t.cdf(t_val, df)
    t_val = -1. * t_val
    if two_tailed:
        if t_val > 0:
             the_p = 2. * the_p
             the_dir = '+'
        else:
            # It is numerically preferable to avoid 
            # y = 2 * (1 - x) in case x is close to zero
            the_p = 2. * t.cdf(t_val, df)
            the_dir = '-'
    return_dict = {'t': t_val, 'p': the_p, 'df': df, 'dir': the_dir}
    return return_dict


def stouffer_p_agg(p_list, **kwargs):
    """Use Stouffer's method to aggregate a list of p-values.
    References :  M.C. Whitlock 2005 J Ev Biology
                  D.V. Zaykin 2011 J Evol Biol
    Arguments:
     p_list: list of pvalues, assumed from a one-sided test with 
      p-values near 0 corresponding to increases in
      expression, but if unless dir_list is given
      then it is assumed p-values are from two-sided tests

    kwargs:
     dir_list: a list of '+' / '-' for the changes
     two_tailed: True (default) / False
                True -> the p-value returned is equivalent
                        to a 2-tailed test for change.  
                False -> the p-value returned is equivalent
                        to a 1-tailed test for increase.
                        Here, the convention is adapted that
                        more positive changes correspond to
                        smaller p-values
     verbose: True (default) / False
     
    Returns:
     A dict of {'p': the_p, 'dir': the_dir}
      p: Returns the two_tailed p-value

    """
    from scipy.stats import norm
    from copy import deepcopy
    continue_flag = True

    verbose = test_kwarg('verbose', kwargs, [True, False])

    p_list = deepcopy(p_list)

    if 'dir_list' in kwargs:
        dir_list = kwargs['dir_list']
        for the_entry in dir_list:
            if continue_flag:
                if the_entry not in ['+','-']:
                    continue_flag = False
                    if verbose:
                        print "Cannot calculate p-value, need '+' or '-' for direction."
    else:
        dir_list = []

    if continue_flag:
        z_agg = stouffer_z_agg(p_list, dir_list = dir_list)
        pval = norm.cdf(-1 * z_agg)
        if pval <= 0.5:
            the_dir = '+'
        else:
            the_dir = '-'
        # We require the returned p-value
        # to be two-tailed
        if the_dir == '+':
            pval = pval * 2
        else:
            pval = 2 * norm.cdf(z_agg)
    else:
        pval = None
        the_dir = None
                
    return {'p': pval, 'dir': the_dir}


def stouffer_z_agg(p_list, **kwargs):
    """ Use Stouffer's method to aggregate a list of p-values and return the z-value.
    References :  M.C. Whitlock 2005 J Ev Biology
                  D.V. Zaykin 2012 J Evol Biol

    Arguments:
     p_list: list of pvalues to aggregate. By default assumed from a 
      one-sided test for increase: values near
      zero are assumed to corespond to increases, near 0.5 no change,
      and near 1 decreases. If the dir_list kwarg is given, then this
      assumed to be the result of a two-tail test (preferred).

    kwargs:
     dir_list: optional list of '+' / '-' corresponding to 
      to the list of p-values to determine directionality.
    Returns:
     z: the z-value, positive numbers imply increases
        and negative numbers imply decreases

    """
    from scipy.stats import norm
    from numpy import nan, inf
    from copy import deepcopy

    continue_flag = True

    verbose = test_kwarg('verbose', kwargs, [True, False])
    
    # avoid overwriting the original list in memory
    p_list = deepcopy(p_list)

    if 'dir_list' in kwargs:
        dir_list = kwargs['dir_list']
        if (len(dir_list) != 0) & (len(dir_list) != len(p_list)):
            continue_flag = False
            if verbose:
                print "Cannot calculate z-score, dir_list must be same length as p_list."
        for the_entry in dir_list:
            if the_entry not in ['+','-']:
                continue_flag = False
                if verbose:
                    print "Cannot calculate z-score, need '+' or '-' for direction."
    else:
        dir_list = []

    for the_p in p_list:
        if ((the_p < 0) | (the_p > 1)):
            continue_flag = False
            if verbose:
                print "Cannot calculate z-score, p values must be between 0 and 1."
  
    if continue_flag:
        if len(dir_list) > 0:
            # Note for the moment we neglect the sign in two-tailed tests 
            # to avoid numerical errors caused by p values just smaller than 1 
            # being rounded up.  A sign is assigned once the z-value is calculated
            p_list = [the_p * 0.5 for the_p in p_list]
        if ((min(p_list) > 0) & (max(p_list) < 1)):
            z = 0.
            for the_index in range(0, len(p_list)):
                # Take advantage of symmetry to avoid numerical issues
                pval = p_list[the_index]
                znew = -1*norm.ppf(pval)
                # Now we correct the sign if this is originating from a two-tail test
                # Force the convention that increases lead to positive z values
                if len(dir_list) > 0:
                    if dir_list[the_index] == '-':
                        znew = znew * -1.
                z += znew
            z = z / (len(p_list)**.5)
        else:
            if len(dir_list) > 0:
                if (min(p_list) <= 0):
                    check_indices = [i for i, x in enumerate(p_list) if x <= 0]
                    true_sum = sum([1 for i in check_indices if dir_list[i] == '+'])
                    if true_sum == len(check_indices):
                        z = inf
                    true_sum = sum([1 for i in check_indices if dir_list[i] == '-'])
                    if true_sum == len(check_indices):
                        z = -inf
                    else:
                        z = nan
            else:
                if ((min(p_list) > 0) & (max(p_list) >= 1)):
                    z = inf
                elif ((min(p_list) <= 0) & (max(p_list) < 1)):
                    z = -inf
                else:
                    z = nan
            
    else:
        z = None
    return z


def get_mapping_list(probe_key_gene_list_values):
    """ Remap a dict of {key: [values]}
    So that we get the mappings.  Here we call
    probe 1 and gene 2 in the code for simplicity
    but this function is obviously
    generic for what we are mapping

    Arguments:
     probe_key_gene_list_values: a dict of {key: [values]}

    Returns:
     mapping_list: a list with each element as
      [[1s], [2s]]
      So we can scan through the set and pick
      out the clear mapping types we like
      and can use in subsequent calculations
      
    """
    gene_set= set([])
    for the_gene_list in probe_key_gene_list_values.values():
        gene_set.update(set(the_gene_list))

    the_probe_set = set(probe_key_gene_list_values.keys())

    mapped_genes = set([])
    
    mapping_list = []

    for the_probe in the_probe_set:
        candidate_gene_set = set(probe_key_gene_list_values[the_probe])
        if len(mapped_genes.intersection(candidate_gene_set)) == 0:
            # We have a new entry and can just update the mapping_set
            mapping_list.append([set([the_probe]), candidate_gene_set])
            mapped_genes.update(candidate_gene_set)
        else:
            # a new probe but there is overlap with an existing gene set
            # means we need to merge these, can no longer distinguish
            # information from that probe
            merge_indices = [i for i, x in enumerate(mapping_list) if len(x[1].intersection(candidate_gene_set)) > 0]
            merge_indices.sort(reverse = True)
            the_new_entry = [set([the_probe]),candidate_gene_set]

            for the_index in merge_indices:
                the_new_entry[0].update(mapping_list[the_index][0])
                the_new_entry[1].update(mapping_list[the_index][1])
            for the_index in merge_indices:
                mapping_list.pop(the_index)
            mapping_list.append(the_new_entry)
            mapped_genes.update(candidate_gene_set)
    return mapping_list


def aggregate_probe_statistics_to_gene(probe_stats, probe_key_gene_list_values, **kwargs):
    """ Calculates p-values for genes based
    on the probe p-values.  Useful for mapping
    array data to values that can
    be used with a model

    Arguments:
     probe_stats: a dict with:
      {probe_id: {'p', 'dir'}}
      where p: p-value from a two-tailed test for change
            dir: '+' or '-'
     gene_key_probe_list_values: a dict with
      {probe_id: [gene_1, gene_2, ...]}

    kwargs:
     aggregation_type: if more than one probe maps
      to a gene, this is how to combine them.  Currently
      supports: ['stouffer' (default), 'max', 'min', 'none']
      One worthwhile resource:
      Fundel, K., Küffner, R., Aigner, T., & Zimmer, R. (2008). 
      Bioinformatics and Biology Insights, 2, 291–305.
      We opt not to include Fisher here since
      it there is no clear inference for directionality
      for the p-value as with these methods.  Directionality
      is useful for assessing increases/decreases in reactions.

    Returns:
     gene_mappings_dict: {mapped: {gene_id:{'p': the_p, 'dir': the_dir}} unmapped:[set(probe ids), set(gene_ids)]}
      'p': p-value for a two-tailed test
      'dir': '+' or '-'

    """

    aggregation_type = test_kwarg('aggregation_type', kwargs, ['stouffer', 'max', 'min', 'none'])

    
    # first get the mappings
    mapping_list = get_mapping_list(probe_key_gene_list_values)

    gene_mappings_dict = {}
    mapped_genes_dict = {}
    unmapped_pairs = []

    for the_pair in mapping_list:
        # We can only meaningfully compile stats
        # if we have 1 gene to 1 or more probe values
        if len(the_pair[1]) == 1:
            the_gene = list(the_pair[1])[0]
            
            if len(the_pair[0]) == 1:
                the_probe = list(the_pair[0])[0]
                the_p = probe_stats[the_probe]['p']
                the_dir = probe_stats[the_probe]['dir']
                mapped_genes_dict[the_gene] = {}
                mapped_genes_dict[the_gene]['p'] = the_p
                mapped_genes_dict[the_gene]['dir'] = the_dir                
            elif (len(the_pair[0]) > 0) & (aggregation_type != 'none'):
                the_p_list = [probe_stats[x]['p'] for x in the_pair[0]]
                the_dir_list = [probe_stats[x]['dir'] for x in the_pair[0]]                
                if aggregation_type == 'stouffer':
                    the_p_dict = stouffer_p_agg(the_p_list, dir_list = the_dir_list)
                elif aggregation_type == 'min':
                    # 'max'/'min' are the other options,
                    # but we have a good reason to prefer Stouffer's over
                    # either of these.
                    the_index = the_p_list.index(min(the_p_list))
                    the_p_dict = {'p': the_p_list[the_index], 'dir': the_dir_list[the_index]}
                else:
                    # max is the only one remaining
                    the_index = the_p_list.index(max(the_p_list))
                    the_p_dict = {'p': the_p_list[the_index], 'dir': the_dir_list[the_index]}
                mapped_genes_dict[the_gene] = the_p_dict
            else:
                unmapped_pairs.append(the_pair)
        else:
            unmapped_pairs.append(the_pair)
                
    gene_mappings_dict['mapped'] = mapped_genes_dict
    gene_mappings_dict['unmapped'] = unmapped_pairs

    return gene_mappings_dict
