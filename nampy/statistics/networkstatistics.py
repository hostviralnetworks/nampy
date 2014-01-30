# This module requires a functional rpy2 and R
        
def mtcorrect(p_value_dict, method = '"BY"'):
    """ Apply MT correction.  This is a wrapper for R's p.adjust function.

    Arguments:
     p_value_dict: a dict with upper key = condition names, next key = probe names, and values of p-values
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
     side_dict: a dict to indicate if the result lies in the upper or lower tail ('up' or 'down')
 
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
                side = 'down'
                M = M2
            else:
                side = 'up'
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
     side_dict: a dictionary of {id: 'up' or 'down'}

    Returns:
     otp_dict
     
    """
    
    otp_dict = {}
    for the_key in ttp_dict.keys():
        the_dir = side_dict[the_key]
        if the_dir == 'up':
            otp_dict[the_key] = ttp_dict[the_key] / 2.
        else:
            otp_dict[the_key] = 1. - ttp_dict[the_key] / 2.
            
    return otp_dict
