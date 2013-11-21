# Defines solvers for models.
# This script is developed from GIM3E solvers,
# Based on COBRApy solvers, which were distributed
# under GNU GENERAL PUBLIC LICENSE Version 3

# Schmidt, B., Ebrahim, A., Metz, T., Adkins, J., Palsson, B., & Hyduke, D. (2013). 
# GIM3E: condition-specific models of cellular metabolism developed from 
# metabolomics and expression data. Bioinformatics, 29(22), 2900–2908.

# Ebrahim, A., Lerman, J. a, Palsson, B. O., & Hyduke, D. R. (2013). 
# COBRApy: COnstraints-Based Reconstruction and Analysis for Python. 
# BMC systems biology, 7(1), 74. doi:10.1186/1752-0509-7-74


# Load these to help determine which solver solutions are OK
acceptable_solution_strings = ['optimal', 'MIP_optimal', 'optimal_tolerance', 'x_bound_infeasible']
# May want to include this as acceptable when running cplex
# if a small violation of the boundaries is OK: 'x_bound_infeasible'
optimal_solution_strings = ['optimal']
# Technically we can set the integer_tolerance to 0 in cplex
# This might lead to some differences in convergence from gurobi, 
# which is limited to 1e-9.
integer_tolerances = {'cplex': 0, 'gurobi': 1e-9, 'glpk': 1e-9}
from  cobra import __version__ as cobra_version

def nam_optimize(cobra_model, solver='cplex', error_reporting=True, **kwargs):
    """ A variation on optimize, constructed for network calculations.

    Arguments:
     cobra_model: the COBRA model object to import
     solver: preferred solver
     error_reporting: whether to ouput optimization
      errors (e.g. due to infeasibilities, etc...)
     kwargs: optional keyword arguments.  Note that
      because glpk requires its tolerances
      whenever it is run, updates are applied to all
      'the_solution'/lp objects for consistency
      for all solvers is performed

    """
    from cobra.flux_analysis.objective import update_objective
    from cobra import solvers
    from copy import deepcopy
    import types

    ### Solver specific parameters
    from cobra.solvers.parameters import status_dict, \
        parameter_mappings, parameter_defaults

    solver = check_solver(solver)

    # For transparency, declare parameter defaults locally
    # These are overwitten by any kwargs as configuration_parameters
    parameter_defaults = deepcopy(parameter_defaults[solver])
    parameter_defaults = {'new_objective': None,
                          'objective_sense': 'maximize',
                          'min_norm': 0,
                          'the_problem': None, 
                          'tolerance_optimality': 1e-4,
                          'tolerance_feasibility': 1e-4,
                          'tolerance_integer': 1e-9, 
                          'error_reporting': None,
                          'print_solver_time': False,
                          'quadratic_component': None}
    
    # Also, pop integer gap defaults if any
    if 'MIP_gap_abs' in parameter_defaults.keys():
        parameter_defaults.pop('MIP_gap_abs')
    elif 'MIP_gap' in parameter_defaults.keys():
        parameter_defaults.pop('MIP_gap')

    # Separate out parameters that should only be set for mips,
    # that give cplex problems if set for lps
    mip_only_parameter_defaults = {'tolerance_barrier': 1E-11}
    mip_only_parameter_dict = {}
    for the_parameter in mip_only_parameter_defaults.keys():
        if the_parameter in kwargs.keys():
            mip_only_parameter_dict[the_parameter] = kwargs[the_parameter]
            kwargs.pop(the_parameter)
        else:
            mip_only_parameter_dict[the_parameter] = mip_only_parameter_defaults[the_parameter]

    # Assume we have addressed potential solver issues in the calling function
    the_solver = solvers.solver_dict[solver]
    if 'the_problem' not in kwargs:
        kwargs['the_problem'] = None

    # Update objectives if they are new.
    # update objective will be looking reaction objects
    # or reaction integer indices as reaction objects or integers
    # Convert id's
    reaction_id_list = []
    if 'new_objective' in kwargs and \
           kwargs['new_objective'] not in ['update problem', None]:
        new_objectives = {}
        if isinstance(kwargs['new_objective'], str):
            reaction_id_list = [x.id for x in cobra_model.reactions]
            if kwargs['new_objective'] in reaction_id_list:
                the_key = cobra_model.reactions.get_by_id(kwargs['new_objective'])
                new_objectives[the_key] = 1
                kwargs['new_objective'] = new_objectives
        elif isinstance(kwargs['new_objective'], dict):
            for the_objective in kwargs['new_objective'].keys():
                if isinstance(the_objective, str):
                    if len(reaction_id_list) == 0:
                        reaction_id_list = [x.id for x in cobra_model.reactions]
                    if the_objective in reaction_id_list:
                        the_key = cobra_model.reactions.get_by_id(the_objective)
                        new_objectives[the_key] = kwargs['new_objective'][the_objective]
                else:
                    new_objectives[the_objective] = kwargs['new_objective'][the_objective]
            kwargs['new_objective'] = new_objectives
        elif isinstance(kwargs['new_objective'], list):
            for the_objective in kwargs['new_objective']:
                if isinstance(the_objective, str):
                    if len(reaction_id_list) == 0:
                        reaction_id_list = [x.id for x in cobra_model.reactions]
                    if the_objective in reaction_id_list:
                        the_key = cobra_model.reactions.get_by_id(the_objective)
                        new_objectives[the_key] = 1.
                else:
                    new_objectives[the_objective] = 1.
            kwargs['new_objective'] = new_objectives     
        update_objective(cobra_model, kwargs['new_objective'])    

    alt_cplex_flag = False
    alt_gurobi_flag = False    

    # If not given, construct the lp items manually so we can
    # supply custom settings.
    # Also note that if supplied, the_problem gets updated here
    # Note the call to update the problem just adjusts
    # the objective and bounds; any additional parameter
    # alterations are made on the call to solve_problem
    if solver == 'cplex':
        the_methods = [1, 2, 3, 4, 5, 6]
        parameter_defaults.update({'lp_method': 1,
                    'tolerance_markowitz': 0.9})
                    # It can be detrimental to set the gap too small
        configuration_parameters = deepcopy(parameter_defaults)
        configuration_parameters.update(kwargs)
        # Allow a little bit of slack in the integer solution
        # based on the selected tolerance
        # As a rule of thumb, small absolute gaps are pretty hard for the solver
        # But seem to be able to find good solutions with 1E-6.       
        if configuration_parameters['the_problem'] == None:
            the_problem = the_solver.create_problem(cobra_model, **configuration_parameters)
            #the_problem.parameters.read.scale.set(-1)
            #the_problem.parameters.emphasis.numerical.set(1)
            #the_problem.parameters.preprocessing.presolve.set(0)
            alt_cplex_flag = True
        else:
            alt_cplex_flag = True            
            the_problem = configuration_parameters['the_problem']
            the_solver.update_problem(the_problem, cobra_model, **configuration_parameters)

    elif solver == 'glpk':
        the_methods = [1, 2, 3]
        configuration_parameters = deepcopy(parameter_defaults)
        configuration_parameters.update(kwargs)        
        # Might want to adjust tm_lim        
        if configuration_parameters['the_problem'] == None:
            # Note most GLPK tolerance parameters are set at execution
            # and not lp construction
            the_problem = the_solver.create_problem(cobra_model, **configuration_parameters)
            # Additional arguments here if needed
        else:
            the_problem = configuration_parameters['the_problem']
            the_solver.update_problem(the_problem, cobra_model, **configuration_parameters)

    # Gurobi and java glpk were not tested extensively.
    # Add your own tweaks here.
    elif solver == 'gurobi':      
        the_methods = [0, 2, 1]
        parameter_defaults.update({'lp_method': 1,
                    'tolerance_markowitz': 1E-4})
        configuration_parameters = deepcopy(parameter_defaults)
        configuration_parameters.update(kwargs)
        # Allow a little bit of slack in the integer solution
        # based on the selected tolerance        
        configuration_parameters['MIP_gap_abs'] = configuration_parameters['tolerance_feasibility']
        configuration_parameters['MIP_gap'] = 0 
        if configuration_parameters['the_problem'] == None:
            the_problem = the_solver.create_problem(cobra_model, **configuration_parameters)
            #the_problem.setParam('presolve', 0)
            #the_problem.setParam('scaleflag', 0)
            alt_gurobi_flag = True            
            # Additional arguments here if needed
        else:
            the_problem = configuration_parameters['the_problem']
            alt_gurobi_flag = True            
            the_solver.update_problem(the_problem, cobra_model, **configuration_parameters)

    # didn't see this parameter in the cplex or gurobi dict
    if solver == 'cplex':
        the_problem.parameters.mip.tolerances.integrality.set(configuration_parameters['tolerance_integer'])
    if solver == 'gurobi':
        the_problem.setParam('intfeastol', configuration_parameters['tolerance_integer'])

    # Only try to set troublesome parameters for a MILP
    # in cplex or it will return an error message
    if solver == 'cplex':
        problem_type = the_problem.problem_type[the_problem.get_problem_type()]
        if problem_type != 'LP':
            the_solver.update_problem(the_problem, cobra_model, **mip_only_parameter_dict)
    else:
        the_solver.update_problem(the_problem, cobra_model, **mip_only_parameter_dict)

    ###Try to solve the problem using other methods if the first method doesn't work
    lp = the_problem
    try:
        lp_method = kwargs['lp_method']
    except:
        lp_method = 1
    if lp_method in the_methods:
        the_methods.remove(lp_method)
    #Start with the user specified method
    the_methods.insert(0, lp_method)
    # lp.set_log_stream("cplex_log.txt", lambda a:  a + " LOG " )
    # lp.set_warning_stream("cplex_warnings.txt", lambda a:  " WARNING " + a)
    # lp.set_error_stream("cplex_errors.txt", lambda a:  " ERROR " + a)
    # lp.set_results_stream("cplex_results.txt", lambda a:  " RESULT " + a)    
    acceptable_solution = None
    acceptable_solution_status = None
    for the_method in the_methods:
        configuration_parameters['lp_method'] = the_method
        try:
            if not (alt_cplex_flag | alt_gurobi_flag):
                status = the_solver.solve_problem(lp, **configuration_parameters)
            elif alt_cplex_flag:
                status = local_solve_cplex_problem(lp, **configuration_parameters)
            else:
                status = local_solve_gurobi_problem(lp, **configuration_parameters)
        except:
            status = 'failed'

        if status in optimal_solution_strings:
            break
        # Prefer clearly optimal solutions but there may be other
        # acceptable results.  Back this up and come back to if we don't get
        # an optimal one.
        elif status in acceptable_solution_strings:
            acceptable_solution = the_solver.format_solution(lp, cobra_model)
            acceptable_solution_status = status 

    if status in optimal_solution_strings:
        the_solution = the_solver.format_solution(lp, cobra_model)
    elif type(acceptable_solution) != types.NoneType:
        the_solution = acceptable_solution
        the_solution.status = status
    elif status != 'infeasible':
        if error_reporting:
            print '%s failed: %s'%(solver, status)
        the_solution = the_solver.format_solution(lp, cobra_model)
        # Keep this in case local solvers have modified the status
        the_solution.status = status
    else:
        # Otherwise we have a solution that didn't fail
        # but isn't optimal or acceptable:
        # could be infeasible
        # e.g. the solver didn't find a solution or
        # it is deemed not worthy, maybe due to quality issues
        # To be safe declare all these as infeasible
        # to force a new solution attempt
        the_solution = the_solver.format_solution(lp, cobra_model)
        the_solution.id = None
        the_solution.f = None
        the_solution.status = 'infeasible'
        if error_reporting:
            print '%s solution is infeasible'%(solver)          

    cobra_model.solution = the_solution
    return lp


def local_solve_cplex_problem(lp, **kwargs):
    """A performance tunable method for solving a problem
    

    """
    from cobra.solvers.parameters import parameter_mappings
    from cobra.solvers.cplex_solver import set_parameter, get_status

    parameter_mappings = parameter_mappings['cplex']

    # Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in kwargs.iteritems() if k in parameter_mappings]

    lp.solve()
    # Strictly enforce the bounds on the solution, if available (not available for LPs)
    if 'x_bound_error_max' in dir(lp.solution.get_quality_metrics()):
        if (lp.solution.get_quality_metrics()).x_bound_error_max <= kwargs['tolerance_feasibility']:
            status = get_status(lp)
        else:
            status = get_status(lp)
            if status not in acceptable_solution_strings:
                # Don't need to modify here
                status = status
            else:
                # Then we need to enforce x_bound_infeasible status
                # so this is not missed
                status = 'x_bound_infeasible'
                
    else:
        status = get_status(lp)
    return status


def local_solve_gurobi_problem(lp, **kwargs):
    """A performance tunable method for solving a problem
    

    """
    from cobra.solvers.parameters import parameter_mappings
    from cobra.solvers.gurobi_solver import set_parameter, get_status

    parameter_mappings = parameter_mappings['gurobi']

    # Update parameter settings if provided
    if kwargs:
        [set_parameter(lp, parameter_mappings[k], v)
         for k, v in kwargs.iteritems() if k in parameter_mappings]
    
    lp.optimize()
    status = get_status(lp)
    return status

def check_solver(solver = 'cplex'):
    """ A simple function to verify the solver to employ.


    """
    from cobra import solvers
    
    available_solvers = solvers.solver_dict.keys()

    if solver not in available_solvers:
        if 'cplex' in available_solvers:
            print("Switched solver from " + solver + " to cplex.")            
            return 'cplex'
        # Gurobi should work well but it is not
        # tested as thoroughly as cplex
        elif 'gurobi' in available_solvers:
            print("Switched solver from " + solver + " to gurobi.")             
            return 'gurobi'
        # glpk works well for LPs but is usually not effective
        # for MILPs
        elif 'glpk' in available_solvers:
            print("Switched solver from " + solver + " to glpk.") 
            return 'glpk'
        # Otherwise, none available.
        # Also note Java_glpk not tested/supported
        else:
            print("No working solvers detected.")
            return None
    else:
        return solver


def get_solver_parameters(**kwargs):
    """ A simple function to return the solver tolerance parameters
    Note the tolerances will preferably be supplied, the defaults
    may not be consistent with those used previously

    """
    if 'solver' in kwargs:
        solver = check_solver(kwargs['solver'])
    else:
        solver = check_solver('glpk')

    if 'tolerance_optimality' not in kwargs:
        tolerance_optimality = 1E-7
    else:
        tolerance_optimality = kwargs['tolerance_optimality']    

    if 'tolerance_feasibility' not in kwargs:
        tolerance_feasibility = 1E-7
    else:
        tolerance_feasibility = kwargs['tolerance_feasibility']

    if 'tolerance_barrier' not in kwargs:
        tolerance_barrier = 1E-7 * 0.0001
    else:
        tolerance_barrier = kwargs['tolerance_barrier']

    if 'tolerance_integer' not in kwargs:
        tolerance_integer = integer_tolerances[solver]
    else:
        tolerance_integer = kwargs['tolerance_integer']
    
    return solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer



def convert_objective_to_constraint(cobra_model, **kwargs):
    """ This function converts the objective to a turnover metabolite and
    reaction, applying the appropriate constraints.  Note if you
    wish to leave the target optimum unbounded and
    only restrict the minimum you should manipulate
    the appropriate bound.

    Arguments:
     cobra_model

    kwargs:

     objective_sense: 'maximize' or 'minimize'
      will also affect how fraction of optimum is applied

     fraction_of_optimum: fraction of the best achievable value to restrict
      the solution to.  Restricted to [0, 1].

     copy_model: whether to modify the model in place
      or return a copy.  Booelean: True/False
     
     new_reaction_name: what to call the new objective

     bound_best_optimum: a bound will always be placed at the worst
      objective value, as determined by fraction_of_optimum.
      This Boolean variable indicates whether to also add
      one to reflect the best.  True cements the best
      bound into the model and is a handy reference.
      If set to False, and the implicit bound (0 or 1000)
      returned when setting up the reaction is
      greater, then the bound will not be further constrained.

     solver parameters: solver, tolerance_feasibility,  tolerance_optimality,  tolerance_integer

    
    """
    from copy import deepcopy
    from cobra.core.Metabolite import Metabolite
    from cobra.core.Reaction import Reaction
    from numpy import array
    from math import ceil
    from math import floor

    if 'objective_sense' in kwargs:
        objective_sense = kwargs['objective_sense']
    else:
        objective_sense = 'maximize'

    if 'fraction_of_optimum' in kwargs:
        fraction_of_optimum = kwargs['fraction_of_optimum']
    else:
        fraction_of_optimum = 1.

    if 'copy_model' in kwargs:
        copy_model = kwargs['copy_model']
    else:
        copy_model = True

    if 'new_reaction_name' in kwargs:
        new_reaction_name = kwargs['new_reaction_name']
    else:
        new_reaction_name = "objective"

    if 'bound_best_optimum' in kwargs:
        bound_best_optimum = kwargs['bound_best_optimum']
    else:
        bound_best_optimum = True

    solver, tolerance_optimality, tolerance_feasibility, tolerance_barrier, tolerance_integer = get_solver_parameters(**kwargs)
    
    continue_conversion = False

    if copy_model:
        if cobra_version == '0.2.0':
            cobra_model = deepcopy(cobra_model)
        else:
            cobra_model = cobra_model.copy()

    nam_optimize(cobra_model, objective_sense = objective_sense, solver = solver,
        tolerance_optimality = tolerance_optimality,
        tolerance_feasibility = tolerance_feasibility,
        tolerance_barrier = tolerance_barrier,
        tolerance_integer = tolerance_integer)
    wt_solution = cobra_model.solution.f

    # Just pick a new name in case we have multiple
    # Calls to this function.
    named_objective_constraint = False
    while named_objective_constraint == False:
        same_name_list = cobra_model.metabolites.query(new_reaction_name)
        if len(same_name_list) < 1:
            named_objective_constraint = True
        else:
            new_reaction_name = new_reaction_name + "_1"

    objective_metabolite = Metabolite(new_reaction_name)
    objective_reaction = Reaction(new_reaction_name)
    old_obj_list = []
    old_obj_coefficients = []
    for x in cobra_model.reactions:
        if x.objective_coefficient != 0:
            x.add_metabolites({objective_metabolite: x.objective_coefficient})
            old_obj_list.append(x)
            old_obj_coefficients.append(x.objective_coefficient)
            x.objective_coefficient = 0
            continue_conversion = True
    if continue_conversion == False:
        print("No objective detected, exiting conversion of objective to constraint...")
    else:
        objective_reaction.add_metabolites({objective_metabolite:-1})
        cobra_model.add_reactions(objective_reaction)
        if objective_sense == 'maximize':
            # Make sure the new objective isn't overly constrained
            if bound_best_optimum | (objective_reaction.upper_bound < wt_solution):
                objective_reaction.upper_bound = wt_solution            
            if wt_solution >= 0:
                objective_reaction.lower_bound = wt_solution * fraction_of_optimum
            else:
                if fraction_of_optimum > 0:
                    # See the "minimize" case for comments.
                    objective_reaction.lower_bound = wt_solution / fraction_of_optimum
                else:
                    objective_reaction.lower_bound =\
                    floor(sum(array([(old_obj_coefficients[index] * max(0, x.upper_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] < 0)])) +\
                    sum(array([(old_obj_coefficients[index] * min(0, x.lower_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] > 0)])))   
        elif objective_sense == 'minimize':
            if bound_best_optimum | (objective_reaction.lower_bound > wt_solution):
                objective_reaction.lower_bound = wt_solution            
            if wt_solution >= 0:
                # There are several options to scale the upper bound  for poor 
                # consistency in this scenario.  We could find the global worst 
                # case and scale over this range, but then the scaling becomes 
                # somewhat arbitrary.  We also could scale by 1 + (1 - 
                # fr_optimum), but then to explore ranges much worse than the 
                # optimum we need negative fr_optimum values.  A ratio is
                # therefore applied in the case fr_optimum is between 0 and 1.
                if fraction_of_optimum > 0:
                    objective_reaction.upper_bound = wt_solution / fraction_of_optimum
                # Otherwise, the biggest value is determined by the upper_bounds
                # (to avoid divide by zero errors)
                else:
                    objective_reaction.upper_bound =\
                    ceil(sum(array([(old_obj_coefficients[index] * max(0, x.upper_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] > 0)])) +\
                    sum(array([(old_obj_coefficients[index] * min(0, x.lower_bound))
                               for index, x in enumerate(old_obj_list) if
                               (old_obj_coefficients[index] < 0)])))    
            else:
                objective_reaction.upper_bound = wt_solution * fraction_of_optimum
        objective_reaction.objective_coefficient = 1.
        # run one optimization to store the best value
        # in cobra_model.solution
        nam_optimize(cobra_model, solver=solver,
            objective_sense=objective_sense,
            # Not supplying a problem and
            # have already updated the
            # model objective coefficients
            # e.g., no new_objective = {objective_reaction: 1.},
            tolerance_optimality = tolerance_optimality,
            tolerance_feasibility = tolerance_feasibility,
            tolerance_barrier = tolerance_barrier,
            tolerance_integer = tolerance_integer)
        objective_reaction.objective_coefficient = 0  

    return cobra_model
