
import itertools
import random
import math
from itertools import chain
import time

from optlang.interface import OPTIMAL

from cobra.flux_analysis.gapfilling import GapFiller
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.core import DictList

from cobra.util.solver import linear_reaction_coefficients

from medusa.core.ensemble import Ensemble
from medusa.core.feature import Feature
from medusa.core.member import Member

# functions for expanding existing models to generate an ensemble

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

def gapfill_to_ensemble(model, iterations=1, universal=None, lower_bound=0.05,
                 penalties=None, exchange_reactions=False,
                 demand_reactions=False, integer_threshold=1e-6):
    """
    Performs gapfilling on model, pulling reactions from universal.
    Any existing constraints on base_model are maintained during gapfilling, so
    these should be set before calling gapfill_to_ensemble (e.g. secretion of
    metabolites, choice of objective function etc.).

    Currently, only iterative solutions are supported with accumulating
    penalties (i.e. after each iteration, the penalty for each reaction
    doubles).

    Parameters
    ----------
    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model
        A universal model with reactions that can be used to complete the
        model.
    lower_bound : float, 0.05
        The minimally accepted flux for the objective in the filled model.
    penalties : dict, None
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively.
    integer_threshold : float, 1e-6
        The threshold at which a value is considered non-zero (aka
        integrality threshold). If gapfilled models fail to validate,
        you may want to lower this value. However, picking a threshold that is
        too low may also result in reactions being added that are not essential
        to meet the imposed constraints.
    exchange_reactions : bool, False
        Consider adding exchange (uptake) reactions for all metabolites
        in the model.
    demand_reactions : bool, False
        Consider adding demand reactions for all metabolites.

    Returns
    -------
    ensemble : medusa.core.Ensemble
        The ensemble object created from the gapfill solutions.
    """
    # copy original objective to reassign after generating ensemble
    original_coefficients = \
            model.objective.get_linear_coefficients(model.variables).copy()
    gapfiller = GapFiller(model, universal=universal,
                          lower_bound=lower_bound, penalties=penalties,
                          demand_reactions=demand_reactions,
                          exchange_reactions=exchange_reactions,
                          integer_threshold=integer_threshold)
    solutions = gapfiller.fill(iterations=iterations)
    print("finished gap-filling. Constructing ensemble...")
    ensemble = _build_ensemble_from_gapfill_solutions(model,solutions,
                                                    universal=universal)
    ensemble.base_model.objective.set_linear_coefficients(original_coefficients)
    return ensemble

def iterative_gapfill_from_binary_phenotypes(model,universal,phenotype_dict,
                                            output_ensemble_size,
                                            gapfill_type="continuous",
                                            iterations_per_condition=1,
                                            lower_bound=0.05,
                                            remove_fraction = 0.0,
                                            remove_exceptions = [],
                                            penalties=None,
                                            exchange_reactions=False,
                                            demand_reactions=False,
                                            inclusion_threshold=1e-6,
                                            exchange_prefix="EX_",
                                            solver='glpk'):
    """
    Performs gapfilling on model, pulling reactions from universal.
    Any existing constraints on base_model are maintained during gapfilling, so
    these should be set before calling gapfill_to_ensemble (e.g. secretion of
    metabolites, choice of objective function etc.).

    Cycles through each key:value pair in phenotype_dict, iterating over every
    condition and performing gapfilling on that condition until the number of
    cycles over all conditions is equal to output_ensemble_size. For each cycle,
    the order of conditions is randomized, which generally leads to unique sets
    of solutions for each cycle.

    Currently only supports a single iteration for each condition within each
    cycle (i.e. for each gapfill in a single condition, only one solution is
    returned). Currently only supports gapfilling to positive growth
    conditions.

    Generally, solutions are easier to find and more likely to exist if users
    ensure that transporters for metabolites exist in the model
    already or at minimum are present in the universal model.

    Parameters
    ----------
    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model
        A universal model with reactions that can be used to complete the
        model.
    phenotype_dict : dict
        A dictionary of condition_name:media_dict, where condition name is a
        unique string describing the condition (such as the name of a single
        carbon source) and media_dict is a dictionary of exchange reactions to
        bounds, as set in cobra.core.model.medium. Exchange reactions are
        provided with reaction.id, not cobra.core.reaction objects.
    output_ensemble_size : int
        Number of cycles over all conditions provided to perform. Equal to the
        number of lists returned in 'solutions'. When the ensemble is
        constructed, the number of members may be lower than
        output_ensemble_size if any duplicate solutions were found across
        cycles.
    gapfill_type : string, "continuous"
        The kind of gapfilling to perform. Only "continuous" is currently
        implemented. Previous version allowed "integer", but this is
        generally much too slow to be practical (but may be reimplemented
        in the future given a compelling use case).
    iterations_per_condition : int, 1
        The number of gapfill solutions to return in each condition within each
        cycle. Currently only supports returning a single solution.
    lower_bound : float, 0.05
        The minimally accepted flux for the objective in the filled model.
    penalties : dict, None
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively.
    remove_fraction : float, 0.0
        The fraction of reactions to remove from the original model prior
        to each gapfilling cycle. The ensemble member generated by each cycle
        will have a random set of reactions removed prior to gapfilling, i.e.
        the reactions removed before gapfilling are NOT the same for all
        members.
    remove_exceptions : list, []
        The reactions to exclude from random reaction removal. The fraction of
        reactions removed is calculated on the set of reactions remaining when
        not considering these reactions.
    inclusion_threshold : float, 1e-6
        The threshold at which a value is considered non-zero (aka
        integrality threshold in the integer formulation, or the flux threshold
        in the continuous formulation). If gapfilled models fail to validate,
        you may want to lower this valu. However, picking a threshold that is
        too low may also result in reactions being added that are not essential
        to meet the imposed constraints.
    exchange_reactions : bool, False
        Consider adding exchange (uptake) reactions for all metabolites
        in the model.
    demand_reactions : bool, False
        Consider adding demand reactions for all metabolites.
    exchange_prefix : string, "EX_"
        the default reaction ID prefix to search for when identifying exchange
        reactions. "EX_" is standard for modelSEED models. This will be
        updated to be more database-agnostic when cobrapy boundary
        determination is finalized for cobrapy version 1.0.
    solver : string, "glpk"
        The solver to use when gapfilling. Any cobrapy-compatible solver that
        is set up to function with cobrapy will work.

    Returns
    -------
    solutions : list
        list of lists; each list contains a gapfill solution for a single
        cycle. Number of lists is equal to output_ensemble_size.
    """
    if gapfill_type not in ["continuous"]:
        raise ValueError("only gapfill types of continuous"
                         "are supported")

    # Check that all exchange reactions exist in base model. If not, raise an
    # error.
    exchanges_in_phenotype_dict = set()
    for condition in phenotype_dict.keys():
        exchanges_in_phenotype_dict = exchanges_in_phenotype_dict | set(
                                        phenotype_dict[condition].keys())

    model_rxns = [rxn.id for rxn in model.reactions]
    for exchange in exchanges_in_phenotype_dict:
        if not exchange in model_rxns:
            print('Could not find '+ exchange)
            raise ValueError('phenotype_dict contains exchange reactions not '
                            'found in model. Add any necessary exchange '
                            'reactions prior to gapfilling.')


    # pre-generate the random orderings of conditions to ensure there are no
    # duplicates.
    cycle_order = []
    while len(cycle_order) < output_ensemble_size:
        condition_list = random.sample(list(phenotype_dict.keys()),
                                        len(phenotype_dict.keys()))
        if condition_list not in cycle_order:
            cycle_order.append(condition_list)

    # Cycle through all conditions and gapfill. After each gapfill iteration,
    # our strategy is to reduce the cost for the reactions returned by the
    # previous solution to 0, such that they are automatically included in
    # the model for the next condition.
    if gapfill_type is "continuous":
        solutions, maintained_reactions, gapfiller = _continuous_iterative_binary_gapfill(model,
                              phenotype_dict,
                              cycle_order,
                              universal=universal,
                              output_ensemble_size=output_ensemble_size,
                              lower_bound=lower_bound,
                              penalties=penalties,
                              remove_fraction=remove_fraction,
                              remove_exceptions=remove_exceptions,
                              demand_reactions=demand_reactions,
                              exchange_reactions=exchange_reactions,
                              flux_cutoff=inclusion_threshold,
                              exchange_prefix='EX_',
                              solver=solver)

    ensemble =_build_ensemble_from_gapfill_solutions(model,
                                                    solutions,
                                                    maintained_reactions,
                                                    universal=gapfiller)
    return ensemble

def _continuous_iterative_binary_gapfill(model,phenotype_dict,cycle_order,
                      universal=None, output_ensemble_size=1,
                      lower_bound=0.05, penalties=None,
                      remove_fraction = 0.0,
                      remove_exceptions = [],
                      demand_reactions=False,
                      exchange_reactions=False,
                      flux_cutoff=1E-8,
                      exchange_prefix='EX_',
                      solver="glpk"):

    if exchange_reactions:
        raise NotImplementedError("Inclusion of new exchange reactions is not"
                            "supported for continuous gapfill")
    if demand_reactions:
        raise NotImplementedError("Inclusion of demand reactions is not"
                            "supported for continuous gapfill")

    solutions = [] # list of lists for gapfill solutions
    maintained_reactions = [] # list of lists for reactions maintained in model
    gapfiller = universal.copy()

    # set the solver for the gapfill object
    gapfiller.solver = solver

    # get the original objective from the model being gapfilled
    model_to_gapfill = model.copy()
    original_objective = linear_reaction_coefficients(model_to_gapfill)
    # convert to IDs to avoid issues with model membership when these reactions
    # are added to gapfiller
    original_objective = {rxn.id:original_objective[rxn] for rxn
                            in original_objective.keys()}

    # remove reactions from the gapfiller that are present in the model
    # to be gapfilled
    gapfiller.remove_reactions([rxn for rxn in gapfiller.reactions if rxn.id in \
                        [rxn.id for rxn in model_to_gapfill.reactions]])

    # add the reactions from the model to the gapfiller
    # Get IDs for later use
    original_model_reaction_ids = [reaction.id for reaction
                                in model_to_gapfill.reactions]
    gapfiller.add_reactions([rxn.copy() for rxn
                                in model_to_gapfill.reactions])

    # Add the pFBA constraints and objective (minimizes sum of fluxes)
    add_pfba(gapfiller)

    # set a constraint on flux through the original objective
    for reaction in original_objective.keys():
        print("Constraining lower bound for " + reaction)
        gapfiller.reactions.get_by_id(reaction).lower_bound = lower_bound

    # Close the lower bound on exchange reactions
    exchange_reactions = [rxn for rxn in gapfiller.reactions if\
                            rxn.id.startswith(exchange_prefix)]
    for rxn in exchange_reactions:
        rxn.bounds = (0,rxn.upper_bound)

    original_coefficients = \
            gapfiller.objective.get_linear_coefficients(gapfiller.variables).copy()

    for cycle_num in range(0,output_ensemble_size):
        print("starting cycle number " + str(cycle_num))
        cycle_reactions = set()

        if remove_fraction > 0.0:
            # sample the original reaction id's and select fraction of them,
            # removing reactions in remove_exceptions and the old objective
            removable = list(set(original_model_reaction_ids) -
                             set(remove_exceptions) -
                             set(original_objective.keys()))

            original_to_remove = random.sample(
                    removable,
                    math.floor(
                    (1.0 - remove_fraction)*len(removable)))
            # add the remove_exceptions and objective reactions back to
            # original_to_remove so that the penalties on them get removed.
            original_to_remove.extend(remove_exceptions)
            original_to_remove.extend(original_objective.keys())
        else:
            #otherwise, remove the pentalty for all reactions in the model
            original_to_remove = original_model_reaction_ids

        # Remove penalty
        _set_coefficients(gapfiller, original_to_remove, 0.0)
        gapfiller.slim_optimize()
        # get reaction ids we need to check for flux after optimizing
        get_fluxes = list(set([r.id for r in gapfiller.reactions]) -
                                set(original_to_remove))

        for condition in cycle_order[cycle_num]:
            # set the medium for this condition.
            _set_medium(gapfiller, phenotype_dict[condition], exchange_prefix)

            # gapfill and get the solution
            iteration_solution = gapfiller.optimize()
            filtered_solution = {rxn:iteration_solution.fluxes[rxn] for rxn in\
               get_fluxes if abs(iteration_solution.fluxes[rxn]) > flux_cutoff}

            # combine solution with those from previous iterations within cycle
            cycle_reactions = cycle_reactions | \
                            set([rxn for rxn in filtered_solution.keys()])
            # Get the reaction objects to be added for model validation
            validate_reactions = [gapfiller.reactions.get_by_id(r).copy() for r in
                                    cycle_reactions]

            # set media conditions for the model to be gapfilled
            _set_medium(model_to_gapfill, phenotype_dict[condition], exchange_prefix)
            if not validate(model_to_gapfill,
                            validate_reactions,
                            original_to_remove,
                            lower_bound):
                raise RuntimeError('Failed to validate gapfilled model, '
                                    'try lowering the flux_cutoff through '
                                    'inclusion_threshold')
            # remove the flux minimization penalty on the gapfilled reactions
            _set_coefficients(gapfiller, filtered_solution.keys(), 0.0)
            #check = gapfiller.slim_optimize() # optimizing might be necessary
            # to update coefficients.

            # reset the media condition
            for ex_rxn in exchange_reactions:
                ex_rxn.lower_bound = 0

        gapfiller.objective.set_linear_coefficients(original_coefficients)
        solutions.append(list(cycle_reactions))
        maintained_reactions.append(original_to_remove)
    return solutions, maintained_reactions, gapfiller


def _integer_iterative_binary_gapfill(model,phenotype_dict,cycle_order,
                      universal=None, output_ensemble_size=0,
                      lower_bound=0.05, penalties=False,
                      remove_fraction = None,
                      demand_reactions=False,
                      exchange_reactions=False,
                      integer_threshold=1E-6,
                      solver="glpk"):

    if remove_fraction:
        raise NotImplementedError("Random removal of reactions is not"
                            "supported for integer gapfill")
    gapfiller = GapFiller(model, universal=universal,
                          lower_bound=lower_bound, penalties=penalties,
                          demand_reactions=demand_reactions,
                          exchange_reactions=exchange_reactions,
                          integer_threshold=integer_threshold)
    original_costs = gapfiller.costs

    solutions = []
    for cycle_num in range(0,output_ensemble_size):
        print("starting cycle number " + str(cycle_num))
        cycle_reactions = set()
        for condition in cycle_order[cycle_num]:
            gapfiller.model.medium = phenotype_dict[condition]
            # gapfill and get the solution. The 0 index is necessary because
            # gapfill will return a list of lists; we are only taking the
            # first (and only) list here.
            gapfilled_reactions = gapfiller.fill()[0]
            cycle_reactions = cycle_reactions | set(gapfilled_reactions)
            # iterate through indicators, find those corresponding to the
            # gapfilled reactions from any iteration within this cycle, and
            # reset their cost to 0. Doing this for all reactions from any
            # iteration within the cycle is necessary because cobrapy's
            # gapfill function performs update_costs, which will reset costs
            # and iteratively increase them; without this manual override
            # performed here, costs for previous conditions within a cycle
            # would revert to 1 instead of the desired 0
            for reaction_indicator in [indicator for indicator
                                        in gapfiller.indicators]:
                if reaction_indicator.rxn_id in cycle_reactions:
                    gapfiller.costs[reaction_indicator] = 0
            gapfiller.model.objective.set_linear_coefficients(gapfiller.costs)
        solutions.append(list(cycle_reactions))
        gapfiller.model.objective.set_linear_coefficients(original_costs)
    return solutions


def _build_ensemble_from_gapfill_solutions(model,
                                           solutions,
                                           maintained_reactions,
                                           universal=None):

    ensemble = Ensemble(identifier=model.id,name=model.name)
    ensemble.base_model = model.copy()

    # generate member identifiers for each solution
    # Convert the solution to a dictlist so we can retrieve reactions by id
    solution_dict = {}
    i = 0
    for solution in solutions:
        solution_id = model.id + '_gapfilled_' + str(i)

        # Add the maintained reactions to the solution.
        solution.extend(maintained_reactions[i])

        # For any reaction in the solution missing from universal,
        # remove it from the model and add it to universal. This is done for
        # reactions that were part of a removal fraction, allowing us to call
        # universal.reactions.get_by_id() rather than constantly checking
        # whether a reaction came fom the model or the universal.
        add_to_universal = []
        remove_from_model = []
        solutions_as_rxn_objs = []
        for rxn_id in solution:
            if rxn_id in [r.id for r in model.reactions]:
                add_to_universal.append(model.reactions.get_by_id(rxn_id).copy())
                remove_from_model.append(rxn_id)
            else:
                solutions_as_rxn_objs.append(universal.reactions.get_by_id(rxn_id).copy())

        universal.add_reactions([rxn for rxn in add_to_universal
                                if rxn.id not in
                                    [r.id for r in universal.reactions]])
        ensemble.base_model.remove_reactions(remove_from_model, remove_orphans = True)

        solutions_as_rxn_objs.extend(add_to_universal)
        solution_dict[solution_id] = DictList() + solutions_as_rxn_objs
        i += 1

    # scan through other members and remove them if they are identical.
    # as long as we're looping, we'll also find reactions that are in
    # all solutions and get the list of reactions that are in any ensemble.

    # first, convert the solution dictionary to the same structure except with
    # reaction ids rather than reaction objects so that we can perform set
    # operations
    solutions_as_ids = {}
    for member_id in solution_dict.keys():
        solutions_as_ids[member_id] = [rxn.id for
                                        rxn in solution_dict[member_id]]

    used_members = []
    duplicate_solutions = []
    all_reactions = set()
    in_all = set()
    for member_id in solutions_as_ids.keys():
        used_members.append(member_id)
        member_solution = set(solutions_as_ids[member_id])
        if in_all:
            in_all = member_solution & in_all
        #if this is the first ensemble member, set intersection will fail
        # because of the empty set, so we need this exception
        else:
            in_all = member_solution
        all_reactions = all_reactions | member_solution
        for other_member in solutions_as_ids.keys():
            if other_member not in used_members:
                other_solution = solutions_as_ids[other_member]
                if set(other_solution) == member_solution:
                    duplicate_solutions.append(other_member)
                    used_members.append(other_member)

    # perform the removal of duplicate solutions on the original solution
    # object which contains reaction objects rather than reaction ids as
    # strings
    for duplicate in duplicate_solutions:
        solution_dict.pop(duplicate,None)

    # Reactions that need features are those that were not in all the gapfill
    # solutions.
    reactions_needing_features = list(all_reactions - in_all)
    reactions_needing_features_objs = [
                    universal.reactions.get_by_id(rxn).copy()
                    for rxn in reactions_needing_features]

    # add reaction objects to the base model for all reactions
    all_reactions_as_objects = [universal.reactions.get_by_id(rxn).copy()
                                for rxn in all_reactions]
    all_reactions_to_add = [rxn for rxn in all_reactions_as_objects
                            if rxn.id not in
                                [r.id for r in ensemble.base_model.reactions]]
    ensemble.base_model.add_reactions(all_reactions_as_objects)

    # add metabolite objects to the base model for all new metabolites from
    # the new reactions -- NOT NECESSARY WITH NEW COBRAPY VERSIONS
    # mets = [x.metabolites for x in all_reactions_as_objects]
    # all_keys = set().union(*(d.keys() for d in mets))
    # ensemble.base_model.add_metabolites(all_keys)

    print('building features...')
    # generate features for the reactions that vary across solutions and add
    # them to the ensemble. assume that all reactions have the same attribute
    # values; if different attribute values are desired for reactions with the
    # same ID, these need to be added to the universal reaction bag prior to
    # gapfilling
    features = DictList()
    for reaction in reactions_needing_features_objs:
        for attribute in REACTION_ATTRIBUTES:
            identifier = reaction.id + "_" + attribute
            name = reaction.name
            feature = Feature(identifier,name)
            feature.ensemble = ensemble
            feature.base_component = (ensemble.base_model.reactions
                                        .get_by_id(reaction.id))
            feature.component_attribute = attribute
            # get the states for this feature as {member.id:value}
            states = {}
            for member_id in solution_dict.keys():
                # get the reaction and its value for the attribute
                if reaction.id in [rxn.id for rxn in solution_dict[member_id]]:
                    rxn_obj = solution_dict[member_id].get_by_id(reaction.id)
                    states[member_id] = getattr(rxn_obj,attribute)
                else:
                    states[member_id] = MISSING_ATTRIBUTE_DEFAULT[attribute]
            feature.states = states
            features += [feature]

    ensemble.features = features

    print('updating members...')
    # update members for the ensemble
    members = DictList()
    for member_id in solution_dict.keys():
        model_states = dict()
        for feature in ensemble.features:
            model_states[feature] = feature.get_model_state(member_id)
        member = Member(ensemble=ensemble,\
                        identifier=member_id,\
                        name=ensemble.name,\
                        states=model_states)
        members += [member]

    ensemble.members = members

    return ensemble

def validate(original_model, reactions, original_to_remove, lower_bound):
    with original_model as model:
        # get any reactions that were inactivated randomly
        # NOTE: original_to_remove are the original reactions that we
        # remove the gapfill penalty for, NOT those that were removed
        # from the model.

        keep = set([r.id for r in model.reactions]) & set(set(original_to_remove) | set([r.id for r in reactions]))
        remove = set([r.id for r in model.reactions]) - keep
        add = set([r.id for r in reactions]) - keep

        model.remove_reactions(remove)
        model.add_reactions([r for r in reactions if r.id in add])
        model.slim_optimize()
        return (model.solver.status == OPTIMAL and
                model.solver.objective.value >= lower_bound)

def _set_coefficients(model, reactions, coefficient):
    """
    Helper function to set the coefficients for a list of reactions.

    Modifies the model in place.

    Parameters
    ----------
    model : cobra.Model
        model that the coefficients should be modified for (e.g. the
        model being gapfilled).
    reactions : list
        list of reaction ids to adjust coefficients for.
    coefficient : float
        Value to set the coefficient to.
    """

    coefficients = (model.objective.
                get_linear_coefficients(model.variables).copy())
    reaction_variables = (((model.reactions.get_by_id(rxn)
                            .forward_variable),
                            (model.reactions.get_by_id(rxn)
                            .reverse_variable))
                             for rxn in reactions)
    variables = chain(*reaction_variables)
    for variable in variables:
         coefficients[variable] = 0.0
    model.objective.set_linear_coefficients(coefficients)

def _set_medium(model, condition, exchange_prefix):
    """
    Helper function to set the media conditions. Slightly faster
    than using model.medium.

    Parameters
    ----------
    model : cobra.Model
        Model to adjust the medium for.
    condition : dict
        Dictionary mapping exchange reactions to the absolute value
        of their uptake rate.
    exchange_prefix : string
        The prefix to search reaction ids for when closing
        other exchange reactions.
    """
    for ex_rxn in [rxn for rxn in model.reactions if \
                    rxn.id.startswith(exchange_prefix)]:
        ex_rxn.bounds = (0,ex_rxn.upper_bound)
    for ex_rxn in condition.keys():
        model.reactions.get_by_id(ex_rxn).bounds = \
            (-1.0*condition[ex_rxn],
            1.0*condition[ex_rxn])
