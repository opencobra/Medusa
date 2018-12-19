
import itertools
import random
from itertools import chain

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
    gapfiller = GapFiller(model, universal=universal,
                          lower_bound=lower_bound, penalties=penalties,
                          demand_reactions=demand_reactions,
                          exchange_reactions=exchange_reactions,
                          integer_threshold=integer_threshold)
    solutions = gapfiller.fill(iterations=iterations)
    print("finished gap-filling. Constructing ensemble...")
    ensemble = _build_ensemble_from_gapfill_solutions(model,solutions,
                                                    universal=universal)

    return ensemble

def iterative_gapfill_from_binary_phenotypes(model,universal,phenotype_dict,
                                            output_ensemble_size,
                                            gapfill_type="continuous",
                                            iterations_per_condition=1,
                                            lower_bound=0.05, penalties=None,
                                            exchange_reactions=False,
                                            demand_reactions=False,
                                            inclusion_threshold=1e-6,
                                            exchange_prefix="EX_"):
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

    Returns
    -------
    solutions : list
        list of lists; each list contains a gapfill solution for a single
        cycle. Number of lists is equal to output_ensemble_size.
    """
    if gapfill_type not in ["integer","continuous"]:
        raise ValueError("only gapfill types of integer and continuous"
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
    if gapfill_type is "integer":
        solutions =  _integer_iterative_binary_gapfill(model,
                              phenotype_dict,
                              cycle_order,
                              universal=universal,
                              output_ensemble_size=output_ensemble_size,
                              lower_bound=lower_bound,
                              penalties=penalties,
                              demand_reactions=demand_reactions,
                              exchange_reactions=exchange_reactions,
                              integer_threshold=inclusion_threshold)
    elif gapfill_type is "continuous":
        solutions = _continuous_iterative_binary_gapfill(model,
                              phenotype_dict,
                              cycle_order,
                              universal=universal,
                              output_ensemble_size=output_ensemble_size,
                              lower_bound=lower_bound,
                              penalties=penalties,
                              demand_reactions=demand_reactions,
                              exchange_reactions=exchange_reactions,
                              flux_cutoff=inclusion_threshold,
                              exchange_prefix='EX_')

    ensemble =_build_ensemble_from_gapfill_solutions(model,solutions,
                                                    universal=universal)
    return ensemble

def _continuous_iterative_binary_gapfill(model,phenotype_dict,cycle_order,
                      universal=None, output_ensemble_size=1,
                      lower_bound=0.05, penalties=None,
                      demand_reactions=False,
                      exchange_reactions=False,
                      flux_cutoff=1E-8,
                      exchange_prefix='EX_'):
    if exchange_reactions:
        raise NotImplementedError("Inclusion of new exchange reactions is not"
                            "supported for continuous gapfill")
    if demand_reactions:
        raise NotImplementedError("Inclusion of demand reactions is not"
                            "supported for continuous gapfill")

    solutions = []
    gapfiller = universal.copy()



    # get the original objective from the model being gapfilled
    model_to_gapfill = model.copy()
    original_objective = linear_reaction_coefficients(model_to_gapfill)
    # convert to IDs to avoid issues with model membership when these reactions
    # are added to gapfiller
    original_objective = {rxn.id:original_objective[rxn] for rxn
                            in original_objective.keys()}

    # get the reactions in the original model, which need to be removed from
    # the universal if present. This cannot catch identical reactions that do
    # not share IDs, so make sure your model and universal are in the same
    # namespace.
    rxns_to_remove = [rxn for rxn in gapfiller.reactions if rxn.id in \
                        [rxn.id for rxn in model_to_gapfill.reactions]]
    gapfiller.remove_reactions(rxns_to_remove)

    # get the list of reactions currently in the gapfiller, which are the ones
    # we will need to check for flux after solving the problem (e.g. these are
    # the reactions we are considering adding to the model)
    get_fluxes = [rxn.id for rxn in gapfiller.reactions]

    # add the reactions from the model to the gapfiller, which are not
    # included in the pFBA formulation, and thus flux is not penalized
    # through them.
    original_model_reactions = [rxn.copy() for rxn
                                in model_to_gapfill.reactions]
    gapfiller.add_reactions(original_model_reactions)
    original_reaction_ids = [reaction.id for reaction
                                in original_model_reactions]


    # Add the pFBA constraints and objective (minimizes sum of fluxes)
    add_pfba(gapfiller)

    # set the linear coefficients for reactions in the original model to 0
    coefficients = (gapfiller.objective
                    .get_linear_coefficients(gapfiller.variables))
    reaction_variables = (((gapfiller.reactions.get_by_id(reaction)
                            .forward_variable),
                           (gapfiller.reactions.get_by_id(reaction)
                            .reverse_variable))
                            for reaction in original_reaction_ids)
    variables = chain(*reaction_variables)
    for variable in variables:
        coefficients[variable] = 0.0
    gapfiller.objective.set_linear_coefficients(coefficients)

    ## set a constraint on flux through the original objective
    for reaction in original_objective.keys():
        print("Constraining lower bound for " + reaction)
        gapfiller.reactions.get_by_id(reaction).lower_bound = lower_bound

    exchange_reactions = [rxn for rxn in gapfiller.reactions if\
                            rxn.id.startswith(exchange_prefix)]
    for rxn in exchange_reactions:
        rxn.lower_bound = 0

    for cycle_num in range(0,output_ensemble_size):
        print("starting cycle number " + str(cycle_num))
        cycle_reactions = set()
        original_coefficients = \
            gapfiller.objective.get_linear_coefficients(gapfiller.variables)

        for condition in cycle_order[cycle_num]:
            # set the medium for this condition.
            for ex_rxn in phenotype_dict[condition].keys():
                gapfiller.reactions.get_by_id(ex_rxn).lower_bound = \
                    -1.0*phenotype_dict[condition][ex_rxn]
                gapfiller.reactions.get_by_id(ex_rxn).upper_bound = \
                    1.0*phenotype_dict[condition][ex_rxn]

            # gapfill and get the solution
            iteration_solution = gapfiller.optimize()

            filtered_solution = {rxn:iteration_solution.x_dict[rxn] for rxn in\
               get_fluxes if abs(iteration_solution.x_dict[rxn]) > flux_cutoff}

            add_rxns = [universal.reactions.get_by_id(rxn).copy() for \
                                            rxn in filtered_solution.keys()]
            # combine the solution from this iteration with all others within
            # the current cycle
            cycle_reactions = cycle_reactions | \
                                set([rxn.id for rxn in add_rxns])

            # validate that the proposed solution restores flux through the
            # objective in the original model
            # set the bounds on the original model to represent media
            # and validate the gapfill solution
            for ex_rxn in [rxn for rxn in model_to_gapfill.reactions if \
                            rxn.id.startswith(exchange_prefix)]:
                ex_rxn.lower_bound = 0
            for ex_rxn in phenotype_dict[condition].keys():
                model_to_gapfill.reactions.get_by_id(ex_rxn).lower_bound = \
                    -1.0*phenotype_dict[condition][ex_rxn]
                model_to_gapfill.reactions.get_by_id(ex_rxn).upper_bound = \
                    1.0*phenotype_dict[condition][ex_rxn]
            if not validate(model_to_gapfill,\
                            [universal.reactions.get_by_id(rxn).copy() for \
                            rxn in cycle_reactions],lower_bound):
                raise RuntimeError('Failed to validate gapfilled model, '
                                    'try lowering the flux_cutoff through '
                                    'inclusion_threshold')

            # remove the flux minimization penalty on the gapfilled reactions
            coefficients = (gapfiller.objective.
                        get_linear_coefficients(gapfiller.variables).copy())
            reaction_variables = (((gapfiller.reactions.get_by_id(rxn)
                                    .forward_variable),
                                    (gapfiller.reactions.get_by_id(rxn)
                                    .reverse_variable))
                                     for rxn in cycle_reactions)
            variables = chain(*reaction_variables)
            for variable in variables:
                 coefficients[variable] = 0.0

            gapfiller.objective.set_linear_coefficients(coefficients)
            check = gapfiller.slim_optimize() # optimizing might be necessary
            # to update coefficients.

            # reset the media condition
            for ex_rxn in exchange_reactions:
                ex_rxn.lower_bound = 0

        gapfiller.objective.set_linear_coefficients(original_coefficients)
        solutions.append(list(cycle_reactions))
    return solutions


def _integer_iterative_binary_gapfill(model,phenotype_dict,cycle_order,
                      universal=None, output_ensemble_size=0,
                      lower_bound=0.05, penalties=False,
                      demand_reactions=False,
                      exchange_reactions=False,
                      integer_threshold=1E-6):


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


def _build_ensemble_from_gapfill_solutions(model,solutions,universal=None):

    ensemble = Ensemble(identifier=model.id,name=model.name)
    ensemble.base_model = model.copy()

    # generate member identifiers for each solution
    # Convert the solution to a dictlist so we can retrieve reactions by id
    solution_dict = {}
    i = 0
    for solution in solutions:
        solution_id = model.id + '_gapfilled_' + str(i)
        solution_as_rxn_objs = [universal.reactions.get_by_id(rxn).copy()
                                    for rxn in solution]
        solution_dict[solution_id] = DictList() + solution_as_rxn_objs
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
    ensemble.base_model.add_reactions(all_reactions_as_objects)

    # add metabolite objects to the base model for all new metabolites from
    # the new reactions
    mets = [x.metabolites for x in all_reactions_as_objects]
    all_keys = set().union(*(d.keys() for d in mets))
    ensemble.base_model.add_metabolites(all_keys)

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
                # get the reaction and it's value for the attribute
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

def validate(original_model, reactions, lower_bound):
        with original_model as model:
            model.add_reactions(reactions)
            mets = [x.metabolites for x in reactions]
            all_keys = set().union(*(d.keys() for d in mets))
            model.add_metabolites(all_keys)
            model.slim_optimize()
            return (model.solver.status == OPTIMAL and
                    model.solver.objective.value >= lower_bound)
