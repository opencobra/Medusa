
import itertools
import random
from itertools import chain

from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

from cobra.flux_analysis.gapfilling import GapFiller
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.flux_analysis import find_blocked_reactions
from cobra import sampling
from cobra.core import DictList

from cobra.util.solver import linear_reaction_coefficients
from cobra.util import solver as sutil

from medusa.core.ensemble import Ensemble
from medusa.core.feature import Feature
from medusa.core.member import Member
# functions for expanding existing models to generate an ensemble

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

GAPFILL_METHODS = ["integer","continuous"]

def gapfill_to_ensemble(model, iterations=1, universal=None, lower_bound=0.05,
                 gapfill_type='continuous', penalties=None, 
                 exchange_reactions=False, demand_reactions=False,
                 integer_threshold=1e-6, sampler_method='optgp'):
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
    sampler_method : String, 'optgp'
        For continuous pFBA-based gapfilling, specifies the sampling method
        to be used. Takes any valid sampling method implemented in cobrapy,
        e.g., 'optgp' or 'achr'.

    Returns
    -------
    ensemble : medusa.core.Ensemble
        The ensemble object created from the gapfill solutions.
    """

    if gapfill_type not in GAPFILL_METHODS:
        raise ValueError("gapfill_type \'{0}\' is not supported."
                         "Must be one of {1}.".format(
                         gapfill_type, GAPFILL_METHODS))

    if gapfill_type == "integer":
        gapfiller = GapFiller(model, universal=universal,
                              lower_bound=lower_bound, penalties=penalties,
                              demand_reactions=demand_reactions,
                              exchange_reactions=exchange_reactions,
                              integer_threshold=integer_threshold)
        solutions = gapfiller.fill(iterations=iterations)
    
    elif gapfill_type == "continuous":
        print("Constructing pFBA gapfiller")
        gapfiller = PfbaGapFiller(model, universal=universal,
                              lower_bound=lower_bound, penalties=penalties,
                              demand_reactions=demand_reactions,
                              exchange_reactions=exchange_reactions)
        solutions = gapfiller.fill(iterations=iterations,
                                   sampler_method = sampler_method)

    print("finished gapfilling. Constructing ensemble...")
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
    if gapfill_type not in GAPFILL_METHODS:
        raise ValueError("gapfill_type \'{0}\' is not supported."
                         "Must be one of {1}.".format(
                         gapfill_type, GAPFILL_METHODS))

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

            filtered_solution = {rxn:iteration_solution.fluxes[rxn] for rxn in\
               get_fluxes if abs(iteration_solution.fluxes[rxn]) > flux_cutoff}

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

class PfbaGapFiller(object):
    def __init__(self, model, universal=None, lower_bound=0.05,
                 penalties=None, exchange_reactions=False,
                 demand_reactions=False, threshold=1e-8):
        print("\tCopying original model and universal model")
        self.original_model = model
        self.lower_bound = lower_bound
        self.model = model.copy()
        #tolerances = self.model.solver.configuration.tolerances
        self.universal = universal.copy() if universal else Model('universal')
        self.penalties = dict(universal=1, exchange=100, demand=1)
        if penalties is not None:
            self.penalties.update(penalties)
        self.extend_model(exchange_reactions, demand_reactions)
        print("\tRemoving blocked reactions from the combined model "
            "and universal")
        self.prune_blocked_reactions()
        #fix_objective_as_constraint(self.model, bound=lower_bound)
        print("\tAdding pFBA objective")
        add_pfba(self.model,lower_bound = lower_bound)
        
        # set the linear coefficients for reactions in the original model to 0
        print("\tRemoving flux penalty from reactions in the original model")
        coefficients = (self.model.objective
                        .get_linear_coefficients(self.model.variables))
        # get the original reactions from the model, excluding those that
        # were always blocked even with the addition of the universal.
        original_reaction_ids = [r.id for r in self.original_model.reactions
                                if r.id in 
                                [rxn.id for rxn in self.model.reactions]]
        reaction_variables = (((self.model.reactions.get_by_id(reaction)
                                .forward_variable),
                               (self.model.reactions.get_by_id(reaction)
                                .reverse_variable))
                                for reaction in original_reaction_ids)
        variables = chain(*reaction_variables)
        for variable in variables:
            coefficients[variable] = 0.0
        self.model.objective.set_linear_coefficients(coefficients)



    def extend_model(self, exchange_reactions=False, demand_reactions=True):
        """Extend gapfilling model.
        Add reactions from universal model and optionally exchange and
        demand reactions for all metabolites in the model to perform
        gapfilling on.
        Parameters
        ----------
        exchange_reactions : bool
            Consider adding exchange (uptake) reactions for all metabolites
            in the model.
        demand_reactions : bool
            Consider adding demand reactions for all metabolites.
        """
        for rxn in self.universal.reactions:
            rxn.gapfilling_type = 'universal'
        new_metabolites = self.universal.metabolites.query(
            lambda metabolite: metabolite not in self.model.metabolites
                                                           )
        self.model.add_metabolites(new_metabolites)
        existing_exchanges = []
        for rxn in self.universal.boundary:
            existing_exchanges = existing_exchanges + \
                [met.id for met in list(rxn.metabolites)]

        for met in self.model.metabolites:
            if exchange_reactions:
                # check for exchange reaction in model already
                if met.id not in existing_exchanges:
                    rxn = self.universal.add_boundary(
                        met, type='exchange_gapfill', lb=-1000, ub=0,
                        reaction_id='EX_{}'.format(met.id))
                    rxn.gapfilling_type = 'exchange'
            if demand_reactions:
                rxn = self.universal.add_boundary(
                    met, type='demand_gapfill', lb=0, ub=1000,
                    reaction_id='DM_{}'.format(met.id))
                rxn.gapfilling_type = 'demand'

        new_reactions = self.universal.reactions.query(
            lambda reaction: reaction not in self.model.reactions
        )
        self.model.add_reactions(new_reactions)

    def prune_blocked_reactions(self):
        blocked = find_blocked_reactions(self.model)
        self.model.remove_reactions(blocked)


    def fill(self, iterations=1, threshold=None, sampler_method='optgp'):
            print("\nBeginning gapfill")
            if not threshold:
                threshold = self.model.tolerance

            # only get the fluxes for reactions that weren't in the original
            # model
            get_fluxes = [r.id for r in self.model.reactions if 
                        r.id not in 
                        [rxn.id for rxn in self.original_model.reactions]]

            # if performing more than one iteration, generate a sampler
            if iterations > 1:
                print("\tGenerating samples")
                samples = sampling.sample(self.model, 
                                 iterations, 
                                 method = sampler_method)
                filtered_samples = samples > threshold

                # reduce the samples to only reactions in get_fluxes
                filtered_samples = filtered_samples[get_fluxes]

                # collapse redundant rows
                filtered_samples = filtered_samples.drop_duplicates()
                # convert to list of lists with active reactions
                filtered_solution = [
                                    filtered_samples.columns[
                                    filtered_samples.loc[index]].tolist()
                                    for index in filtered_samples.index]
                for solution in filtered_solution:
                    if not validate(self.original_model,[
                        self.model.reactions.get_by_id(rxn).copy() for
                        rxn in solution]):
                        raise RuntimeError('Failed to validate gapfilled model, '
                                    'try lowering the flux_cutoff through '
                                    'inclusion_threshold')

            else:
                iter_sol = self.model.optimize()
                filtered_solution = [rxn for rxn in get_fluxes if 
                                    abs(iter_sol.fluxes[rxn]) > threshold]
                if not validate(self.original_model,
                                    [self.model.reactions.get_by_id(rxn).copy()
                                    for rxn in filtered_solution.keys()]):
                    raise RuntimeError('Failed to validate gapfilled model, '
                                    'try lowering the flux_cutoff through '
                                    'inclusion_threshold')
                filtered_solution = [filtered_solution] # list of lists
                # to be consistent with iterations > 1 case.
            
            return filtered_solution



def add_pfba(model, objective=None, fraction_of_optimum=1.0, lower_bound=None):
    """Add pFBA objective modified for medusa
    Add objective to minimize the summed flux of all reactions to the
    current objective.
    See Also
    -------
    pfba
    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum. Ignored if lower_bound is passed.
    lower_bound : float
        Lower bound on the former objective that must be maintained. The
        original objective is constrained in the same manner as when
        frasction_of_optimum is used. If used, fraction_of_optimum is ignored.
    """
    if objective is not None:
        model.objective = objective
    if model.solver.objective.name == '_pfba_objective':
        raise ValueError('The model already has a pFBA objective.')
    if lower_bound:
        sutil.fix_objective_as_constraint(model, bound=lower_bound)
    else:
        sutil.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
    reaction_variables = ((rxn.forward_variable, rxn.reverse_variable)
                          for rxn in model.reactions)
    variables = chain(*reaction_variables)
    model.objective = model.problem.Objective(
        Zero, direction='min', sloppy=True, name="_pfba_objective")
    model.objective.set_linear_coefficients({v: 1.0 for v in variables})

def validate(original_model, reactions, lower_bound):
        with original_model as model:
            model.add_reactions(reactions)
            mets = [x.metabolites for x in reactions]
            all_keys = set().union(*(d.keys() for d in mets))
            model.add_metabolites(all_keys)
            model.slim_optimize()
            return (model.solver.status == OPTIMAL and
                    model.solver.objective.value >= lower_bound)
