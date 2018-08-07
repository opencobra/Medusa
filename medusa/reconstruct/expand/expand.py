
import itertools

from cobra.flux_analysis.gapfilling import GapFiller
from cobra.core import DictList

from medusa.core.ensemble import Ensemble
from medusa.core.feature import Feature
from medusa.core.member import Member
# functions for expanding existing models to generate an ensemble

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

def gapfill_to_ensemble(model, iterations=1, universal=None, lower_bound=0.05,
                 penalties=None, exchange_reactions=False,
                 demand_reactions=True, integer_threshold=1e-6):
    """
    Performs gapfilling on model, pulling reactions from universal.
    Any existing constraints on base_model are maintained during gapfilling, so
    these should be set before calling gapfill_to_ensemble (e.g. secretion of
    metabolites, choice of objective function etc.).

    Currently, only iterative solutions are supported with accumulating penalties
    (i.e. after each iteration, the penalty for each reaction doubles).

    model : cobra.Model
        The model to perform gap filling on.
    universal : cobra.Model
        A universal model with reactions that can be used to complete the
        model.
    lower_bound : float
        The minimally accepted flux for the objective in the filled model.
    penalties : dict, None
        A dictionary with keys being 'universal' (all reactions included in
        the universal model), 'exchange' and 'demand' (all additionally
        added exchange and demand reactions) for the three reaction types.
        Can also have reaction identifiers for reaction specific costs.
        Defaults are 1, 100 and 1 respectively.
    integer_threshold : float
        The threshold at which a value is considered non-zero (aka
        integrality threshold). If gapfilled models fail to validate,
        you may want to lower this value. However, picking a threshold that is
        too low may also result in reactions being added that are not essential
        to meet the imposed constraints.
    exchange_reactions : bool
        Consider adding exchange (uptake) reactions for all metabolites
        in the model.
    demand_reactions : bool
        Consider adding demand reactions for all metabolites.
    """
    gapfiller = GapFiller(model, universal=universal,
                          lower_bound=lower_bound, penalties=penalties,
                          demand_reactions=demand_reactions,
                          exchange_reactions=exchange_reactions,
                          integer_threshold=integer_threshold)
    solutions = gapfiller.fill(iterations=iterations)

    print("finished gap-filling. Constructing ensemble...")

    ensemble = Ensemble(identifier=model.identifier,name=model.name)
    ensemble.base_model = model.copy()

    # generate member identifiers for each solution
    # Convert the solution to a dictlist so we can retrieve reactions by id
    solution_dict = {}
    i = 0
    for solution in solutions:
        solution_id = model.identifier + '_gapfilled_' + str(i)
        solution_dict[solution_id] = DictList() + solution
        i += 1

    # remove non-unique solutions
    solutions_as_ids = {}
    for member_id in solution_dict.keys():
        solutions_as_ids[member_id] = [rxn.id for rxn in solution_dict[member_id]]

    #TODO: scan through other members and remove them if they are identical
    for member_id in solutions_as_ids:
        for other_member in
        unique_rxns = 0
        for rxn in sultion_dict[member_id]:

    solutions = list(solutions for solutions,_ in itertools.groupby(solutions))


    # parse the solutions to find reactions that are always added and those
    # that vary.
    all_reactions = set()
    in_all = set()
    for solution in solutions:
        all_reactions = all_reactions & solution
        in_all = in_all + all_reactions

    reactions_needing_features = all_reactions - in_all

    # add reaction objects to the base model for all reactions
    ensemble.base_model.add_reactions([rxn.copy for rxn in all_reactions])

    print('building features...')
    # generate features for the variable reactions and add them to the ensemble.
    # assume that all reactions have the same attribute values; if different
    # attribute values are desired for reactions with the same ID, these need
    # to be added to the universal reaction bag prior to gapfilling
    features = DictList()
    for reaction in reactions_needing_features:
        for attribute in REACTION_ATTRIBUTES:
            identifier = reaction.id + attribute
            name = reaction.name
            feature = Feature(identifier,name)
            feature.ensemble = ensemble
            feature.base_component = ensemble.base_model.reactions.get_by_id(reaction.id)
            feature.component_attribute = attribute
            # get the states for this feature in {member.id:value}
            states = {}
            for member_id in solution_dict.keys():
                # get the reaction and it's value for the attribute
                if reaction in solution_dict[solution]:
                    rxn_obj = solution_dict[member_id].get_by_id(reaction.id)
                    states[member_id] = getattr(rxn_obj,attribute)
            feature.states = states
            features += [feature]

    print('updating members...')
    # update members for the ensemble
    ensemble.members = DictList()
    for member_id in solution_dict.keys():
        model_states = dict()
        for feature in self.features:
            model_states[feature] = feature.get_model_state(member_id)
        member = Member(ensemble=ensemble,\
                        identifier=member_id,\
                        name=ensemble.name,\
                        states=model_states)
        ensemble.members += [member]

    return ensemble
