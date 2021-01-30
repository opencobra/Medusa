# Import packages
import mackinac
import cobra
import pandas as pd
import json
import os
import numpy as np
import medusa
from pickle import load
from termcolor import colored
import matplotlib.pyplot as plt
import pickle
from os.path import join 
import os.path
import itertools
import random
from itertools import chain
import os, json
from optlang.interface import OPTIMAL

from cobra.flux_analysis.gapfilling import GapFiller
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.core import DictList

from cobra.util.solver import linear_reaction_coefficients

from medusa.core.ensemble import Ensemble
from medusa.core.feature import Feature
from medusa.core.member import Member
# from probanno import probanno
REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}


## Genral Functions
def gapfill_to_ensemble(model, iterations=1, universal=None, lower_bound=0.05,
                 penalties= None, exchange_reactions=False,
                 demand_reactions=False, integer_threshold=1e-6):
    gapfiller = GapFiller(model, universal=universal,
                          lower_bound=lower_bound, penalties= 1-reaction_probability,
                          demand_reactions=demand_reactions,
                          exchange_reactions=exchange_reactions,
                          integer_threshold=integer_threshold)
    # update the linear coefficients for the added reaction to the penalitties of the reaction
    new_coefficients =  {}
    for reaction in [rxn.id for rxn in gapfiller.reactions]:
        if reaction in reaction_probability.keys():
            new_coefficients[reaction] = 1-reaction_probability[reaction]
    for react_id in new_coefficients.keys():
        reaction = gapfiller.reactions.get_by_id(react_id)
        for_var = reaction.forward_variable
        rev_var = reaction.reverse_variable
        if coefficients[for_var]>0.0:
            coefficients[for_var] = new_coefficients[react_id]
        if coefficients[rev_var]>0.0:
            coefficients[rev_var] = new_coefficients[react_id]
    gapfiller.objective.set_linear_coefficients(coefficients)
    solutions = gapfiller.fill(iterations=iterations)
    print("finished gap-filling. Constructing ensemble...")
    ensemble = _build_ensemble_from_gapfill_solutions(model,solutions,
                                                    universal=universal)
    return ensemble


def _continuous_iterative_binary_gapfill(model = None, 
                      universal=None, output_ensemble_size=1,
                      lower_bound=0.05, reaction_probability=None,
                      demand_reactions=False,
                      exchange_reactions=False,
                      flux_cutoff=1E-8,
                      exchange_prefix='EX_'):
    """ 
    This function takes the gapfilled removed_model, penality, universla model,
    demand_reaction, exchange reaction, flux_cutff and exchanges_prefix. It returns 
    the undated a solution with as set of coefficients of the added reactions
    """
    
    if exchange_reactions:
        raise NotImplementedError("Inclusion of new exchange reactions is not"
                            "supported for continuous gapfill")
    if demand_reactions:
        raise NotImplementedError("Inclusion of demand reactions is not"
                            "supported for continuous gapfill")
    solutions = []
    gapfiller = universal.copy()
    model.reactions.EX_cpd00007_e.lower_bound = 0
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

    # set the linear coefficients for reactions in the original model to 0 and 
    
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
#     # update the linear coefficients for the added reaction to the penalitties of the reaction
#     new_coefficients =  {}
#     for reaction in [rxn.id for rxn in gapfiller.reactions]:
#         if reaction in reaction_probability.keys():
#             new_coefficients[reaction] = 1-reaction_probability[reaction]
#     for react_id in new_coefficients.keys():
#         reaction = gapfiller.reactions.get_by_id(react_id)
#         for_var = reaction.forward_variable
#         rev_var = reaction.reverse_variable
#         if coefficients[for_var]>0.0:
#             coefficients[for_var] = new_coefficients[react_id]
#         if coefficients[rev_var]>0.0:
#             coefficients[rev_var] = new_coefficients[react_id]
#     gapfiller.objective.set_linear_coefficients(coefficients)

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
    # Generate the coefficients list and Randomly varying the coefficients 
    iou = []
    for y in coefficients.values():
        iou.append(y)
    solutions = []
    t = 0.05
    for n in range(100):
        new_coeff = iou + np.random.normal(0, t, len(iou))
        coefficients = dict(zip(coefficients.keys(), abs(new_coeff)))
                           
        gapfiller.objective.set_linear_coefficients(coefficients)
        for reaction in original_objective.keys():
            print("Constraining lower bound for " + reaction)
            gapfiller.reactions.get_by_id(reaction).lower_bound = lower_bound
        cycle_reactions = set()
        iteration_solution = gapfiller.optimize()

        filtered_solution = {rxn:iteration_solution.fluxes[rxn] for rxn in\
            get_fluxes if abs(iteration_solution.fluxes[rxn]) > flux_cutoff}

        add_rxns = [universal.reactions.get_by_id(rxn).copy() for \
                                            rxn in filtered_solution.keys()]

        cycle_reactions = cycle_reactions | \
                                set([rxn.id for rxn in add_rxns])
        if not validate(model_to_gapfill,\
                        [universal.reactions.get_by_id(rxn).copy() for \
                        rxn in cycle_reactions],lower_bound):
            raise RuntimeError('Failed to validate gapfilled model, '
                                    'try lowering the flux_cutoff through '
                                    'inclusion_threshold')
        solutions.append(cycle_reactions)

    return solutions

def _build_ensemble_from_gapfill_solutions(model,solutions,universal=None):
    
#     model.reactions.EX_cpd00007_e.lower_bound = 0 # out of the function on medusa
    ensemble = Ensemble(identifier=model.id,name=model.name)
    ensemble.base_model = model.copy()
    """
    This function takes the model, solution and universla model 
    and it will generate the ensembles of gapfilled metabolic model
    model: the gapfilled_removed model
    universal: the universal model for that class of bacteria
    solution: The set of reactions that generated to gapfilled models
    """

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
    # update members for the ensemble and save the ensembles as pkl file
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
