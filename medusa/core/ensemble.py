# ensemble_model class

from __future__ import absolute_import

import cobra
import random
import pandas as pd

class Ensemble:
    """
    Ensemble of metabolic models. Contains a base_model of type
    cobra.core.Model which consists of all reactions and genes
    present in any model within the ensemble. Also contains a
    reaction_diff of type dictionary which specifies reaction
    bounds for each model for every reaction that contains a difference
    in any model in the ensemble.

    reaction_diff: {ensemble_member_id:{reaction_id:{'ub':ub,'lb':lb}}}
    base_model: Cobra.core.Model object.
    """

    def __init__(self,model_list=None,base_id=None):
        """
        model_list: list of models to merge into a base model.
        base_id: id that will be used for the new base_model.
        """

        if model_list:
            self.base_model = self._create_base_model(model_list,base_id=base_id+'_base_model')
            self.reaction_diffs = self._create_reaction_diffs(model_list)
        else:
            self.base_model = cobra.core.Model(base_id+'_base_model')
            self.reaction_diffs = {}
        self.id=base_id


    def add_models(self,model_list):
        '''
        Adds an additional model to an existing ensemble. Modifies base_model to
        include any unique reactions and updates the reaction diff.

        When creating an ensemble for more than hundreds of models, creating an
        empty base model, then adding the models individually using add_models
        can prevent the need to load all models simultaneously (e.g. you can
        load and add one at a time, then clear the individual model from memory).
        This methods conserves memory relative to Ensemble_model.init, but may take
        slightly longer to run as the reaction diff is recalculated for each
        addition and more calls to base_model.add_reactions() are made.

        WARNING: currently requires starting from an ensemble object which
        contains at least one base model. If an empty ensemble is used with
        add_models(), the diff is considered to be every reaction. We should
        refactor to enable starting with an empty model, since this makes
        code a little cleaner for the user.
        '''
        new_base_model = self._create_base_model(model_list=model_list)
        new_reaction_diffs = self._create_reaction_diffs(model_list=model_list)

        old_base_rxns = [rxn.id for rxn in self.base_model.reactions]
        new_base_rxns = [rxn.id for rxn in new_base_model.reactions]
        # currently grabs the reaction list from the first model. Requires update
        # if reaction_diff is refactored.
        if len(self.reaction_diffs.keys()) > 0: # check if the ensemble was empty
            old_diff_rxns = [rxn for rxn in self.reaction_diffs[list(self.reaction_diffs.keys())[0]]]
        else:
            old_diff_rxns = []
        new_diff_rxns = [rxn for rxn in new_reaction_diffs[list(new_reaction_diffs.keys())[0]]]

        # find reactions in the new base that aren't in the old_base
        new_base_notin_old_base = list(set(new_base_rxns) - set(old_base_rxns))
        # find reactions in the old base that aren't in the new base
        old_base_notin_new_base = list(set(old_base_rxns) - set(new_base_rxns))
        # find reactions present in both new and old bases.
        in_both_bases = list(set(old_base_rxns) & set(new_base_rxns))

        rxns_to_add_to_old_base = []
        for reaction in new_base_notin_old_base:
            rxns_to_add_to_old_base.append(new_base_model.reactions.get_by_id(reaction))
            for model in self.reaction_diffs.keys():
                # close the reaction in the old reaction diffs
                self.reaction_diffs[model][reaction] = {'lb':0.0,'ub':0.0}
            if reaction not in new_diff_rxns:
                # if the reaction wasn't in the new_diffs, we need to add it with
                # the bounds from the new base since the reaction is absent from the old models
                lb_in_new_base = new_base_model.reactions.get_by_id(reaction).lower_bound
                ub_in_new_base = new_base_model.reactions.get_by_id(reaction).upper_bound
                for model in new_reaction_diffs.keys():
                    new_reaction_diffs[model][reaction] = \
                        {'lb':lb_in_new_base,'ub':ub_in_new_base}

        for reaction in old_base_notin_new_base:
            if reaction not in old_diff_rxns:
                # if the reaction isn't in the old diff, we need to add it with
                # the bounds from old_reaction_diff since it isn't in the new base.
                # (e.g. the reaction is now different in some models)
                lb_in_old_base = self.base_model.reactions.get_by_id(reaction).lower_bound
                ub_in_old_base = self.base_model.reactions.get_by_id(reaction).upper_bound
                for model in self.reaction_diffs.keys():
                    self.reaction_diffs[model][reaction] = \
                        {'lb':lb_in_old_base,'ub':ub_in_old_base}
                # add the reaction to the new diffs, with it closed since it
                # wasn't in the new base.
                for model in new_reaction_diffs.keys():
                    new_reaction_diffs[model][reaction] = {'lb':0.0,'ub':0.0}
            else:
                #if the reaction was in the old_diff, we must add it to new_diff since it's
                # variable in a set of models
                if not reaction in new_diff_rxns: # this should always be true since these reactions aren't in the new base, but just in case.
                    for model in new_reaction_diffs.keys():
                        new_reaction_diffs[model][reaction] = {'lb':0.0,'ub':0.0}

        for reaction in in_both_bases:
            if reaction in old_diff_rxns:
                if reaction not in new_diff_rxns:
                    # add to new diff with bounds from new_base
                    lb_in_new_base = new_base_model.reactions.get_by_id(reaction).lower_bound
                    ub_in_new_base = new_base_model.reactions.get_by_id(reaction).upper_bound
                    for model in new_reaction_diffs.keys():
                        new_reaction_diffs[model][reaction] = \
                            {'lb':lb_in_new_base,'ub':ub_in_new_base}
            if reaction in new_diff_rxns:
                if reaction not in old_diff_rxns:
                    # add to old diff with bound from old_base
                    lb_in_old_base = self.base_model.reactions.get_by_id(reaction).lower_bound
                    ub_in_old_base = self.base_model.reactions.get_by_id(reaction).upper_bound
                    for model in self.reaction_diffs.keys():
                        self.reaction_diffs[model][reaction] = \
                            {'lb':lb_in_old_base,'ub':ub_in_old_base}
            else:
                # if the reaction is in neither diff, check whether the bounds
                # are the same in both old and new base, then update diffs
                lb_in_new_base = new_base_model.reactions.get_by_id(reaction).lower_bound
                ub_in_new_base = new_base_model.reactions.get_by_id(reaction).upper_bound
                lb_in_old_base = self.base_model.reactions.get_by_id(reaction).lower_bound
                ub_in_old_base = self.base_model.reactions.get_by_id(reaction).upper_bound
                if lb_in_new_base != lb_in_old_base or ub_in_new_base != ub_in_old_base:
                    for model in self.reaction_diffs.keys():
                        self.reaction_diffs[model][reaction] = \
                            {'lb':lb_in_old_base,'ub':ub_in_old_base}
                    for model in new_reaction_diffs.keys():
                        new_reaction_diffs[model][reaction] = \
                            {'lb':lb_in_new_base,'ub':ub_in_new_base}
        # update the base_model and reaction_diff
        self.base_model.add_reactions(rxns_to_add_to_old_base)
        self.base_model.repair()
        for model in new_reaction_diffs.keys():
            self.reaction_diffs[model] = new_reaction_diffs[model]

    def _apply_diffs(self,model_name):
        '''
        apply for the reaction diffs for the given model_name. Used as a helper
        function to more easily utilize utility functions from cobra (e.g.
        single reaction/gene deletions).

        '''
        diffs = self.reaction_diffs[model_name]
        for reaction in diffs.keys():
            rxn = self.base_model.reactions.get_by_id(reaction)
            rxn.lower_bound = diffs[reaction]['lb']
            rxn.upper_bound = diffs[reaction]['ub']

    def _create_base_model(self,model_list,base_id='ensemble base'):
        '''
        Creates a base model from a list of universal models. Base models
        contain all reactions present in a list of models. A reaction_diff
        object is also returned, indicating the bounds for any reaction that
        is different across any model in the ensemble.
        '''
        base_model = cobra.Model(base_id)
        for model in model_list:
            rxns_to_add = []
            for reaction in model.reactions:
                if reaction.id not in [rxn.id for rxn in base_model.reactions]:
                    rxns_to_add.append(reaction.copy())
            base_model.add_reactions(rxns_to_add)
            base_model.repair()

        return base_model

    def _create_reaction_diffs(self,model_list):
        '''


        '''
        # Get list of every reaction present in all models
        reaction_set = set()
        for model in model_list:
            model_reactions = set([reaction.id for reaction in model.reactions])
            if len(reaction_set) > 0:
                reaction_set = reaction_set & model_reactions
            else:
                reaction_set = model_reactions
        reaction_list = list(reaction_set)

        # Get list of every reaction that isn't common to all models
        variable_reactions = set()
        for model in model_list:
            reactions_different = [rxn.id for rxn in model.reactions if rxn.id not in reaction_list]
            if len(variable_reactions) > 0:
                variable_reactions = variable_reactions | set(reactions_different)
            else:
                variable_reactions = set(reactions_different)
        variable_reaction_list = list(variable_reactions)


        # For every reaction that isn't common to all models, find the
        # upper and lower bound. Reactions that aren't present in a model
        # have upper and lower bounds of 0 to prevent any flux.
        reaction_diffs = {}
        for model in model_list:
            reaction_diffs[model.id] = {}
            for reaction in variable_reaction_list:
                if reaction in [rxn.id for rxn in model.reactions]:
                    rxn_in_model = model.reactions.get_by_id(reaction)
                    lb = rxn_in_model.lower_bound
                    ub = rxn_in_model.upper_bound
                    reaction_diffs[model.id][reaction] = {'lb':lb,'ub':ub}
                else:
                    # if the reaction wasn't in the model, it should be closed
                    # in both directions for this member of the ensemble
                    reaction_diffs[model.id][reaction] = {'lb':0.0,'ub':0.0}
        return reaction_diffs

    def ensemble_single_reaction_deletion(self,optimal_only=True,num_models=[]):
        '''
        IN PROGRESS, not functional

        Performs single reaction deletions for num_models in the ensemble. Reports
        absolute objective value, and can be modified to return only optimal solutions.

        Currently does not return the solver status for all members of an ensemble,
        so infeasible and 0 flux values can't be discriminated unless optimal_only=True
        is passed (in which case feasible solutions with flux through the objective
        of 0 are returned, but infeasible solutions are not)
        '''
        if not num_models:
            # if not specified, use all models
            num_models = len(self.reaction_diffs.keys())

        # for each model, perform all the deletions, then advance to the next model.
        for model in random.sample(list(self.reaction_diffs.keys()),num_models):
            # set the correct bounds for the model
            self._apply_diffs(model)
            # perform single reaction deletion for all reactions in the model
            # TODO: make exception for biomass and the demand reaction.
            print('performing deletions for ' + model)
            model_deletion_results = cobra.single_reaction_deletion(self.base_model,self.base_model.reactions)

        # if optimal_only, filter the output dataframe for each model to exclude infeasibles,
        # then append to a master dataframe


    def leak_test(self,metabolites_to_test=[],\
                 exchange_prefix='EX_',verbose=False,num_models=[],**kwargs):
        '''
        Checks for leaky metabolites in every member of the ensemble by opening
        and optimizing a demand reaction while all exchange reactions are closed.

        By default, checks for leaks for every metabolite for all models.
        '''

        if not num_models:
            # if the number of models wasn't specified, test all
            num_models = len(self.reaction_diffs.keys())

        if not metabolites_to_test:
            metabolites_to_test = [met for met in self.base_model.metabolites]

        old_objective = self.base_model.objective
        dm_rxns = []
        for met in metabolites_to_test:
            rxn = cobra.Reaction(id='leak_DM_' + met.id)
            rxn.lower_bound = 0.0
            rxn.upper_bound = 0.0
            rxn.add_metabolites({met:-1})
            dm_rxns.append(rxn)
        self.base_model.add_reactions(dm_rxns)
        self.base_model.repair()

        leaks = {}

        for rxn in dm_rxns:
            rxn.upper_bound = 1000.0
            self.base_model.objective = self.base_model.reactions.get_by_id(rxn.id)

            if verbose:
                print('checking leak for ' + rxn.id)
            solutions = self.optimize_ensemble(return_flux=[rxn.id],num_models=num_models,**kwargs)
            leaks[rxn.id.split('_DM_')[1]] = {}
            for model in solutions.keys():
                leaks[rxn.id.split('_DM_')[1]][model] = solutions[model][rxn.id] > 0.0001
            #rxn.objective_coefficient = 0.0
            rxn.upper_bound = 0.0

        # remove the demand reactions and restore the original objective
        self.base_model.remove_reactions(dm_rxns,remove_orphans=True)
        self.base_model.repair()
        self.base_model.objective = old_objective

        return leaks
