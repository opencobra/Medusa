

def leak_test(ensemble,metabolites_to_test=[],\
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

    old_objective = ensemble.base_model.objective
    dm_rxns = []
    for met in metabolites_to_test:
        rxn = cobra.Reaction(id='leak_DM_' + met.id)
        rxn.lower_bound = 0.0
        rxn.upper_bound = 0.0
        rxn.add_metabolites({met:-1})
        dm_rxns.append(rxn)
    ensemble.base_model.add_reactions(dm_rxns)
    ensemble.base_model.repair()

    leaks = {}

    for rxn in dm_rxns:
        rxn.upper_bound = 1000.0
        ensemble.base_model.objective = ensemble.base_model.reactions.get_by_id(rxn.id)

        if verbose:
            print('checking leak for ' + rxn.id)
        solutions = ensemble.optimize_ensemble(return_flux=[rxn.id],num_models=num_models,**kwargs)
        leaks[rxn.id.split('_DM_')[1]] = {}
        for model in solutions.keys():
            leaks[rxn.id.split('_DM_')[1]][model] = solutions[model][rxn.id] > 0.0001
        #rxn.objective_coefficient = 0.0
        rxn.upper_bound = 0.0

    # remove the demand reactions and restore the original objective
    ensemble.base_model.remove_reactions(dm_rxns,remove_orphans=True)
    ensemble.base_model.repair()
    ensemble.base_model.objective = old_objective

    return leaks
