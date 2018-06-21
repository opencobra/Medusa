

def ensemble_single_reaction_deletion(ensemble,optimal_only=True,num_models=[]):
    '''
    IN PROGRESS, not functional. Need to refactor for features/states.

    Performs single reaction deletions for num_models in the ensemble. Reports
    absolute objective value, and can be modified to return only optimal solutions.

    Currently does not return the solver status for all members of an ensemble,
    so infeasible and 0 flux values can't be discriminated unless optimal_only=True
    is passed (in which case feasible solutions with flux through the objective
    of 0 are returned, but infeasible solutions are not)
    '''
    if not num_models:
        # if not specified, use all models
        num_models = len(ensemble.reaction_diffs.keys())

    # for each model, perform all the deletions, then advance to the next model.
    for model in random.sample(list(ensemble.reaction_diffs.keys()),num_models):
        # set the correct bounds for the model
        ensemble._apply_diffs(model)
        # perform single reaction deletion for all reactions in the model
        # TODO: make exception for biomass and the demand reaction.
        print('performing deletions for ' + model)
        model_deletion_results = cobra.single_reaction_deletion(self.base_model,self.base_model.reactions)

    # if optimal_only, filter the output dataframe for each model to exclude infeasibles,
    # then append to a master dataframe
