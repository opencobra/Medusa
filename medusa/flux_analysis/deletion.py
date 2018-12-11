from cobra.flux_analysis import single_reaction_deletion, single_gene_deletion
from random import sample

def ensemble_single_reaction_deletion(ensemble, num_models=None,
                                        specific_models=[]):
    '''
    Performs single reaction deletions on models within an ensemble and
    returns the objective value after optimization with each reaction removed.

    Parameters
    ----------
    ensemble: medusa.core.Ensemble
        The ensemble with which to perform reaction deletions
    num_models: int, optional
        Number of models for which reaction deletions will be performed. The
        number of models indicated will be randomly sampled and reaction
        deletions will be performed on the sampled models. If None, all models
        will be selected (default), or the models specified by specific_models
        will be selected. Cannot be passed concurrently with specific_models.
    specific_models: list of str, optional
        List of member.id corresponding to the models for which reaction
        deletions will be performed. If None, all models will be selected
        (default), or num_models will be randomly sampled and selected.
        Cannot be passed concurrently with num_models.

    Returns
    -------
    pandas.DataFrame
        A dataframe in which each row (index) represents a model within the
        ensemble, and each column represents a reaction for which values of
        objective when the reaction is deleted are returned.
    '''
    if not num_models:
        num_models = len(ensemble.members)

    if specific_models:
        model_list = specific_models
    else:
        model_list = sample(ensemble.members,num_models)

    deletion_results = {}
    with ensemble.base_model:
        for model in model_list:
            print('performing deletions for ' + model.id)
            ensemble.set_state(model)
            deletion_result = single_reaction_deletion(ensemble.base_model)
            deletion_results[model.id] = deletion_result

    return deletion_results

def ensemble_single_gene_deletion(ensemble, num_models=None,
                                        specific_models=[],
                                        specific_genes=[]):
    '''
    Performs single reaction deletions on models within an ensemble and
    returns the objective value after optimization with each reaction removed.

    Parameters
    ----------
    ensemble: medusa.core.Ensemble
        The ensemble with which to perform reaction deletions
    num_models: int, optional
        Number of models for which reaction deletions will be performed. The
        number of models indicated will be randomly sampled and reaction
        deletions will be performed on the sampled models. If None, all models
        will be selected (default), or the models specified by specific_models
        will be selected. Cannot be passed concurrently with specific_models.
    specific_models: list of str, optional
        List of member.id corresponding to the models for which reaction
        deletions will be performed. If None, all models will be selected
        (default), or num_models will be randomly sampled and selected.
        Cannot be passed concurrently with num_models.
    specific_genes: list of str, optionsl
        List of gene.id corresponding to the genes for which deletions
        should be performed. If none, all genes will be selected (default).
        We recommend identifying genes that are essential in all ensemble
        members first, then excluding those genes from specific_genes.
        This will generally speed up computation.

    Returns
    -------
    pandas.DataFrame
        A dataframe in which each row (index) represents a model within the
        ensemble, and each column represents a reaction for which values of
        objective when the reaction is deleted are returned.
    '''
    if not num_models:
        num_models = len(ensemble.members)

    if specific_models:
        model_list = specific_models
    else:
        model_list = sample(ensemble.members,num_models)

    deletion_results = {}
    with ensemble.base_model:
        for model in model_list:
            print('performing deletions for ' + model.id)
            ensemble.set_state(model)
            deletion_result = single_gene_deletion(ensemble.base_model,specific_genes)
            deletion_results[model.id] = deletion_result

    return deletion_results
