from pandas import DataFrame
from random import sample
from cobra.flux_analysis.variability import (
    flux_variability_analysis, find_blocked_reactions,
    find_essential_genes, find_essential_reactions)

def ensemble_fva(ensemble, reaction_list, num_models=[],specific_models=None,
                 fraction_of_optimum=1.0, loopless=False,**solver_args):
    '''
    Performs FVA on num_models. If num_models is not passed, performs FVA
    on every model in the ensemble. If the model is a community model,
    num_models must be passed.

    Performs flux variability analysis (FVA) on models within an ensemble.

    Parameters
    ----------
    ensemble: medusa.core.Ensemble
        The ensemble on which FVA is to be performed.
    reaction_list: str or list of str, optional
        List of reaction ids (cobra.core.reaction.id), or a single reaction id,
        for which to return flux ranges. If None, all reaction fluxes are
        returned (default).
    num_models: int, optional
        Number of models for which FVA will be performed. The number of models
        indicated will be randomly sampled and FVA will be performed on the
        sampled models. If None, all models will be selected (default), or the
        models specified by specific_models will be selected. Cannot be passed
        concurrently with specific_models.
    specific_models: list of str, optional
        List of ensemble_member.id corresponding to the models for which FVA
        will be performed. If None, all models will be selected (default), or
        num_models will be randomly sampled and selected. Cannot be passed
        concurrently with num_models.
    fraction_of_optimum: float, optional
        fraction of the optimum objective value, set as a constraint such that
        the objective never falls below the provided fraction when assessing
        variability of each reaction.
    loopless: boolean, optional
        Whether or not to perform loopless FVA. This is much slower. See
        cobrapy.flux_analysis.variability for details.

    Returns
    -------
    pandas.DataFrame
        A dataframe in which each row (index) represents a model within the
        ensemble and the lower or upper value of flux ranges, and each column
        represents a reaction and its lower or upper value of its flux range.
        Based on this formatting, each model is present in two rows, one of
        which contains the lower flux value and the other of which contains
        the upper flux value.

    '''
    if not num_models:
        # if not specified, use all models
        num_models = len(ensemble.members)

    if isinstance(reaction_list,str):
        reaction_list = [reaction_list]

    # initialize dataframe to store results. max and min for each model will
    # each take a single row, and the columns will be reactions.
    all_fva_results = DataFrame()

    if specific_models:
        model_list = specific_models
    else:
        model_list = sample(ensemble.members,num_models)
        model_list = [model.id for model in model_list]
    with ensemble.base_model:
        for model in model_list:
            ensemble.set_state(model)
            fva_result = flux_variability_analysis(
                    ensemble.base_model,reaction_list=reaction_list,
                    fraction_of_optimum=fraction_of_optimum,
                    loopless=loopless,**solver_args)
            fva_result.columns = ['maximum_'+model,
                                    'minimum_'+model]
            fva_result = fva_result.T
            # add the model source as a column
            fva_result['model_source'] = [model,model]
            all_fva_results = all_fva_results.append(fva_result)
    return all_fva_results
