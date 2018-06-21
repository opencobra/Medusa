
from __future__ import absolute_import

from pandas import DataFrame
from random import sample

def optimize_ensemble(ensemble,return_flux=None,num_models=None,specific_models=None,**kwargs):
    '''
    Performs flux balance analysis (FBA) on models within an ensemble.

    Parameters
    ----------
    ensemble: medusa.core.Ensemble
        The ensemble on which FBA is to be performed.
    return_flux: str or list of str, optional
        List of reaction ids (cobra.core.reaction.id), or a single reaction id,
        for which to return flux values. If None, all reaction fluxes are
        returned (default).
    num_models: int, optional
        Number of models for which FBA will be performed. The number of models
        indicated will be randomly sampled and FBA will be performed on the
        sampled models. If None, all models will be selected (default), or the
        models specified by specific_models will be selected. Cannot be passed
        concurrently with specific_models.
    specific_models: list of str, optional
        List of ensemble_member.id corresponding to the models for which FBA
        will be performed. If None, all models will be selected (default), or
        num_models will be randomly sampled and selected. Cannot be passed
        concurrently with num_models.

    Returns
    -------
    pandas.DataFrame
        A dataframe in which each row (index) represents a model within the
        ensemble, and each column represents a reaction for which flux values
        are returned.
    '''

    if not num_models:
        num_models = len(ensemble.members)

    if isinstance(return_flux,str):
        return_flux = [return_flux]

    flux_dict = {}
    return_vals = {}
    if specific_models:
        model_list = specific_models
    else:
        model_list = sample(ensemble.members,num_models)
    with ensemble.base_model:
        for model in model_list:
            ensemble.set_state(model)
            ensemble.base_model.optimize(**kwargs)
            flux_dict[model] = {rxn.id:rxn.flux for rxn in ensemble.base_model.reactions}
    if return_flux:
        for model in flux_dict.keys():
            return_vals[model] = {rxn_id:flux_dict[model][rxn_id] for rxn_id in return_flux}
    else:
        return_vals = flux_dict

    return_vals = DataFrame(return_vals).transpose()
    return return_vals
