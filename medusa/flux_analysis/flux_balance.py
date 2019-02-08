
from __future__ import absolute_import
import multiprocessing

from pandas import DataFrame
from random import sample

from builtins import dict, map
from functools import partial

from cobra import Reaction

from medusa.core.member import Member


def _optimize_ensemble(ensemble, return_flux, member_id, **kwargs):
    ensemble.set_state(member_id)
    ensemble.base_model.optimize(**kwargs)
    flux_dict = {rxn:ensemble.base_model.reactions.get_by_id(rxn).flux
                        for rxn in return_flux}
    return (member_id, flux_dict, ensemble.base_model.solver.status)


def _optimize_ensemble_worker(member_id):
    global _ensemble
    global _return_flux
    return _optimize_ensemble(_ensemble, _return_flux, member_id)


def _init_worker(ensemble, return_flux):
    global _ensemble
    global _return_flux
    _ensemble = ensemble
    _return_flux = return_flux


def optimize_ensemble(ensemble, return_flux = None, num_models = None,
                        specific_models = None, num_processes = None,
                        **kwargs):
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
    num_processes : int, optional
        An integer corresponding to the number of processes (i.e. cores) to
        use. Using more cores will speed up computation, but will have a larger
        memory footprint because the ensemble object must be temporarily
        copied for each additional core used. If None, one core is used.

    Returns
    -------
    pandas.DataFrame
        A dataframe in which each row (index) represents a model within the
        ensemble, and each column represents a reaction for which flux values
        are returned.
    '''
    if not num_models:
        num_models = len(ensemble.members)

    if not return_flux:
        return_flux = [rxn.id for rxn in ensemble.base_model.reactions]

    if isinstance(return_flux,str):
        return_flux = [return_flux]

    if isinstance(return_flux[0],Reaction):
        return_flux = [rxn.id for rxn in return_flux]

    if num_processes is None:
        num_processes = 1

    if specific_models:
        # If member objects were passed, convert to member.id
        if isinstance(specific_models[0],Member):
            model_list = [member.id for member in specific_models]
        else:
            model_list = specific_models
    elif len(ensemble.members) > num_models:
        model_list = sample([member.id for member in ensemble.members],
                            num_models)
    else:
        model_list = [member.id for member in ensemble.members]


    # Can't have fewer ensemble members than processes
    num_processes = min(num_processes, num_models)

    def extract_results(result_iter):
        return {member_id:flux_dict for
                (member_id, flux_dict, status) in result_iter}

    if num_processes > 1:
        # create worker
        worker = _optimize_ensemble_worker

        # determine chunk size
        chunk_size = num_models // num_processes

        pool = multiprocessing.Pool(
            num_processes,
            initializer = _init_worker,
            initargs = (ensemble, return_flux)
        )

        results = extract_results(pool.imap_unordered(
            worker,
            model_list,
            chunksize = chunk_size
        ))
        pool.close()
        pool.join()
    else:
        worker = _optimize_ensemble
        results = extract_results(map(
            partial(worker, ensemble, return_flux), model_list))

    return_vals = DataFrame(results).transpose()
    return return_vals
