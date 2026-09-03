
from __future__ import absolute_import
import multiprocessing
import warnings

from pandas import DataFrame
from numpy import nan
from random import sample

from builtins import dict, map
from functools import partial

from cobra import Reaction
from cobra.exceptions import OptimizationError
from cobra.util.solver import OPTIMAL

from medusa.core.member import Member

# Policies for what to do with an ensemble member whose solve did not return
# an optimal status. See optimize_ensemble for the rationale.
INFEASIBLE_POLICIES = ('warn', 'nan', 'raise')

# Number of member ids named in an aggregated message before truncating.
_MAX_IDS_IN_MESSAGE = 10


def _summarize_failures(statuses):
    """Format the non-optimal members for a message, capping how many list."""
    failed = {member_id: status for member_id, status in statuses.items()
              if status != OPTIMAL}
    summary = ', '.join(
        '%s (%s)' % (member_id, failed[member_id])
        for member_id in list(failed)[:_MAX_IDS_IN_MESSAGE])
    if len(failed) > _MAX_IDS_IN_MESSAGE:
        summary += ' ... (%i total)' % len(failed)
    return failed, summary


def _apply_infeasible_policy(statuses, infeasible):
    """Warn or raise about non-optimal members according to `infeasible`.

    Called once for the whole ensemble rather than once per member. cobrapy
    warns per solve, which is the right granularity when a human is looking
    at a single model but is unusable across hundreds of members.
    """
    failed, summary = _summarize_failures(statuses)
    if not failed:
        return
    if infeasible == 'raise':
        raise OptimizationError(
            "%i of %i ensemble members did not solve to optimality: %s. "
            "Pass infeasible='warn' or infeasible='nan' to return NaN for "
            "these members instead of raising."
            % (len(failed), len(statuses), summary))
    if infeasible == 'warn':
        warnings.warn(
            "%i of %i ensemble members did not solve to optimality and their "
            "fluxes were set to NaN: %s. Member statuses are available on the "
            "returned DataFrame as .attrs['member_status']."
            % (len(failed), len(statuses), summary),
            UserWarning)


def _optimize_ensemble(ensemble, return_flux, member_id, **kwargs):

    ensemble.set_state(member_id)
    # cobrapy warns once per non-optimal solve and then returns whatever
    # primal the solver left behind. Because every member is solved against
    # the same base_model, that stale primal is usually the *previously
    # solved member's* solution, so an infeasible member silently inherits a
    # neighbour's fluxes. The status is returned alongside the fluxes and the
    # values are discarded below; the per-solve warning is suppressed here and
    # re-emitted once, aggregated, by _apply_infeasible_policy.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore', message='Solver status is', category=UserWarning)
        ensemble.base_model.optimize(**kwargs)

    status = ensemble.base_model.solver.status
    if status == OPTIMAL:
        flux_dict = {rxn_id: ensemble.base_model.reactions.get_by_id(
            rxn_id).flux for rxn_id in return_flux}
    else:
        flux_dict = {rxn_id: nan for rxn_id in return_flux}
    return (member_id, flux_dict, status)


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
                        infeasible = 'warn', **kwargs):
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
    infeasible : {'warn', 'nan', 'raise'}, optional
        What to do about members whose solve does not return an optimal
        status. In every case the fluxes reported for such a member are NaN;
        this argument controls only how loudly that is announced.

        - 'warn' (default): emit a single UserWarning naming the affected
          members once every member has been solved.
        - 'nan': stay silent. Appropriate for large ensembles in which some
          members are expected to be infeasible under the tested condition.
        - 'raise': raise cobra.exceptions.OptimizationError naming the
          affected members.

        This deliberately differs from cobrapy, which warns and then returns
        the solver's stale primal values as though they were fluxes. That is
        defensible for a single model being inspected by hand, but not for an
        ensemble: the per-solve warning is lost among hundreds of members,
        nothing in the returned table distinguishes a real flux from a stale
        one, and because all members share one base_model the stale values are
        typically the previously-solved member's solution rather than anything
        to do with the member that failed.

    Returns
    -------
    pandas.DataFrame
        A dataframe in which each row (index) represents a model within the
        ensemble, and each column represents a reaction for which flux values
        are returned. Rows are ordered to match the members that were
        requested, regardless of num_processes. Members that did not solve to
        optimality are all-NaN rows, so they can be located with
        ``results.isna().all(axis=1)``.

        The solver status for every member is attached as
        ``results.attrs['member_status']``, a dict of {member_id: status}.
        Note that pandas does not preserve ``.attrs`` through most operations,
        so read it before manipulating the frame.
    '''
    if infeasible not in INFEASIBLE_POLICIES:
        raise ValueError(
            "infeasible must be one of %s; got %r."
            % (', '.join(repr(p) for p in INFEASIBLE_POLICIES), infeasible))

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

    if not model_list:
        raise ValueError(
            "No ensemble members were selected for optimization; the ensemble "
            "has %i members." % len(ensemble.members))

    # Chunking and the process cap must follow the number of members actually
    # being solved, not the num_models default of len(ensemble.members), which
    # is still that full count when specific_models selects a smaller subset.
    num_models = len(model_list)

    # Can't have fewer ensemble members than processes
    num_processes = min(num_processes, num_models)

    def extract_results(result_iter):
        fluxes = {}
        statuses = {}
        for (member_id, flux_dict, status) in result_iter:
            fluxes[member_id] = flux_dict
            statuses[member_id] = status
        return fluxes, statuses

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

        results, statuses = extract_results(pool.imap_unordered(
            worker,
            model_list,
            chunksize = chunk_size
        ))
        pool.close()
        pool.join()
    else:
        worker = _optimize_ensemble
        # set_state mutates base_model. The sibling functions in
        # medusa.flux_analysis wrap their loop in the model's context manager
        # so that the ensemble is left as it was found. Without it the
        # base_model retains the last-solved member's bounds, which silently
        # changes the result of the *next* simulation run on the same ensemble
        # and, when that last member was infeasible, leaves the base_model
        # unsolvable altogether.
        with ensemble.base_model:
            results, statuses = extract_results(map(
                partial(worker, ensemble, return_flux), model_list))

    _apply_infeasible_policy(statuses, infeasible)

    return_vals = DataFrame(results).transpose()
    # imap_unordered yields members in completion order, so without this the
    # row order of the result depends on num_processes and on solver timing.
    return_vals = return_vals.reindex(model_list)
    return_vals.attrs['member_status'] = statuses
    return return_vals
