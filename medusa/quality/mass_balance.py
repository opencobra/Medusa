"""Quality assessments for the genome-scale networks in an ensemble.

Every function here answers a question about whether the *network* is
physically sensible, independent of whether it fits any particular data:

- ``leak_test``: can a metabolite be produced from nothing?
- ``energy_generating_cycles``: can ATP be produced from nothing?
- ``mass_charge_balance``: is every internal reaction elementally and
  charge balanced?

All three run per ensemble member and return the member id on the index, so
results are directly comparable across members. Reported results always
distinguish "checked and clean" from "could not be checked": a member whose
linear program does not solve to optimality is reported as NaN, never as a
pass.
"""
from __future__ import absolute_import

import warnings

from numpy import nan
from pandas import DataFrame, Series

from cobra import Reaction
from cobra.util.solver import OPTIMAL

from medusa.flux_analysis._selection import resolve_member_ids

# Default tolerance for calling a flux nonzero. Solvers routinely return
# values a few orders of magnitude below this for a reaction that is truly
# blocked, so a tolerance well above solver noise is required to avoid
# reporting every model as leaky.
DEFAULT_TOLERANCE = 1e-6


def _boundary_reactions(model, exchange_prefix):
    """Return every reaction that can act as a free source or sink.

    The union of cobrapy's own boundary detection and an id-prefix match. The
    prefix is honoured because boundary detection depends on a reaction having
    exactly one metabolite, which is not true of every namespace's exchange
    conventions; the cobrapy set is included because demand and sink reactions
    are just as capable of feeding a spurious leak as an exchange is.
    """
    boundary = set(model.boundary)
    if exchange_prefix:
        prefixed = {rxn for rxn in model.reactions
                    if rxn.id.startswith(exchange_prefix)}
        if not prefixed:
            warnings.warn(
                "No reaction id starts with exchange_prefix %r. Falling back "
                "to cobrapy's boundary detection alone, which found %i "
                "reactions. If this model uses a different namespace (for "
                "example 'R_EX_'), pass the matching exchange_prefix, or "
                "pass exchange_prefix=None to silence this warning."
                % (exchange_prefix, len(boundary)), UserWarning)
        boundary |= prefixed
    return sorted(boundary, key=lambda rxn: rxn.id)


def _close_reactions(reactions):
    for reaction in reactions:
        reaction.bounds = (0.0, 0.0)


def _relax_forced_fluxes(model):
    """Widen every bound so that the all-zero flux distribution is feasible.

    A leak or energy-cycle test asks whether the network *can* produce
    something from nothing, which is only a meaningful question if it is also
    allowed to produce nothing. Maintenance reactions routinely carry a
    positive lower bound (ATPM is 8.39 in the E. coli core model), and once
    the boundary reactions are closed such a bound makes the problem
    infeasible rather than optimal at zero, so every metabolite would be
    reported as unmeasurable instead of clean. Relaxing forced fluxes toward
    zero is the standard treatment and does not weaken the test: a network
    that cannot create mass from nothing still cannot do so with wider bounds.
    """
    for reaction in model.reactions:
        lower, upper = reaction.bounds
        if lower > 0.0 or upper < 0.0:
            reaction.bounds = (min(lower, 0.0), max(upper, 0.0))


def leak_test(ensemble, metabolites_to_test=None, exchange_prefix='EX_',
              specific_models=None, num_models=None,
              tolerance=DEFAULT_TOLERANCE, verbose=False):
    '''
    Checks for leaky metabolites in members of the ensemble by opening and
    maximizing a demand reaction while every boundary reaction is closed.

    A metabolite is "leaky" if the network can produce it from nothing. That
    is always a reconstruction error: it means some reaction, or combination
    of reactions, creates mass out of nothing.

    Parameters
    ----------
    ensemble : medusa.core.ensemble.Ensemble
        The ensemble to test.
    metabolites_to_test : iterable of str or cobra.core.metabolite.Metabolite, optional
        The metabolites to test, as ids or Metabolite objects. If None, every
        metabolite in the base model is tested (default).
    exchange_prefix : str, optional
        Id prefix identifying exchange reactions, used in addition to
        cobrapy's own boundary detection. Defaults to 'EX_'. Pass None to rely
        on cobrapy's detection alone. Note that in earlier versions this
        argument was accepted and then never used, so no boundary reaction was
        ever actually closed and the test could not have detected a leak.
    specific_models : iterable of str or medusa.core.member.Member, optional
        The members to test. If None, all members are tested (default).
    num_models : int, optional
        Randomly sample this many members instead. Cannot be passed
        concurrently with specific_models.
    tolerance : float, optional
        Demand flux above which a metabolite counts as leaking.
    verbose : bool, optional
        Print progress per metabolite.

    Returns
    -------
    pandas.DataFrame
        Rows are member ids, columns are metabolite ids, values are booleans
        (True where the metabolite leaked). A member/metabolite pair whose
        linear program did not solve to optimality is NaN, which makes the
        column dtype object rather than bool for that column.

    Notes
    -----
    This solves one linear program per member per metabolite, so the cost is
    ``len(members) * len(metabolites_to_test)`` solves. Restrict
    metabolites_to_test or num_models on large reconstructions.
    '''
    model_list = resolve_member_ids(ensemble, specific_models, num_models)
    base_model = ensemble.base_model

    if metabolites_to_test is None:
        metabolite_ids = [met.id for met in base_model.metabolites]
    else:
        if isinstance(metabolites_to_test, str):
            metabolites_to_test = [metabolites_to_test]
        metabolite_ids = [met if isinstance(met, str) else met.id
                          for met in metabolites_to_test]
        unknown = [met_id for met_id in metabolite_ids
                   if met_id not in base_model.metabolites]
        if unknown:
            raise KeyError(
                "metabolites_to_test refers to metabolites that are not in "
                "the base model: %s" % ', '.join(repr(m) for m in unknown[:10]))

    if not metabolite_ids:
        raise ValueError("metabolites_to_test is empty; nothing to test.")

    results = {}
    # Everything below mutates base_model. The context manager reverts the
    # demand reactions, the closed bounds and the objective on exit, so the
    # ensemble is left exactly as it was found even if a solve raises.
    with base_model:
        demand_reactions = []
        for met_id in metabolite_ids:
            metabolite = base_model.metabolites.get_by_id(met_id)
            demand = Reaction(id='leak_DM_' + met_id)
            demand.lower_bound = 0.0
            demand.upper_bound = 0.0
            demand.add_metabolites({metabolite: -1})
            demand_reactions.append(demand)
        base_model.add_reactions(demand_reactions)

        for member_id in model_list:
            ensemble.set_state(member_id)
            # Boundaries are closed *after* set_state: a feature may vary the
            # bounds of an exchange reaction, which would otherwise reopen a
            # free source and mask every leak for that member.
            _close_reactions(_boundary_reactions(base_model, exchange_prefix))
            _relax_forced_fluxes(base_model)

            member_results = {}
            for demand in demand_reactions:
                met_id = demand.id.split('leak_DM_', 1)[1]
                if verbose:
                    print('checking leak for ' + met_id +
                          ' in member ' + member_id)
                demand.upper_bound = 1000.0
                base_model.objective = demand
                value = base_model.slim_optimize()
                if base_model.solver.status == OPTIMAL:
                    member_results[met_id] = bool(value > tolerance)
                else:
                    member_results[met_id] = nan
                demand.upper_bound = 0.0
            results[member_id] = member_results

    return DataFrame(results).transpose().reindex(
        index=model_list, columns=metabolite_ids)


def energy_generating_cycles(ensemble, atp_hydrolysis_id='ATPM',
                             exchange_prefix='EX_', specific_models=None,
                             num_models=None, tolerance=DEFAULT_TOLERANCE):
    '''
    Checks whether members can generate energy from nothing.

    Closes every boundary reaction and maximizes an ATP hydrolysis reaction.
    With no inputs available the maximum must be zero; anything above
    `tolerance` means the network contains a thermodynamically impossible
    energy-generating cycle. This is the single most common serious defect in
    an automatically generated reconstruction, and unlike a metabolite leak it
    inflates growth yields on every medium rather than only on some.

    Parameters
    ----------
    ensemble : medusa.core.ensemble.Ensemble
    atp_hydrolysis_id : str, optional
        Id of the ATP hydrolysis / non-growth-associated maintenance reaction
        to maximize. Defaults to 'ATPM'.
    exchange_prefix : str, optional
        See leak_test.
    specific_models : iterable of str or medusa.core.member.Member, optional
    num_models : int, optional
    tolerance : float, optional
        Flux above which the cycle is reported.

    Returns
    -------
    pandas.DataFrame
        Indexed by member id, with columns 'max_flux' (float, NaN where the
        member did not solve to optimality) and 'has_cycle' (object; True,
        False, or NaN where max_flux is NaN).
    '''
    model_list = resolve_member_ids(ensemble, specific_models, num_models)
    base_model = ensemble.base_model

    if atp_hydrolysis_id not in base_model.reactions:
        raise KeyError(
            "atp_hydrolysis_id %r is not a reaction in the base model. Pass "
            "the id of this model's ATP hydrolysis or non-growth-associated "
            "maintenance reaction." % atp_hydrolysis_id)

    max_flux = {}
    with base_model:
        atp_reaction = base_model.reactions.get_by_id(atp_hydrolysis_id)
        for member_id in model_list:
            ensemble.set_state(member_id)
            _close_reactions(_boundary_reactions(base_model, exchange_prefix))
            _relax_forced_fluxes(base_model)
            # _relax_forced_fluxes drops the maintenance reaction's positive
            # lower bound to zero; make sure it can still carry flux upward,
            # since that is the quantity being maximized.
            if atp_reaction.upper_bound <= 0.0:
                atp_reaction.upper_bound = 1000.0
            base_model.objective = atp_reaction
            value = base_model.slim_optimize()
            if base_model.solver.status == OPTIMAL:
                max_flux[member_id] = value
            else:
                max_flux[member_id] = nan

    flux_series = Series(max_flux, dtype=float).reindex(model_list)
    has_cycle = Series(
        [nan if value != value else bool(value > tolerance)
         for value in flux_series],
        index=flux_series.index, dtype=object)
    return DataFrame({'max_flux': flux_series, 'has_cycle': has_cycle})


def mass_charge_balance(ensemble, specific_models=None, num_models=None,
                        skip_boundary=True):
    '''
    Reports internal reactions that are not elementally and charge balanced.

    Parameters
    ----------
    ensemble : medusa.core.ensemble.Ensemble
    specific_models : iterable of str or medusa.core.member.Member, optional
    num_models : int, optional
    skip_boundary : bool, optional
        Skip exchange, demand and sink reactions, which are unbalanced by
        construction. True by default; setting it False is rarely useful.

    Returns
    -------
    dict
        ``{member_id: {reaction_id: imbalance}}``, where `imbalance` is the
        dict returned by cobrapy's ``Reaction.check_mass_balance`` (element
        symbols and 'charge' mapped to the size of the imbalance). Reactions
        that balance are omitted, so a member that maps to an empty dict is
        fully balanced.

    Notes
    -----
    Mass balance is a property of stoichiometry and of metabolite formulas,
    neither of which varies across members unless a Feature targets the
    ``metabolites`` component_attribute. For ensembles that vary only reaction
    bounds, every member therefore returns an identical result, and running
    this on a single member is sufficient. It is computed per member anyway so
    that ensembles varying biomass or other reaction stoichiometry are handled
    correctly.

    A reaction whose metabolites are missing formulas cannot be checked. Those
    are reported under the key ``'medusa_uncheckable'`` with the exception
    text as the value, so they are never silently counted as balanced.
    '''
    model_list = resolve_member_ids(ensemble, specific_models, num_models)
    base_model = ensemble.base_model

    imbalances = {}
    with base_model:
        for member_id in model_list:
            ensemble.set_state(member_id)
            member_imbalances = {}
            for reaction in base_model.reactions:
                if skip_boundary and reaction.boundary:
                    continue
                try:
                    imbalance = reaction.check_mass_balance()
                except Exception as error:
                    member_imbalances[reaction.id] = {
                        'medusa_uncheckable': str(error)}
                    continue
                if imbalance:
                    member_imbalances[reaction.id] = imbalance
            imbalances[member_id] = member_imbalances

    return imbalances
