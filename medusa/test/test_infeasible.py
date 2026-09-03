"""Tests for how optimize_ensemble reports members that do not solve.

cobrapy's contract, verified against cobra 0.29: a non-optimal status is in
``has_primals``, so ``Model.optimize()`` with the default ``raise_error=False``
emits a UserWarning and then returns whatever primal values the solver left
behind. ``Reaction.flux`` happily reports those as though they were fluxes.

For a single model inspected by hand that is defensible. For an ensemble it is
not, for a reason specific to medusa: every member is solved against the same
``base_model``, so the stale primal an infeasible member picks up is usually
the *previously solved member's* solution. The failing member therefore
reports a plausible-looking flux distribution belonging to someone else, and
nothing in the returned table marks it.
"""

import warnings

import pytest
from cobra.exceptions import OptimizationError
from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.flux_analysis.flux_balance import optimize_ensemble


def _textbook(name, glc=-10.0, atpm_lb=8.39):
    model = load_model("textbook")
    model.id = name
    model.reactions.ATPM.lower_bound = atpm_lb
    model.reactions.EX_glc__D_e.lower_bound = glc
    return model


@pytest.fixture
def mixed_ensemble():
    """'fed' solves; 'starved' has no glucose but ATPM forced at 100."""
    return Ensemble(
        list_of_models=[_textbook("fed", glc=-10.0, atpm_lb=8.39),
                        _textbook("starved", glc=0.0, atpm_lb=100.0)],
        identifier="mixed")


FLUXES = ["Biomass_Ecoli_core", "ATPM"]


def test_infeasible_member_is_nan_not_a_stale_flux(mixed_ensemble):
    results = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                infeasible="nan")
    assert results.loc["starved"].isna().all(), (
        "the infeasible member reported fluxes: %s"
        % results.loc["starved"].to_dict())
    # A biomass flux is a growth rate; the stale primal used to come back
    # negative, which is not a physically meaningful value at all.
    assert results.loc["fed", "Biomass_Ecoli_core"] > 0.1


def test_feasible_members_are_unaffected(mixed_ensemble):
    with_failure = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                     infeasible="nan")
    only_good = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                  specific_models=["fed"])
    assert with_failure.loc["fed", "Biomass_Ecoli_core"] == pytest.approx(
        only_good.loc["fed", "Biomass_Ecoli_core"], abs=1e-6)


def test_member_status_is_reported(mixed_ensemble):
    results = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                infeasible="nan")
    statuses = results.attrs["member_status"]
    assert statuses["fed"] == "optimal"
    assert statuses["starved"] != "optimal"


def test_failed_members_are_locatable_from_the_frame(mixed_ensemble):
    results = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                infeasible="nan")
    failed = results.index[results.isna().all(axis=1)].tolist()
    assert failed == ["starved"]


def test_warn_policy_emits_one_warning_counting_the_failures(mixed_ensemble):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                          infeasible="warn")
    messages = [str(w.message) for w in caught
                if issubclass(w.category, UserWarning)
                and "did not solve to optimality" in str(w.message)]
    assert len(messages) == 1, (
        "expected exactly one aggregated warning, got %i: %s"
        % (len(messages), messages))
    assert "1 of 2" in messages[0]


def test_warn_policy_does_not_list_member_ids(mixed_ensemble):
    """An ensemble can hold thousands of members, so the warning counts them.

    The ids stay available on the result, which is where a caller who needs
    them should look.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        results = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                    infeasible="warn")
    message = [str(w.message) for w in caught
               if "did not solve to optimality" in str(w.message)][0]
    assert "starved" not in message
    assert "member_status" in message
    assert results.attrs["member_status"]["starved"] != "optimal"


def test_nan_policy_is_silent(mixed_ensemble):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                          infeasible="nan")
    assert not [w for w in caught
                if "did not solve to optimality" in str(w.message)]


def test_raise_policy_raises_naming_the_member(mixed_ensemble):
    with pytest.raises(OptimizationError, match="starved"):
        optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                          infeasible="raise")


def test_all_feasible_ensemble_does_not_warn_or_raise():
    ensemble = Ensemble(
        list_of_models=[_textbook("a", glc=-10.0), _textbook("b", glc=-5.0)],
        identifier="healthy")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        results = optimize_ensemble(ensemble, return_flux=FLUXES,
                                    infeasible="raise")
    assert not results.isna().any().any()
    assert not [w for w in caught
                if "did not solve to optimality" in str(w.message)]


def test_unknown_policy_is_rejected(mixed_ensemble):
    with pytest.raises(ValueError, match="infeasible must be one of"):
        optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                          infeasible="explode")


def test_row_order_follows_the_requested_members(mixed_ensemble):
    requested = ["starved", "fed"]
    results = optimize_ensemble(mixed_ensemble, return_flux=FLUXES,
                                specific_models=requested, infeasible="nan")
    assert results.index.tolist() == requested
