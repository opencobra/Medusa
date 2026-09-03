"""Tests for defects in medusa.reconstruct.expand.

Four separate problems, all silent or fatal rather than merely inaccurate:

- Continuous gapfilling, the default gapfill_type, died in ``add_pfba`` with a
  TypeError out of optlang's GLPK backend.
- ``gapfill_to_ensemble`` looked up cobrapy's Reaction objects as though they
  were ids, so it could never return an ensemble.
- ``exchange_prefix`` was accepted and then dropped, hardcoding 'EX_'.
- The intersection guard in ``_build_ensemble_from_gapfill_solutions`` treated
  a legitimately empty intersection as "first member", so reactions that were
  not in every solution got no Feature and stayed switched on for everyone.
"""

import pytest
from cobra.io import load_model

from medusa.reconstruct import expand


def _universal():
    return load_model("textbook")


def _knockout_model(name="target", knocked=("PGI", "PGK", "TPI")):
    model = load_model("textbook")
    model.id = name
    model.remove_reactions([model.reactions.get_by_id(r) for r in knocked])
    return model


# --------------------------------------------------------------------------
# The intersection guard.
# --------------------------------------------------------------------------

def _ensemble_from_solutions(solutions):
    return expand._build_ensemble_from_gapfill_solutions(
        _knockout_model(), solutions, universal=_universal())


def test_reactions_absent_from_some_solutions_become_features():
    """With no reaction common to every solution, the intersection is empty.

    Treating that as "first member" reset it to the current member's whole
    solution, so only the first reaction ever got a Feature.
    """
    ensemble = _ensemble_from_solutions([["PGI"], ["PGK"], ["PGK", "TPI"]])
    feature_ids = {feature.id for feature in ensemble.features}
    for reaction_id in ("PGI", "PGK", "TPI"):
        assert reaction_id + "_lower_bound" in feature_ids, (
            "%s varies across solutions but got no feature" % reaction_id)


def test_members_only_carry_their_own_gapfilled_reactions():
    ensemble = _ensemble_from_solutions([["PGI"], ["PGK"], ["PGK", "TPI"]])
    members = [member.id for member in ensemble.members]

    ensemble.set_state(members[0])
    assert ensemble.base_model.reactions.PGI.bounds != (0, 0)
    assert ensemble.base_model.reactions.PGK.bounds == (0, 0), (
        "this member's solution did not include PGK, but it is switched on")
    assert ensemble.base_model.reactions.TPI.bounds == (0, 0)

    ensemble.set_state(members[1])
    assert ensemble.base_model.reactions.PGI.bounds == (0, 0)
    assert ensemble.base_model.reactions.PGK.bounds != (0, 0)
    assert ensemble.base_model.reactions.TPI.bounds == (0, 0)


def test_a_reaction_in_every_solution_gets_no_feature():
    """The intersection still does its job when it is genuinely non-empty."""
    ensemble = _ensemble_from_solutions([["PGI", "PGK"], ["PGI", "TPI"]])
    feature_ids = {feature.id for feature in ensemble.features}
    assert "PGI_lower_bound" not in feature_ids
    assert "PGK_lower_bound" in feature_ids
    assert "TPI_lower_bound" in feature_ids


def test_feature_order_is_deterministic():
    ensemble = _ensemble_from_solutions([["PGI"], ["PGK"], ["PGK", "TPI"]])
    feature_ids = [feature.id for feature in ensemble.features]
    assert feature_ids == sorted(feature_ids)


# --------------------------------------------------------------------------
# gapfill_to_ensemble.
# --------------------------------------------------------------------------

def test_gapfill_to_ensemble_returns_an_ensemble():
    """cobrapy's GapFiller.fill returns Reaction objects, not ids."""
    ensemble = expand.gapfill_to_ensemble(
        _knockout_model(), iterations=2, universal=_universal(),
        lower_bound=0.05)
    assert len(ensemble.members) >= 1
    for member in ensemble.members:
        ensemble.set_state(member.id)
        assert ensemble.base_model.slim_optimize() > 0.04


# --------------------------------------------------------------------------
# exchange_prefix must reach the worker.
# --------------------------------------------------------------------------

def test_exchange_prefix_is_forwarded(monkeypatch):
    captured = {}

    def spy(*args, **kwargs):
        captured.update(kwargs)
        raise RuntimeError("stop here")

    monkeypatch.setattr(expand, "_continuous_iterative_binary_gapfill", spy)
    universal = _universal()
    medium = dict(load_model("textbook").medium)
    with pytest.raises(RuntimeError, match="stop here"):
        expand.iterative_gapfill_from_binary_phenotypes(
            _knockout_model(), universal, {"cond": medium}, 1,
            gapfill_type="continuous", exchange_prefix="R_EX_")
    assert captured.get("exchange_prefix") == "R_EX_", (
        "the caller's exchange_prefix was replaced with a hardcoded value; "
        "any non-ModelSEED namespace silently gapfills with no exchange "
        "handling at all")


# --------------------------------------------------------------------------
# The continuous path must run at all.
# --------------------------------------------------------------------------

def test_continuous_gapfill_completes():
    """Regression for the optlang/GLPK TypeError raised inside add_pfba."""
    universal = _universal()
    default_medium = dict(load_model("textbook").medium)
    conditions = {}
    for carbon in ("EX_glc__D_e", "EX_fru_e"):
        medium = dict(default_medium)
        medium.pop("EX_glc__D_e", None)
        medium[carbon] = 10.0
        conditions[carbon] = medium

    ensemble = expand.iterative_gapfill_from_binary_phenotypes(
        _knockout_model(), universal, conditions, 2,
        gapfill_type="continuous", lower_bound=0.05,
        inclusion_threshold=1e-10, exchange_reactions=False,
        demand_reactions=False, exchange_prefix="EX_")
    assert len(ensemble.members) >= 1
