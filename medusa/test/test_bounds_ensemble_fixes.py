"""Tests for defects in medusa.boundsEnsemble.

- The ensemble used the caller's model as its base model without copying, so
  set_state permanently rewrote the caller's reactions.
- A reaction whose requested bounds were constant across members but different
  from the base model got no feature and no change, so the request vanished.
- Member ids were taken from whichever reaction's DataFrame came first, so a
  member present only in a later frame was dropped without a word.
- reversibility=True was implemented as a test for a bound of exactly zero, so
  an irreversible reaction with a nonzero minimum flux was still given
  negative lower bounds.
- The full-factorial "active" option scaled each bound by its own sign, which
  turns a (8.39, 1000) reaction into (1000, 1000): forced flux, not active.
"""

import pandas as pd
import pytest
from cobra.io import load_model

from medusa.boundsEnsemble import (boundsEnsemble, _setBounds,
                                   _setBoundsFullFactorial, _setBoundsRandom)

# ATPM in the E. coli core model is irreversible with a positive minimum flux,
# which is what makes it the interesting case for all of these.
ATPM_BOUNDS = (8.39, 1000.0)


def _textbook():
    return load_model("textbook")


def _frame(ids, lower, upper):
    return pd.DataFrame({"lower_bound": lower, "upper_bound": upper},
                        index=ids)


# --------------------------------------------------------------------------
# The caller's model must survive.
# --------------------------------------------------------------------------

def test_caller_model_is_not_mutated():
    model = _textbook()
    assert model.reactions.ATPM.bounds == ATPM_BOUNDS

    bounds = _setBounds(model, ["ATPM"], method="onOff")
    ensemble = boundsEnsemble(model, bounds)
    for member in ensemble.members:
        ensemble.set_state(member.id)

    assert model.reactions.ATPM.bounds == ATPM_BOUNDS, \
        "the caller's model was left holding a member's state"
    assert ensemble.base_model is not model


def test_on_off_ensemble_can_be_switched_back_on():
    model = _textbook()
    ensemble = boundsEnsemble(model, _setBounds(model, ["ATPM"],
                                                method="onOff"))
    seen = set()
    for member in ensemble.members:
        ensemble.set_state(member.id)
        seen.add(ensemble.base_model.reactions.ATPM.bounds)
    assert (0.0, 0.0) in seen
    assert ATPM_BOUNDS in seen


# --------------------------------------------------------------------------
# Constant-but-different bounds must be applied.
# --------------------------------------------------------------------------

def test_constant_request_that_differs_from_base_is_applied():
    model = _textbook()
    assert model.reactions.PGI.bounds == (-1000.0, 1000.0)
    bounds = {
        # constant across members, but not the base model's value
        "PGI": _frame(["m1", "m2"], [-500.0, -500.0], [500.0, 500.0]),
        # genuinely varying, so this one becomes a feature
        "PGK": _frame(["m1", "m2"], [-1000.0, 0.0], [1000.0, 1000.0]),
    }
    ensemble = boundsEnsemble(model, bounds)

    for member in ensemble.members:
        ensemble.set_state(member.id)
        assert ensemble.base_model.reactions.PGI.bounds == (-500.0, 500.0), \
            "the requested constant bound for PGI was silently ignored"


def test_constant_request_does_not_create_a_feature():
    model = _textbook()
    bounds = {
        "PGI": _frame(["m1", "m2"], [-500.0, -500.0], [500.0, 500.0]),
        "PGK": _frame(["m1", "m2"], [-1000.0, 0.0], [1000.0, 1000.0]),
    }
    ensemble = boundsEnsemble(model, bounds)
    feature_ids = {feature.id for feature in ensemble.features}
    assert "PGI_lower_bound" not in feature_ids
    assert "PGK_lower_bound" in feature_ids


# --------------------------------------------------------------------------
# Member selection must not silently drop anyone.
# --------------------------------------------------------------------------

def test_inconsistent_member_index_raises():
    model = _textbook()
    bounds = {
        "PGI": _frame(["s1", "s2"], [-1000.0, 0.0], [1000.0, 1000.0]),
        "PGK": _frame(["s1", "s2", "s3"], [-1000.0, 0.0, -250.0],
                      [1000.0, 1000.0, 1000.0]),
    }
    with pytest.raises(ValueError, match="same member ids"):
        boundsEnsemble(model, bounds)


def test_missing_column_raises():
    model = _textbook()
    bounds = {"PGI": pd.DataFrame({"lower_bound": [-1000.0, 0.0]},
                                  index=["s1", "s2"])}
    with pytest.raises(ValueError, match="missing required column"):
        boundsEnsemble(model, bounds)


def test_empty_boundsdict_raises():
    with pytest.raises(ValueError, match="boundsDict is empty"):
        boundsEnsemble(_textbook(), {})


# --------------------------------------------------------------------------
# reversibility must mean reversibility.
# --------------------------------------------------------------------------

def test_random_respects_irreversibility_of_a_nonzero_lower_bound():
    model = _textbook()
    assert not model.reactions.ATPM.reversibility
    bounds = _setBoundsRandom(model, ["ATPM"], bound=1000,
                              reversibility=True, n_models=8)
    lower = bounds["ATPM"]["lower_bound"]
    assert (lower >= 0).all(), (
        "an irreversible reaction was given negative lower bounds: %s"
        % lower.tolist())


def test_random_still_allows_reversibility_when_not_requested():
    model = _textbook()
    bounds = _setBoundsRandom(model, ["ATPM"], bound=1000,
                              reversibility=False, n_models=8)
    assert (bounds["ATPM"]["lower_bound"] < 0).any()


def test_full_factorial_active_option_does_not_force_flux():
    model = _textbook()
    bounds = _setBoundsFullFactorial(model, ["ATPM"], bound=1000,
                                     reversibility=True)
    pairs = set(zip(bounds["ATPM"]["lower_bound"],
                    bounds["ATPM"]["upper_bound"]))
    assert (1000.0, 1000.0) not in pairs, (
        "the 'active' option forces flux at exactly the bound rather than "
        "allowing anything up to it")
    assert (0.0, 0.0) in pairs
    assert (0.0, 1000.0) in pairs
