"""Tests for Ensemble construction and for set_state.

These cover four defects that all produced a wrong answer quietly:

- ``set_state`` wrote lower_bound and upper_bound as separate assignments, so
  switching an off reaction back on raised whenever it had a positive minimum
  flux.
- ``set_state`` applied a 'metabolites' state with ``combine=False``, which
  leaves metabolites a previous member introduced in place, so a member's
  stoichiometry depended on which members had been visited before it.
- The constructor handed live Reaction objects from ``list_of_models[1:]`` to
  ``add_reactions``, which takes ownership of them, so building an ensemble
  mutated the caller's input models.
- ``features=`` rewired the caller's Feature objects in place, so passing one
  feature list to two ensembles left the first pointing at the second.
"""

import pytest
from cobra import Metabolite, Reaction
from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.core.feature import Feature


def _textbook(name):
    model = load_model("textbook")
    model.id = name
    return model


# --------------------------------------------------------------------------
# set_state must not raise on an ordinary on/off ensemble.
# --------------------------------------------------------------------------

def _on_off_ensemble():
    """ATPM is on at (8.39, 1000) in one member and fully off in the other."""
    on = _textbook("on")
    off = _textbook("off")
    off.reactions.ATPM.bounds = (0.0, 0.0)
    return Ensemble(list_of_models=[on, off], identifier="onoff")


def test_set_state_survives_toggling_a_reaction_off_and_on():
    ensemble = _on_off_ensemble()
    atpm = ensemble.base_model.reactions.ATPM

    ensemble.set_state("off")
    assert atpm.bounds == (0.0, 0.0)

    # Writing lower_bound first would try to set 8.39 while the upper bound is
    # still 0, which cobra rejects.
    ensemble.set_state("on")
    assert atpm.bounds == (pytest.approx(8.39), pytest.approx(1000.0))


def test_set_state_is_order_independent_for_bounds():
    ensemble = _on_off_ensemble()
    atpm = ensemble.base_model.reactions.ATPM

    ensemble.set_state("on")
    first = atpm.bounds
    ensemble.set_state("off")
    ensemble.set_state("on")
    assert atpm.bounds == first


# --------------------------------------------------------------------------
# set_state must not leak 'metabolites' between members.
# --------------------------------------------------------------------------

def _biomass_ensemble():
    """Two biomass compositions; only one of them uses g6p_c."""
    model = _textbook("bof")
    return Ensemble.from_reaction_states(
        model,
        "Biomass_Ecoli_core",
        {"wt": {model.metabolites.atp_c: -59.81},
         "alt": {model.metabolites.atp_c: 0.0,
                 model.metabolites.g6p_c: -59.81}},
        component_attribute="metabolites",
        allow_new_metabolites=True,
        identifier="bof_ensemble")


def test_metabolite_state_does_not_leak_between_members():
    ensemble = _biomass_ensemble()
    biomass = ensemble.base_model.reactions.Biomass_Ecoli_core
    g6p = ensemble.base_model.metabolites.g6p_c

    ensemble.set_state("wt")
    wt_first = biomass.metabolites.get(g6p, 0.0)

    ensemble.set_state("alt")
    ensemble.set_state("wt")
    wt_second = biomass.metabolites.get(g6p, 0.0)

    assert wt_second == pytest.approx(wt_first), (
        "revisiting 'wt' gave a different g6p_c coefficient (%r then %r); "
        "the 'alt' member's metabolite was left behind"
        % (wt_first, wt_second))


def test_member_growth_rate_is_reproducible_across_visits():
    ensemble = _biomass_ensemble()

    ensemble.set_state("wt")
    first = ensemble.base_model.slim_optimize()

    ensemble.set_state("alt")
    ensemble.set_state("wt")
    second = ensemble.base_model.slim_optimize()

    assert second == pytest.approx(first, abs=1e-6), (
        "the same member reported two different growth rates depending on "
        "which member was visited before it")


# --------------------------------------------------------------------------
# Building an ensemble must not mutate the caller's models.
# --------------------------------------------------------------------------

def _model_with_extra_reaction(name):
    model = _textbook(name)
    extra = Reaction(id="EXTRA_RXN", lower_bound=-1000.0, upper_bound=1000.0)
    extra.add_metabolites({model.metabolites.atp_c: -1,
                           Metabolite("novel_c", compartment="c"): 1})
    model.add_reactions([extra])
    return model


def test_constructor_does_not_steal_reactions_from_input_models():
    first = _textbook("first")
    second = _model_with_extra_reaction("second")

    before_bounds = second.reactions.EXTRA_RXN.bounds
    ensemble = Ensemble(list_of_models=[first, second], identifier="steal")
    ensemble.set_state(ensemble.members[0].id)

    assert "EXTRA_RXN" in second.reactions, \
        "the reaction was removed from the caller's model"
    assert second.reactions.EXTRA_RXN.bounds == before_bounds, \
        "the caller's model was mutated by set_state"
    assert second.reactions.EXTRA_RXN.model is second, \
        "the caller's reaction was reassigned to the ensemble's base model"
    assert second.reactions.EXTRA_RXN is not \
        ensemble.base_model.reactions.EXTRA_RXN


def test_input_model_stays_structurally_consistent():
    first = _textbook("first")
    second = _model_with_extra_reaction("second")
    Ensemble(list_of_models=[first, second], identifier="consistent")

    reaction = second.reactions.EXTRA_RXN
    for metabolite in reaction.metabolites:
        assert metabolite is second.metabolites.get_by_id(metabolite.id), (
            "%s on the caller's reaction is not the caller's own metabolite"
            % metabolite.id)


# --------------------------------------------------------------------------
# Feature order must be reproducible.
# --------------------------------------------------------------------------

def test_feature_order_is_sorted_and_therefore_reproducible():
    first = _textbook("first")
    second = _textbook("second")
    second.reactions.PGI.lower_bound = 0.0
    second.reactions.PGK.lower_bound = 0.0
    second.reactions.ACALD.lower_bound = 0.0
    ensemble = Ensemble(list_of_models=[first, second], identifier="ordered")

    feature_ids = [feature.id for feature in ensemble.features]
    assert feature_ids == sorted(feature_ids), (
        "features are ordered by set iteration, which varies with string "
        "hash randomization across processes: %s" % feature_ids)


# --------------------------------------------------------------------------
# Prebuilt features must not be hijacked.
# --------------------------------------------------------------------------

def test_prebuilt_features_are_not_rebound_to_a_later_ensemble():
    first = _textbook("first")
    second = _textbook("second")
    feature = Feature(
        identifier="PGI_lower_bound",
        name="PGI",
        base_component=first.reactions.PGI,
        component_attribute="lower_bound",
        states={"a": -1000.0, "b": 0.0},
    )

    ensemble_one = Ensemble(list_of_models=[first], features=[feature],
                            identifier="one")
    ensemble_two = Ensemble(list_of_models=[second], features=[feature],
                            identifier="two")

    assert ensemble_one.features[0] is not ensemble_two.features[0]
    assert ensemble_one.features[0].ensemble is ensemble_one
    assert ensemble_two.features[0].ensemble is ensemble_two

    ensemble_one.set_state("b")
    assert first.reactions.PGI.lower_bound == 0.0, \
        "setting state on the first ensemble did not affect its own model"
    assert second.reactions.PGI.lower_bound == -1000.0, \
        "setting state on the first ensemble mutated the second's model"


def test_caller_feature_object_is_left_alone():
    model = _textbook("solo")
    feature = Feature(
        identifier="PGI_lower_bound",
        base_component=model.reactions.PGI,
        component_attribute="lower_bound",
        states={"a": -1000.0, "b": 0.0},
    )
    Ensemble(list_of_models=[model], features=[feature], identifier="one")
    assert feature.ensemble is None, \
        "the caller's Feature was rebound to the ensemble"


def test_default_list_of_models_is_not_shared():
    one = Ensemble(identifier="empty_one")
    two = Ensemble(identifier="empty_two")
    assert one.base_model is not two.base_model
