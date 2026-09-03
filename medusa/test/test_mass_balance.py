"""Tests for medusa.quality.mass_balance.

The module previously could not run at all: it had no imports, referenced
``self`` twice inside a module-level function, called ``cobra.Reaction``
without importing cobra, called ``ensemble.optimize_ensemble`` which is not a
method of Ensemble, and indexed the result by reaction id where member ids
were meant. It also accepted an ``exchange_prefix`` argument and never used
it, so no boundary reaction was ever closed and the leak test could not have
detected a leak even had it run.

Each check below is tested both ways: it must stay quiet on a curated model
and it must fire on a model with a deliberately planted defect.
"""

import warnings

import pytest
from cobra import Reaction
from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.quality.mass_balance import (leak_test, energy_generating_cycles,
                                         mass_charge_balance)


def _textbook(name, glc=-10.0):
    model = load_model("textbook")
    model.id = name
    model.reactions.EX_glc__D_e.lower_bound = glc
    return model


def _leaky(name, glc=-10.0):
    """Net creation of g6p_c from nothing, with the co-product balanced.

    FAKE_A:  -> g6p_c + pyr_c
    FAKE_B:  pyr_c -> g6p_c

    Both carry two metabolites, so neither is a cobra boundary reaction and
    neither gets closed by the test. At equal flux the pair nets g6p_c out of
    nothing while consuming its own co-product.
    """
    model = _textbook(name, glc=glc)
    a = Reaction(id="FAKE_A", lower_bound=0.0, upper_bound=1000.0)
    a.add_metabolites({model.metabolites.g6p_c: 1,
                       model.metabolites.pyr_c: 1})
    b = Reaction(id="FAKE_B", lower_bound=0.0, upper_bound=1000.0)
    b.add_metabolites({model.metabolites.pyr_c: -1,
                       model.metabolites.g6p_c: 1})
    model.add_reactions([a, b])
    return model


def _free_atp(name, glc=-10.0):
    """A source of atp_c + h2o_c and a sink for adp_c + h_c + pi_c, which
    together let ATPM carry flux with every boundary reaction closed."""
    model = _textbook(name, glc=glc)
    source = Reaction(id="FAKE_ATP_SRC", lower_bound=0.0, upper_bound=1000.0)
    source.add_metabolites({model.metabolites.atp_c: 1,
                            model.metabolites.h2o_c: 1})
    sink = Reaction(id="FAKE_SINK", lower_bound=0.0, upper_bound=1000.0)
    sink.add_metabolites({model.metabolites.adp_c: -1,
                          model.metabolites.h_c: -1,
                          model.metabolites.pi_c: -1})
    model.add_reactions([source, sink])
    return model


def _ensemble(builder, identifier):
    return Ensemble(
        list_of_models=[builder("one", glc=-10.0), builder("two", glc=-5.0)],
        identifier=identifier)


SOME_METABOLITES = ["g6p_c", "pyr_c", "atp_c"]


# --------------------------------------------------------------------------
# leak_test
# --------------------------------------------------------------------------

def test_curated_model_has_no_leaks():
    ensemble = _ensemble(_textbook, "clean")
    results = leak_test(ensemble, metabolites_to_test=SOME_METABOLITES)
    assert results.shape == (2, 3)
    assert not results.to_numpy().any(), (
        "the E. coli core model should not leak: %s" % results.to_dict())


def test_leak_is_detected_and_localized():
    ensemble = _ensemble(_leaky, "leaky")
    results = leak_test(ensemble, metabolites_to_test=SOME_METABOLITES)
    for member_id in results.index:
        assert results.loc[member_id, "g6p_c"] is True or \
            bool(results.loc[member_id, "g6p_c"]), \
            "the planted g6p_c leak was not detected for %s" % member_id
    # The check must be specific, not just alarmed.
    assert not bool(results["atp_c"].any())
    assert not bool(results["pyr_c"].any())


def test_leak_test_solves_rather_than_reporting_nan():
    """A forced maintenance flux must not make every metabolite unmeasurable.

    ATPM carries a lower bound of 8.39 in this model. Closing the boundary
    reactions without relaxing that bound makes the problem infeasible, which
    would report every metabolite as NaN instead of clean.
    """
    ensemble = _ensemble(_textbook, "clean")
    results = leak_test(ensemble, metabolites_to_test=SOME_METABOLITES)
    assert not results.isna().to_numpy().any()


def test_leak_test_accepts_a_single_metabolite_id():
    ensemble = _ensemble(_textbook, "clean")
    results = leak_test(ensemble, metabolites_to_test="atp_c")
    assert results.shape == (2, 1)


def test_leak_test_unknown_metabolite_raises():
    ensemble = _ensemble(_textbook, "clean")
    with pytest.raises(KeyError, match="not in the base model"):
        leak_test(ensemble, metabolites_to_test=["not_a_metabolite"])


def test_leak_test_warns_when_exchange_prefix_matches_nothing():
    ensemble = _ensemble(_textbook, "clean")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        leak_test(ensemble, metabolites_to_test=["atp_c"],
                  exchange_prefix="R_EX_")
    assert any("exchange_prefix" in str(w.message) for w in caught), (
        "a prefix that matches no reaction should be reported, since it "
        "silently disables the caller's intended boundary selection")


def test_leak_test_restores_the_base_model():
    ensemble = _ensemble(_textbook, "clean")
    before = {rxn.id: rxn.bounds for rxn in ensemble.base_model.reactions}
    n_reactions = len(ensemble.base_model.reactions)
    growth = ensemble.base_model.slim_optimize()

    leak_test(ensemble, metabolites_to_test=SOME_METABOLITES)

    assert len(ensemble.base_model.reactions) == n_reactions, \
        "demand reactions were left behind"
    assert {rxn.id: rxn.bounds
            for rxn in ensemble.base_model.reactions} == before
    assert ensemble.base_model.slim_optimize() == pytest.approx(growth,
                                                                abs=1e-6)


# --------------------------------------------------------------------------
# energy_generating_cycles
# --------------------------------------------------------------------------

def test_curated_model_has_no_energy_generating_cycle():
    ensemble = _ensemble(_textbook, "clean")
    results = energy_generating_cycles(ensemble)
    assert list(results.columns) == ["max_flux", "has_cycle"]
    assert not results["has_cycle"].any()
    assert results["max_flux"].abs().max() < 1e-6


def test_energy_generating_cycle_is_detected():
    ensemble = _ensemble(_free_atp, "egc")
    results = energy_generating_cycles(ensemble)
    assert results["has_cycle"].all(), (
        "free ATP production was not detected: %s" % results.to_dict())
    assert results["max_flux"].min() > 1.0


def test_energy_generating_cycles_unknown_reaction_raises():
    ensemble = _ensemble(_textbook, "clean")
    with pytest.raises(KeyError, match="atp_hydrolysis_id"):
        energy_generating_cycles(ensemble, atp_hydrolysis_id="NOPE")


def test_energy_generating_cycles_restores_the_base_model():
    ensemble = _ensemble(_textbook, "clean")
    before = {rxn.id: rxn.bounds for rxn in ensemble.base_model.reactions}
    energy_generating_cycles(ensemble)
    assert {rxn.id: rxn.bounds
            for rxn in ensemble.base_model.reactions} == before


# --------------------------------------------------------------------------
# mass_charge_balance
# --------------------------------------------------------------------------

def test_mass_charge_balance_reports_per_member():
    ensemble = _ensemble(_textbook, "clean")
    results = mass_charge_balance(ensemble)
    assert set(results) == {member.id for member in ensemble.members}


def test_only_biomass_is_unbalanced_in_the_curated_model():
    """Biomass reactions are unbalanced by convention; nothing else should be."""
    ensemble = _ensemble(_textbook, "clean")
    results = mass_charge_balance(ensemble, specific_models=["one"])
    assert list(results["one"]) == ["Biomass_Ecoli_core"]


def test_planted_imbalance_is_reported():
    ensemble = _ensemble(_leaky, "leaky")
    results = mass_charge_balance(ensemble, specific_models=["one"])
    assert "FAKE_A" in results["one"]
    assert "FAKE_B" in results["one"]


def test_boundary_reactions_are_skipped_by_default():
    ensemble = _ensemble(_textbook, "clean")
    results = mass_charge_balance(ensemble, specific_models=["one"])
    exchanges = {rxn.id for rxn in ensemble.base_model.boundary}
    assert not (set(results["one"]) & exchanges)
