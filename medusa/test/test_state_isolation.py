"""Tests that a simulation leaves the ensemble as it found it.

Every entry point in medusa.flux_analysis works by mutating a single shared
``ensemble.base_model`` through ``Ensemble.set_state`` and then solving. If a
function does not undo those mutations, the ensemble is left holding the
last-solved member's state rather than its own. The next thing done with that
ensemble then silently answers a different question.

``optimize_ensemble`` had exactly this defect: its two siblings wrapped their
loop in ``with ensemble.base_model:`` and it did not, so it left the base
model in the final member's state. When that member was infeasible, the base
model was left unsolvable and every later use of it returned NaN.

The tests here are deliberately written as "run one simulation, then run a
second and check it against a value known in advance", because that is the
form in which a user meets the bug.
"""

import pytest
from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.flux_analysis.flux_balance import optimize_ensemble
from medusa.flux_analysis.variability import ensemble_fva
from medusa.flux_analysis.deletion import (ensemble_single_reaction_deletion,
                                           ensemble_single_gene_deletion)

TOL = 1e-6


def _textbook(name, glc=-10.0, atpm_lb=8.39):
    model = load_model("textbook")
    model.id = name
    model.reactions.ATPM.lower_bound = atpm_lb
    model.reactions.EX_glc__D_e.lower_bound = glc
    return model


def healthy_ensemble():
    """Two feasible members differing in how much glucose they may take up."""
    return Ensemble(
        list_of_models=[_textbook("plenty", glc=-10.0),
                        _textbook("scarce", glc=-5.0)],
        identifier="healthy")


def ensemble_with_infeasible_member():
    """The second member cannot solve: no glucose, but ATPM forced at 100.

    Ordering matters. 'starved' is added second so that it is the state the
    base model is left holding if a simulation fails to clean up after itself.
    """
    return Ensemble(
        list_of_models=[_textbook("fed", glc=-10.0, atpm_lb=8.39),
                        _textbook("starved", glc=0.0, atpm_lb=100.0)],
        identifier="mixed")


def bounds_snapshot(model):
    return {reaction.id: reaction.bounds for reaction in model.reactions}


# --------------------------------------------------------------------------
# The base model must survive a simulation unchanged.
# --------------------------------------------------------------------------

def test_optimize_ensemble_restores_base_model_bounds():
    ensemble = healthy_ensemble()
    before = bounds_snapshot(ensemble.base_model)
    optimize_ensemble(ensemble, return_flux=["Biomass_Ecoli_core"])
    after = bounds_snapshot(ensemble.base_model)
    changed = {rxn_id: (before[rxn_id], after[rxn_id])
               for rxn_id in before if before[rxn_id] != after[rxn_id]}
    assert changed == {}, (
        "optimize_ensemble left the base model in a member's state: %s"
        % changed)


def test_optimize_ensemble_restores_base_model_objective_value():
    ensemble = healthy_ensemble()
    before = ensemble.base_model.slim_optimize()
    optimize_ensemble(ensemble, return_flux=["Biomass_Ecoli_core"])
    after = ensemble.base_model.slim_optimize()
    assert before == pytest.approx(after, abs=TOL)


def test_second_simulation_gives_known_result_only_if_base_model_was_reset():
    """The regression test for the corruption, in the form a user meets it.

    The base model's growth rate is recorded before anything is simulated.
    An ensemble-wide FBA is then run, which internally walks the base model
    through every member's state, ending on the infeasible one. Re-running FBA
    on the base model afterwards can only reproduce the recorded growth rate
    if the simulation put the base model back; if it did not, the base model
    is still starved with ATPM forced at 100 and the solve returns NaN.
    """
    ensemble = ensemble_with_infeasible_member()
    known_growth = ensemble.base_model.slim_optimize()
    assert known_growth > 0.1, "precondition: the base model should grow"

    optimize_ensemble(ensemble, return_flux=["Biomass_Ecoli_core"],
                      infeasible="nan")

    second_growth = ensemble.base_model.slim_optimize()
    assert second_growth == pytest.approx(known_growth, abs=TOL), (
        "the base model no longer reproduces its own growth rate after a "
        "simulation; it was left holding a member's state")


def test_repeated_simulations_are_identical():
    """Running the same simulation twice on one ensemble must not drift."""
    ensemble = healthy_ensemble()
    first = optimize_ensemble(ensemble, return_flux=["Biomass_Ecoli_core"])
    second = optimize_ensemble(ensemble, return_flux=["Biomass_Ecoli_core"])
    assert first["Biomass_Ecoli_core"].tolist() == pytest.approx(
        second["Biomass_Ecoli_core"].tolist(), abs=TOL)


def test_simulation_on_used_ensemble_matches_fresh_ensemble():
    """A second, different simulation must not inherit the first one's state.

    The comparison is against the same call made on an ensemble that has never
    been simulated, which is the value the caller is entitled to expect.
    """
    used = ensemble_with_infeasible_member()
    optimize_ensemble(used, return_flux=["Biomass_Ecoli_core"],
                      infeasible="nan")
    after_use = optimize_ensemble(used, return_flux=["ATPM"],
                                  specific_models=["fed"])

    fresh = ensemble_with_infeasible_member()
    on_fresh = optimize_ensemble(fresh, return_flux=["ATPM"],
                                 specific_models=["fed"])

    assert after_use.loc["fed", "ATPM"] == pytest.approx(
        on_fresh.loc["fed", "ATPM"], abs=TOL)


def test_fba_then_fva_matches_fva_on_fresh_ensemble():
    """Cross-function: FVA after an FBA must equal FVA on a fresh ensemble."""
    used = healthy_ensemble()
    optimize_ensemble(used, return_flux=["Biomass_Ecoli_core"])
    after_use = ensemble_fva(used, reaction_list=["PGI"],
                             specific_models=["plenty"])

    fresh = healthy_ensemble()
    on_fresh = ensemble_fva(fresh, reaction_list=["PGI"],
                            specific_models=["plenty"])

    assert after_use["PGI"].tolist() == pytest.approx(
        on_fresh["PGI"].tolist(), abs=1e-5)


# --------------------------------------------------------------------------
# The siblings already used the context manager. These guard that.
# --------------------------------------------------------------------------

def test_ensemble_fva_restores_base_model():
    ensemble = healthy_ensemble()
    before = bounds_snapshot(ensemble.base_model)
    ensemble_fva(ensemble, reaction_list=["PGI"])
    assert bounds_snapshot(ensemble.base_model) == before


def test_reaction_deletion_restores_base_model():
    ensemble = healthy_ensemble()
    before = bounds_snapshot(ensemble.base_model)
    ensemble_single_reaction_deletion(ensemble, specific_models=["plenty"])
    assert bounds_snapshot(ensemble.base_model) == before


def test_gene_deletion_restores_base_model():
    ensemble = healthy_ensemble()
    before = bounds_snapshot(ensemble.base_model)
    gene_ids = [gene.id for gene in ensemble.base_model.genes[:3]]
    ensemble_single_gene_deletion(ensemble, specific_models=["plenty"],
                                  specific_genes=gene_ids)
    assert bounds_snapshot(ensemble.base_model) == before
