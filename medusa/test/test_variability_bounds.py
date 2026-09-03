from cobra.io import load_model
from cobra.flux_analysis.variability import flux_variability_analysis

from medusa.core.ensemble import Ensemble
from medusa.flux_analysis.variability import ensemble_fva

import pytest


def construct_textbook_ensemble():
    model1 = load_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = load_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = load_model("textbook")
    model3.remove_reactions(model3.reactions[2:5])
    model3.id = 'third_textbook'
    return Ensemble(list_of_models=[model1, model2, model3],
                    identifier='textbook_ensemble')


def test_ensemble_fva_maximum_exceeds_minimum():
    """Regression: the 'maximum_' and 'minimum_' rows used to be swapped.

    ensemble_fva relabelled cobrapy's FVA output positionally with
    ['maximum_<id>', 'minimum_<id>'], but cobrapy returns its columns in the
    order ['minimum', 'maximum'], so the row labelled 'maximum_<id>' actually
    held that member's minima. Anything that read the two rows by label -- a
    flux range, a plot, an interval -- got them inverted.
    """
    ensemble = construct_textbook_ensemble()
    reaction_list = [rxn.id for rxn in ensemble.base_model.reactions
                     if rxn.id.startswith('EX_')]
    # A fraction below 1.0 guarantees the ranges are actually wide, so the
    # assertion below has something to bite on.
    fva_fluxes = ensemble_fva(ensemble, reaction_list=reaction_list,
                              fraction_of_optimum=0.9)

    flux_columns = [column for column in fva_fluxes.columns
                    if column != 'model_source']
    assert set(flux_columns) == set(reaction_list)

    strictly_wider = 0
    for member in ensemble.members:
        maximum_row = 'maximum_' + member.id
        minimum_row = 'minimum_' + member.id
        assert maximum_row in fva_fluxes.index
        assert minimum_row in fva_fluxes.index
        for reaction_id in flux_columns:
            maximum = fva_fluxes.loc[maximum_row, reaction_id]
            minimum = fva_fluxes.loc[minimum_row, reaction_id]
            assert maximum >= minimum - 1e-6, (
                "maximum_%s < minimum_%s for %s (%r < %r)"
                % (member.id, member.id, reaction_id, maximum, minimum))
            if maximum > minimum + 1e-6:
                strictly_wider += 1

    # Guard the guard: if every range were degenerate the assertion above
    # would pass trivially and would not have caught the original bug.
    assert strictly_wider > 0


def test_ensemble_fva_labels_match_cobrapy():
    """The labelled rows must hold the values cobrapy calls min and max."""
    ensemble = construct_textbook_ensemble()
    reaction_list = [rxn.id for rxn in ensemble.base_model.reactions
                     if rxn.id.startswith('EX_')]

    member_id = ensemble.members[0].id
    fva_fluxes = ensemble_fva(ensemble, reaction_list=reaction_list,
                              specific_models=[member_id],
                              fraction_of_optimum=0.9)

    with ensemble.base_model:
        ensemble.set_state(member_id)
        reference = flux_variability_analysis(
            ensemble.base_model, reaction_list=reaction_list,
            fraction_of_optimum=0.9)

    for reaction_id in reaction_list:
        assert fva_fluxes.loc['maximum_' + member_id, reaction_id] == \
            pytest.approx(reference.loc[reaction_id, 'maximum'], abs=1e-6)
        assert fva_fluxes.loc['minimum_' + member_id, reaction_id] == \
            pytest.approx(reference.loc[reaction_id, 'minimum'], abs=1e-6)
