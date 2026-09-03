"""Tests for medusa.flux_analysis.deletion.

The gene-deletion entry point previously defaulted ``specific_genes`` to an
empty list and forwarded it to cobrapy's ``gene_list``. cobrapy treats only
``None`` as "every gene"; an empty list means "no genes". The default call
therefore deleted nothing and returned an empty result for every member, with
no error and no warning, while the docstring promised that all genes would be
selected. ``specific_genes`` is now required, and every way of expressing an
empty selection raises rather than silently doing nothing.
"""

import pytest
from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.flux_analysis.deletion import (ensemble_single_gene_deletion,
                                           ensemble_single_reaction_deletion)


def _textbook(name, glc=-10.0):
    model = load_model("textbook")
    model.id = name
    model.reactions.EX_glc__D_e.lower_bound = glc
    return model


@pytest.fixture
def ensemble():
    return Ensemble(
        list_of_models=[_textbook("plenty", glc=-10.0),
                        _textbook("scarce", glc=-5.0)],
        identifier="deletion_ensemble")


# --------------------------------------------------------------------------
# Empty and missing selections must raise, not silently return nothing.
# --------------------------------------------------------------------------

def test_gene_deletion_without_specific_genes_raises(ensemble):
    """The regression test: the no-argument call used to return empty results."""
    with pytest.raises(ValueError, match="specific_genes is required"):
        ensemble_single_gene_deletion(ensemble)


def test_gene_deletion_explicit_none_raises(ensemble):
    with pytest.raises(ValueError, match="specific_genes is required"):
        ensemble_single_gene_deletion(ensemble, specific_genes=None)


@pytest.mark.parametrize("empty", [[], (), set(), ""])
def test_gene_deletion_empty_selection_raises(ensemble, empty):
    with pytest.raises(ValueError, match="specific_genes is empty"):
        ensemble_single_gene_deletion(ensemble, specific_genes=empty)


def test_gene_deletion_blank_entry_in_list_raises(ensemble):
    gene_id = ensemble.base_model.genes[0].id
    with pytest.raises(ValueError, match="empty or unidentifiable"):
        ensemble_single_gene_deletion(ensemble, specific_genes=[gene_id, ""])


def test_gene_deletion_non_iterable_raises(ensemble):
    with pytest.raises(TypeError, match="must be a gene id"):
        ensemble_single_gene_deletion(ensemble, specific_genes=7)


def test_gene_deletion_unknown_gene_raises(ensemble):
    with pytest.raises(KeyError, match="not in the ensemble's base model"):
        ensemble_single_gene_deletion(ensemble,
                                      specific_genes=["not_a_real_gene"])


# --------------------------------------------------------------------------
# Valid selections do real work.
# --------------------------------------------------------------------------

def test_gene_deletion_returns_one_row_per_requested_gene(ensemble):
    gene_ids = [gene.id for gene in ensemble.base_model.genes[:3]]
    results = ensemble_single_gene_deletion(ensemble,
                                            specific_genes=gene_ids)
    assert set(results) == {member.id for member in ensemble.members}
    for member_id, frame in results.items():
        assert frame.shape[0] == len(gene_ids), (
            "member %s got %i rows for %i requested genes"
            % (member_id, frame.shape[0], len(gene_ids)))


def test_gene_deletion_accepts_a_single_gene_id(ensemble):
    gene_id = ensemble.base_model.genes[0].id
    results = ensemble_single_gene_deletion(ensemble, specific_genes=gene_id)
    for frame in results.values():
        assert frame.shape[0] == 1


def test_gene_deletion_accepts_gene_objects(ensemble):
    genes = list(ensemble.base_model.genes[:2])
    results = ensemble_single_gene_deletion(ensemble, specific_genes=genes)
    for frame in results.values():
        assert frame.shape[0] == len(genes)


def test_gene_deletion_all_genes_is_expressible(ensemble):
    """The documented replacement for the old implicit 'all genes' default."""
    gene_ids = [gene.id for gene in ensemble.base_model.genes]
    results = ensemble_single_gene_deletion(
        ensemble, specific_models=["plenty"], specific_genes=gene_ids)
    assert results["plenty"].shape[0] == len(gene_ids)
    assert len(gene_ids) > 1


# --------------------------------------------------------------------------
# Member selection is consistent across entry points.
# --------------------------------------------------------------------------

def test_deletion_accepts_member_ids_as_strings(ensemble):
    results = ensemble_single_reaction_deletion(ensemble,
                                                specific_models=["plenty"])
    assert list(results) == ["plenty"]


def test_deletion_accepts_member_objects(ensemble):
    member = ensemble.members[0]
    results = ensemble_single_reaction_deletion(ensemble,
                                                specific_models=[member])
    assert list(results) == [member.id]


def test_deletion_empty_specific_models_raises(ensemble):
    with pytest.raises(ValueError, match="specific_models is empty"):
        ensemble_single_reaction_deletion(ensemble, specific_models=[])


def test_deletion_unknown_member_raises(ensemble):
    with pytest.raises(KeyError, match="not members of this ensemble"):
        ensemble_single_reaction_deletion(ensemble,
                                          specific_models=["nope"])


def test_deletion_rejects_both_selectors(ensemble):
    with pytest.raises(ValueError, match="cannot be passed concurrently"):
        ensemble_single_reaction_deletion(ensemble, num_models=1,
                                          specific_models=["plenty"])
