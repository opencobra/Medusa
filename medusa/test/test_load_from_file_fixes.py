"""Tests for defects in medusa.reconstruct.load_from_file.

- The batching arithmetic in ``batch_load_from_files`` raised IndexError for
  ordinary combinations of file count and batch size, and could leave a batch
  holding a single model, which produces a featureless ensemble.
- ``add_ensembles`` excluded reactions from feature creation by reaction id
  rather than feature id, so a reaction that already had a lower_bound feature
  could never receive the upper_bound feature it also needed. The difference
  in upper bounds between the two ensembles was dropped silently.
- ``add_ensembles`` rewired both of its inputs, despite appearing to deep-copy.
"""

import pytest
from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.reconstruct.load_from_file import _batch_indices, add_ensembles


def _textbook(name, pgi_lower=-1000.0, pgi_upper=1000.0):
    model = load_model("textbook")
    model.id = name
    model.reactions.PGI.bounds = (pgi_lower, pgi_upper)
    return model


# --------------------------------------------------------------------------
# Batching.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("total,batchsize", [
    (2, 5), (5, 5), (6, 5), (7, 5), (5, 2), (6, 2), (7, 3), (10, 4), (3, 2),
])
def test_batches_cover_everything_exactly_once(total, batchsize):
    batches = _batch_indices(total, batchsize)
    flattened = [index for batch in batches for index in batch]
    assert sorted(flattened) == list(range(total))
    assert len(flattened) == len(set(flattened))


@pytest.mark.parametrize("total,batchsize", [
    (2, 5), (5, 5), (6, 5), (7, 5), (5, 2), (6, 2), (7, 3), (10, 4), (3, 2),
])
def test_no_batch_holds_a_single_model(total, batchsize):
    """A one-model batch builds an ensemble with no features, which cannot
    then be merged."""
    for batch in _batch_indices(total, batchsize):
        assert len(batch) >= 2, _batch_indices(total, batchsize)


def test_batchsize_of_one_is_rejected():
    with pytest.raises(ValueError, match="batchsize must be at least 2"):
        _batch_indices(5, 1)


def test_single_file_is_rejected():
    with pytest.raises(ValueError, match="At least 2 models"):
        _batch_indices(1, 5)


# --------------------------------------------------------------------------
# add_ensembles must not lose an attribute difference.
# --------------------------------------------------------------------------

def _two_ensembles():
    """Both ensembles vary PGI's lower bound; they disagree on its upper."""
    e1 = Ensemble(
        list_of_models=[_textbook("a1", pgi_lower=-1000.0, pgi_upper=1000.0),
                        _textbook("a2", pgi_lower=0.0, pgi_upper=1000.0)],
        identifier="e1")
    e2 = Ensemble(
        list_of_models=[_textbook("b1", pgi_lower=-1000.0, pgi_upper=500.0),
                        _textbook("b2", pgi_lower=0.0, pgi_upper=500.0)],
        identifier="e2")
    return e1, e2


def test_upper_bound_difference_is_not_dropped():
    e1, e2 = _two_ensembles()
    merged = add_ensembles(e1, e2)

    feature_ids = {feature.id for feature in merged.features}
    assert "PGI_lower_bound" in feature_ids
    assert "PGI_upper_bound" in feature_ids, (
        "the two ensembles disagree on PGI's upper bound, but the reaction "
        "already had a lower_bound feature so the difference was dropped")

    seen = {}
    for member in merged.members:
        merged.set_state(member.id)
        seen[member.id] = merged.base_model.reactions.PGI.upper_bound
    assert seen["a1"] == pytest.approx(1000.0)
    assert seen["b1"] == pytest.approx(500.0), (
        "member b1 should keep its own upper bound of 500; got %r"
        % seen["b1"])


# --------------------------------------------------------------------------
# add_ensembles must leave its inputs alone.
# --------------------------------------------------------------------------

def test_inputs_are_not_rewired():
    e1, e2 = _two_ensembles()
    e1_members = [member.id for member in e1.members]
    e2_members = [member.id for member in e2.members]

    merged = add_ensembles(e1, e2)

    for member in e1.members:
        assert member.ensemble is e1, \
            "a member of e1 was rebound to the merged ensemble"
    for member in e2.members:
        assert member.ensemble is e2, \
            "a member of e2 was rebound to the merged ensemble"
    assert [member.id for member in e1.members] == e1_members
    assert [member.id for member in e2.members] == e2_members
    assert merged is not e1 and merged is not e2


def test_input_features_still_belong_to_their_own_ensemble():
    e1, e2 = _two_ensembles()
    add_ensembles(e1, e2)
    for feature in e1.features:
        assert feature.ensemble is e1
    for feature in e2.features:
        assert feature.ensemble is e2


def test_inputs_remain_usable_after_merging():
    e1, e2 = _two_ensembles()
    before = {}
    for member in e1.members:
        e1.set_state(member.id)
        before[member.id] = e1.base_model.reactions.PGI.bounds

    add_ensembles(e1, e2)

    for member in e1.members:
        e1.set_state(member.id)
        assert e1.base_model.reactions.PGI.bounds == before[member.id], \
            "e1 no longer reproduces its own states after being merged"


def test_overlapping_member_ids_are_rejected():
    e1 = Ensemble(
        list_of_models=[_textbook("same1"), _textbook("same2", pgi_lower=0.0)],
        identifier="e1")
    e2 = Ensemble(
        list_of_models=[_textbook("same1"), _textbook("same2", pgi_lower=0.0)],
        identifier="e2")
    with pytest.raises(ValueError, match="share member id"):
        add_ensembles(e1, e2)


def test_merged_feature_order_is_deterministic():
    e1, e2 = _two_ensembles()
    merged = add_ensembles(e1, e2)
    ids = [feature.id for feature in merged.features]
    assert len(ids) == len(set(ids)), "duplicate features in the merge"
