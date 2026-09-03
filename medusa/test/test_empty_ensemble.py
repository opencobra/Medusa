"""Ensemble must expose features/members even when it has none of them.

``Ensemble.__init__`` assigned `features` and `members` only on the path that
takes two or more models. An ensemble built from zero models (the documented
way to create an empty one) or from a single model left both attributes
unassigned, so ordinary attribute access raised AttributeError rather than
returning an empty DictList.

Both shapes are real: `Ensemble(identifier=...)` with no models is how the
class documents creating an empty ensemble, and the single-model form is what
`boundsEnsemble` and `medusa.reconstruct.expand` build before populating the
features themselves.
"""

from cobra.io import load_model

from medusa.core.ensemble import Ensemble


def test_ensemble_with_no_models_has_empty_features_and_members():
    empty = Ensemble(identifier='empty')
    assert len(empty.features) == 0
    assert len(empty.members) == 0


def test_ensemble_with_one_model_has_empty_features_and_members():
    single = Ensemble(list_of_models=[load_model("textbook")],
                      identifier='single')
    assert len(single.features) == 0
    assert len(single.members) == 0


def test_empty_ensemble_features_and_members_are_independent():
    """Each ensemble needs its own DictLists, not a shared default."""
    one = Ensemble(identifier='one')
    two = Ensemble(identifier='two')
    assert one.features is not two.features
    assert one.members is not two.members
