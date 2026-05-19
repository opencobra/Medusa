from cobra.core.model import Model
from cobra.io import load_model
from medusa.core.ensemble import Ensemble

# TODO update later if boundsEnsemble gets moved
from medusa.boundsEnsemble import boundsEnsemble
from medusa.boundsEnsemble import _setBounds

import pytest

def test_boundsEnsemble():

    textbook = load_model("textbook")

    # based on the textbook model, generate an ensemble in which the members
    # only differ w.r.t. their bounds for specified reactions
    boundsDict = _setBounds(textbook, 
                            [rxn.id for rxn in textbook.reactions[0:5]], 
                            method='random', 
                            reversibility=None, 
                            bound=1000, 
                            n_models=10)
    test_ensemble = boundsEnsemble(textbook, boundsDict)

    assert len(test_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_ensemble.base_model.metabolites) == len(textbook.metabolites)

    # the ensemble should have 1 feature and 11 members
    assert len(test_ensemble.features) == 10
    assert len(test_ensemble.members) == 11

    # each member in the ensemble should have 10 features and values in their states
    # each member should have a reference to the correct ensemble object
    for member in test_ensemble.members:
        assert len(member.states) == 10
        assert member.ensemble == test_ensemble

    # each feature should reference a reaction contained in the original model
    assert [feature.base_component in test_ensemble.base_model.reactions for feature in test_ensemble.features]

    # it should be possible to extract an individual member
    model_from_id = test_ensemble.extract_member('model_0')
    assert type(model_from_id) is Model