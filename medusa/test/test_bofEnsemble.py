from cobra.core.model import Model
from cobra.io import load_model
from medusa.core.ensemble import Ensemble

# TODO update later if bofEnsemble gets moved
from medusa.bofEnsemble import bofEnsemble
from medusa.bofEnsemble import _getBofDf

from pickle import load

import pytest

def test_bofEnsemble():
    textbook = load_model("textbook")

    # based on the textbook model, generate an ensemble in which the members
    # only differ w.r.t. their BOF metabolite coefficients
    BofDf = _getBofDf(textbook, BofId = 'Biomass_Ecoli_core', n_models = 10)
    test_ensemble = bofEnsemble(textbook, BofDf, 'Biomass_Ecoli_core')

    assert len(test_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_ensemble.base_model.metabolites) == len(textbook.metabolites)

    # the ensemble should have 1 feature and 11 members
    assert len(test_ensemble.features) == 1
    assert len(test_ensemble.members) == 11

    # each member in the ensemble should have 1 feature and values in their states
    # each member should have a reference to the correct ensemble object
    for member in test_ensemble.members:
        assert len(member.states) == 1
        assert member.ensemble == test_ensemble

    # the BOF should reference a reaction contained in the ensemble
    assert test_ensemble.features[0].base_component in test_ensemble.base_model.reactions

    # it should be possible to extract an individual member
    model_from_id = test_ensemble.extract_member('model_0')
    assert type(model_from_id) is Model