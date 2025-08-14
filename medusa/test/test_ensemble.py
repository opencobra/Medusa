from cobra.io import load_model
from medusa.core.ensemble import Ensemble

from pickle import load

import pytest

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

def construct_textbook_ensemble():
    # create two identical models and make an ensemble
    model1 = load_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = load_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    textbook_ensemble = Ensemble(list_of_models=[model1,model2],
                                        identifier='textbook_ensemble')
    return textbook_ensemble

def construct_mixed_ensemble():
    # create 4 models, which have reactions removed and a bound difference.
    model1 = load_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = load_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = load_model("textbook")
    model3.remove_reactions(model3.reactions[5:7])
    model3.id = 'third_textbook'
    model4 = model3.copy()
    model4.id = 'dual_features'
    model4.reactions[1].lower_bound = 0
    model_list = [model1,model2,model3,model4]
    mixed_ensemble = Ensemble(list_of_models=model_list,identifier='textbook_ensemble')
    return(mixed_ensemble)

def test_ensemble_creation():
    # test whether ensemble components are properly generated in test ensembles
    test_ensemble = construct_textbook_ensemble()
    # The base model should have the same number of reactions and metabolites
    # as the original model, since we only remove/modify reactions.
    textbook = load_model("textbook")
    assert len(test_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_ensemble.base_model.metabolites) == len(textbook.metabolites)

    # the ensemble should have 8 features and 2 members
    assert len(test_ensemble.features) == 8
    assert len(test_ensemble.members) == 2

    # each member in the ensemble should have 8 features and values in their states
    # each member should have a reference to the correct ensemble object
    for member in test_ensemble.members:
        assert len(member.states) == 8
        assert member.ensemble == test_ensemble

    # each feature should reference a reaction contained in the ensemble
    # each feature should have a component_attribute in the list of allowable
    # attributes
    # each feature should have at least two unique state values across all models
    for feature in test_ensemble.features:
        assert feature.base_component in test_ensemble.base_model.reactions
        assert feature.component_attribute in REACTION_ATTRIBUTES
        assert len(set(feature.states.values())) > 1

def test_mixed_ensemble_creation():
    # Same as basic test, but with a member that had a bound change rather than
    # reaction removal
    test_ensemble = construct_mixed_ensemble()
    # The base model should have the same number of reactions and metabolites
    # as the original model, since we only remove/modify reactions.
    textbook = load_model("textbook")
    assert len(test_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_ensemble.base_model.metabolites) == len(textbook.metabolites)

    # the ensemble should have 10 features and 4 members
    assert len(test_ensemble.features) == 10
    assert len(test_ensemble.members) == 4

    # each member in the ensemble should have 8 features and values in their states
    # each member should have a reference to the correct ensemble object
    for member in test_ensemble.members:
        assert len(member.states) == 10
        assert member.ensemble == test_ensemble

    # each feature should reference a reaction contained in the ensemble
    # each feature should have a component_attribute in the list of allowable
    # attributes
    # each feature should have at least two unique state values across all models
    for feature in test_ensemble.features:
        assert feature.base_component in test_ensemble.base_model.reactions
        assert feature.component_attribute in REACTION_ATTRIBUTES
        assert len(set(feature.states.values())) > 1

def test_extract_member():
    test_ensemble = construct_textbook_ensemble()
    original_model = load_model("textbook")

    extracted_member = test_ensemble.extract_member(test_ensemble.members[0])

    # check that the original reactions were removed
    unique_feature_comps = set([
            feat.base_component for feat in test_ensemble.features])

    for feature in test_ensemble.members[0].states.keys():
        if test_ensemble.members[0].states[feature] == 0:
            assert feature.base_component.id not in [
                rxn.id for rxn in extracted_member.reactions]

def test_update_member_id():
    # updating the id on a member should update the index in members,
    # the id in feature.states. Also, attempting to set the member id
    # to an existing member should raise an error.

    test_ensemble = construct_textbook_ensemble()
    member1 = test_ensemble.members[0]
    member2 = test_ensemble.members[1]

    new_id = 'first_with_mod_id'
    member1.id = new_id
    assert member1.id == new_id
    assert test_ensemble.members[0].id == new_id
    assert test_ensemble.members.get_by_id(new_id)
    assert new_id in test_ensemble.features[0].states.keys()

    duplicate_id = member2.id
    with pytest.raises(ValueError):
        member1.id = duplicate_id
    
    



def test_pickle():
    test_ensemble = construct_mixed_ensemble()

    # pickle and unpickle the ensemble, then rerun test_mixed_ensemble_creation
    save_loc = 'test_pickle.pickle'
    test_ensemble.to_pickle(save_loc)

    with open(save_loc,'rb') as infile:
        unpickled = load(infile)

    test_ensemble = unpickled
    textbook = load_model("textbook")
    assert len(test_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_ensemble.base_model.metabolites) == len(textbook.metabolites)

    # the ensemble should have 10 features and 4 members
    assert len(test_ensemble.features) == 10
    assert len(test_ensemble.members) == 4

    # each member in the ensemble should have 8 features and values in their states
    # each member should have a reference to the correct ensemble object
    for member in test_ensemble.members:
        assert len(member.states) == 10
        assert member.ensemble == test_ensemble

    # each feature should reference a reaction contained in the ensemble
    # each feature should have a component_attribute in the list of allowable
    # attributes
    # each feature should have at least two unique state values across all models
    for feature in test_ensemble.features:
        assert feature.base_component in test_ensemble.base_model.reactions
        assert feature.component_attribute in REACTION_ATTRIBUTES
        assert len(set(feature.states.values())) > 1
