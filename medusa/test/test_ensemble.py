from cobra.core.model import Model
from cobra.io import load_model
from medusa.core.ensemble import Ensemble
from medusa.core.feature import Feature

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
    
    



def test_from_reaction_states_metabolites():
    textbook = load_model("textbook")
    biomass_id = "Biomass_Ecoli_core"
    biomass_rxn = textbook.reactions.get_by_id(biomass_id)

    # Pick one metabolite from biomass to vary across two members. Members
    # only specify the metabolite(s) they want to change; unspecified ones
    # remain at their baseline coefficient via add_metabolites(combine=False).
    target_met = next(iter(biomass_rxn.metabolites))
    states = {
        "alt_a": {target_met.id: -0.5},
        "alt_b": {target_met.id: -2.0},
    }

    ensemble = Ensemble.from_reaction_states(
        textbook,
        biomass_id,
        states,
        identifier="bof_test",
    )

    assert len(ensemble.features) == 1
    feature = ensemble.features[0]
    assert feature.id == f"{biomass_id}_metabolites"
    assert feature.component_attribute == "metabolites"
    assert feature.ensemble is ensemble
    assert feature.base_component in ensemble.base_model.reactions

    assert {m.id for m in ensemble.members} == {"alt_a", "alt_b"}
    for member in ensemble.members:
        assert member.ensemble is ensemble
        assert len(member.states) == 1

    # set_state must mutate the metabolite coefficient on the base_model
    # reaction (exercises the metabolites branch in Ensemble.set_state).
    base_biomass = ensemble.base_model.reactions.get_by_id(biomass_id)
    ensemble.set_state("alt_a")
    assert base_biomass.metabolites[target_met] == -0.5
    ensemble.set_state("alt_b")
    assert base_biomass.metabolites[target_met] == -2.0

    extracted = ensemble.extract_member("alt_a")
    assert isinstance(extracted, Model)


def test_from_reaction_states_rejects_unknown_metabolites():
    # By default, referencing a metabolite that is not currently in the
    # target reaction is a construction-time error.
    textbook = load_model("textbook")
    biomass_id = "Biomass_Ecoli_core"
    biomass = textbook.reactions.get_by_id(biomass_id)
    biomass_met_ids = {m.id for m in biomass.metabolites}

    # Find a metabolite in the model that is NOT in biomass.
    foreign_met = next(
        m for m in textbook.metabolites if m.id not in biomass_met_ids
    )

    with pytest.raises(ValueError):
        Ensemble.from_reaction_states(
            textbook,
            biomass_id,
            {"alt": {foreign_met.id: -0.1}},
        )


def test_from_reaction_states_allow_new_metabolites():
    # With allow_new_metabolites=True, members may introduce metabolites
    # that are not currently in the baseline reaction (e.g. alternative
    # energy carriers like swapping ATP for an analogue).
    textbook = load_model("textbook")
    biomass_id = "Biomass_Ecoli_core"
    biomass = textbook.reactions.get_by_id(biomass_id)
    biomass_met_ids = {m.id for m in biomass.metabolites}

    foreign_met = next(
        m for m in textbook.metabolites if m.id not in biomass_met_ids
    )

    ensemble = Ensemble.from_reaction_states(
        textbook,
        biomass_id,
        {"alt": {foreign_met: -0.1}},
        allow_new_metabolites=True,
    )

    ensemble.set_state("alt")
    base_biomass = ensemble.base_model.reactions.get_by_id(biomass_id)
    base_foreign = ensemble.base_model.metabolites.get_by_id(foreign_met.id)
    assert base_biomass.metabolites[base_foreign] == -0.1


def test_from_reaction_states_validation():
    textbook = load_model("textbook")

    with pytest.raises(ValueError):
        Ensemble.from_reaction_states(
            textbook, "not_a_real_rxn", {"a": {}})

    with pytest.raises(ValueError):
        Ensemble.from_reaction_states(
            textbook, "Biomass_Ecoli_core", {})


def test_features_kwarg_validation():
    textbook = load_model("textbook")
    biomass_id = "Biomass_Ecoli_core"
    biomass = textbook.reactions.get_by_id(biomass_id)
    feature = Feature(
        identifier=f"{biomass_id}_metabolites",
        base_component=biomass,
        component_attribute="metabolites",
        states={"a": {}, "b": {}},
    )

    # Empty features list — caller opted in but provided nothing.
    with pytest.raises(ValueError):
        Ensemble(list_of_models=[textbook], features=[])

    # features supplied alongside the wrong number of models.
    with pytest.raises(AttributeError):
        Ensemble(list_of_models=[textbook, textbook], features=[feature])
    with pytest.raises(AttributeError):
        Ensemble(list_of_models=[], features=[feature])


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
    reaction_ids = {rxn.id for rxn in test_ensemble.base_model.reactions}
    for feature in test_ensemble.features:
        assert feature.base_component.id in reaction_ids
        assert feature.component_attribute in REACTION_ATTRIBUTES
        assert len(set(feature.states.values())) > 1
