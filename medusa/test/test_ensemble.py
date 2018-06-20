from cobra.test import create_test_model
from medusa.core.ensemble import Ensemble

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

def construct_textbook_ensemble():
    # create two identical models and make an ensemble
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    textbook_ensemble = Ensemble(list_of_models=[model1,model2],
                                        identifier='textbook_ensemble')
    return textbook_ensemble

def construct_mixed_ensemble():
    # create 4 models, which have reactions removed and a bound difference.
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = create_test_model("textbook")
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
    textbook = create_test_model("textbook")
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
    textbook = create_test_model("textbook")
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




# def test_concurrent_join():
#     # test creation of an ensemble from a list of models with all models loaded
#     # in memory
#     model1 = create_test_model("textbook")
#     model1.remove_reactions(model1.reactions[1:3])
#     model1.id = 'first_textbook'
#     model2 = create_test_model("textbook")
#     model2.remove_reactions(model2.reactions[4:6])
#     model2.id = 'second_textbook'
#     model3 = create_test_model("textbook")
#     model3.remove_reactions(model3.reactions[2:5])
#     model3.id = 'third_textbook'
#     textbook_ensemble = Ensemble(model_list=[model1,model2,model3],
#                                         base_id='textbook_ensemble')
#     # are there three members in the ensemble?
#     num_models = 3
#     assert len(textbook_ensemble.states.index) == num_models
#     # are there 5 reactions in the states and features?
#     assert len(textbook_ensemble.states.columns) == 5
#     assert len(textbook_ensemble.features.index) == 5

# def test_iterative_join():
#     model1 = create_test_model("textbook")
#     model1.remove_reactions(model1.reactions[1:3])
#     model1.id = 'first_textbook'
#     model2 = create_test_model("textbook")
#     model2.remove_reactions(model2.reactions[4:6])
#     model2.id = 'second_textbook'
#     model3 = create_test_model("textbook")
#     model3.remove_reactions(model3.reactions[2:5])
#     model3.id = 'third_textbook'
#     textbook_ensemble = Ensemble(model_list=[model1,model2,model3],\
#                                 base_id='textbook_ensemble',\
#                                 join_method='iterative')
#     # are there three members in the ensemble?
#     num_models = 3
#     assert len(textbook_ensemble.states.index) == num_models
#     # are there 5 reactions in the states and features?
#     assert len(textbook_ensemble.states.columns) == 5
#     assert len(textbook_ensemble.features.index) == 5

# def test_iterative_join_size():
#     model1 = create_test_model("textbook")
#     model1.remove_reactions(model1.reactions[1:3])
#     model1.id = 'first_textbook'
#     model2 = create_test_model("textbook")
#     model2.remove_reactions(model2.reactions[4:6])
#     model2.id = 'second_textbook'
#     model3 = create_test_model("textbook")
#     model3.remove_reactions(model3.reactions[2:5])
#     model3.id = 'third_textbook'
#     textbook_ensemble = Ensemble(model_list=[model1,model2,model3],\
#                                 base_id='textbook_ensemble',\
#                                 join_method='iterative',join_size=2)
#
#     # are there three members in the ensemble?
#     num_models = 3
#     assert len(textbook_ensemble.states.index) == num_models
#     # are there 5 reactions in the states and features?
#     assert len(textbook_ensemble.states.columns) == 5
#     assert len(textbook_ensemble.features.index) == 5


# def test_empty_join():
#     # test creation of an ensemble from a list of models when starting with an
#     # empty ensemble object
#     model1 = create_test_model("textbook")
#     model1.remove_reactions(model1.reactions[1:3])
#     model1.id = 'first_textbook'
#     model2 = create_test_model("textbook")
#     model2.remove_reactions(model2.reactions[4:6])
#     model2.id = 'second_textbook'
#     model3 = create_test_model("textbook")
#     model3.remove_reactions(model3.reactions[2:5])
#     model3.id = 'third_textbook'
#     # do it the normal way first
#     textbook_ensemble = Ensemble(model_list=[model1,model2,model3],
#                                         base_id='textbook_ensemble')
#     from_empty = Ensemble(base_id='ensemble_from_empty')
#     update_features_states(from_empty,model_list=[model1,model2,model3])
#
#     # are there the same number of models in the ensemble?
#     assert len(textbook_ensemble.states.index) == \
#         len(from_empty.states.index)
#     # are there the same number of reactions in the states and features?
#     assert len(textbook_ensemble.states.columns) == len(from_empty.states.columns)
#     assert len(textbook_ensemble.features.index) == len(from_empty.features.index)

#def test_from_json():
    # test functions to load an ensemble from JSON files.
