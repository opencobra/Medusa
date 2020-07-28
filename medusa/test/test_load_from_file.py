import cobra
from cobra.test import create_test_model
from medusa.core.ensemble import Ensemble
from medusa.reconstruct.load_from_file import batch_load_from_files

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

def construct_mixed_ensemble_2():
    # create 4 models, which have reactions removed and a bound difference.
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    cobra.io.save_json_model(model1, "model1.json")
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    cobra.io.save_json_model(model2, "model2.json")
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[5:7])
    model3.id = 'third_textbook'
    model3.reactions.get_by_id('SUCCt2_2').lower_bound = -1000
    cobra.io.save_json_model(model3, "model3.json")
    model4 = model3.copy() # missing rxns [5:7] and has lower_bound=-1000 for SUCCt2_2
    model4.id = 'dual_features'
    model4.reactions[1].lower_bound = 0
    cobra.io.save_json_model(model4, "model4.json")
    model_list = [model1,model2,model3,model4]
    mixed_ensemble = Ensemble(list_of_models=model_list,identifier='textbook_ensemble')
    return(mixed_ensemble)

def construct_mixed_batch_ensemble():
    
    # UPDATE PATHS TO SAVE AND LOAD MODELS
    
    # create 4 models, which have reactions removed and a bound difference.
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    cobra.io.save_json_model(model1, "model1.json")
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    cobra.io.save_json_model(model2, "model2.json")
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[5:7])
    model3.id = 'third_textbook'
    model3.reactions.get_by_id('SUCCt2_2').lower_bound = -1000
    cobra.io.save_json_model(model3, "model3.json")
    model4 = model3.copy() # missing rxns [5:7] and has lower_bound=-1000 for SUCCt2_2
    model4.id = 'dual_features'
    model4.reactions[1].lower_bound = 0
    cobra.io.save_json_model(model4, "model4.json")
    
    model_file_names = ["model1.json","model2.json","model3.json","model4.json"]
    batch_mixed_ensemble = batch_load_from_files(model_file_names,identifier='textbook_ensemble',batchsize = 2)
    return(batch_mixed_ensemble)

def test_batch_load_vs_innate():
    # Same as basic test, but with a member that had a bound change rather than
    # reaction removal
    test_ensemble = construct_mixed_ensemble_2()
    test_batch_ensemble = construct_mixed_batch_ensemble()
    # The base model should have the same number of reactions and metabolites
    # as the original model, since we only remove/modify reactions.
    textbook = create_test_model("textbook")
    assert len(test_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_ensemble.base_model.metabolites) == len(textbook.metabolites)
    assert len(test_batch_ensemble.base_model.reactions) == len(textbook.reactions)
    assert len(test_batch_ensemble.base_model.metabolites) == len(textbook.metabolites)

    # the ensemble should have 10 features and 4 members
    assert len(test_ensemble.features) == 11
    assert len(test_ensemble.members) == 4
    assert len(test_batch_ensemble.features) == 11
    assert len(test_batch_ensemble.members) == 4

    # each member in the ensemble should have 11 features and values in their states
    # each member should have a reference to the correct ensemble object
    for member in test_ensemble.members:
        assert len(member.states) == 11
        assert member.ensemble == test_ensemble
    for member in test_batch_ensemble.members:
        assert len(member.states) == 11
        assert member.ensemble == test_batch_ensemble

    # each feature should reference a reaction contained in the ensemble
    # each feature should have a component_attribute in the list of allowable
    # attributes
    # each feature should have at least two unique state values across all models
    for feature in test_ensemble.features:
        assert feature.base_component in test_ensemble.base_model.reactions
        assert feature.component_attribute in REACTION_ATTRIBUTES
        assert len(set(feature.states.values())) > 1
    for feature in test_batch_ensemble.features:
        assert feature.base_component in test_batch_ensemble.base_model.reactions
        assert feature.component_attribute in REACTION_ATTRIBUTES
        assert len(set(feature.states.values())) > 1
        
def test_all_attributes_in_batch_load_model():
    test_ensemble = construct_mixed_ensemble_2()
    test_batch_ensemble = construct_mixed_batch_ensemble()
    
    # Test number of features and members are equal
    assert len(test_ensemble.features) == len(test_batch_ensemble.features)
    assert len(test_ensemble.members) == len(test_batch_ensemble.members)
    
    # Test that all the feature states are the same
    for feat_obj in test_batch_ensemble.features:
        comp_feat = test_ensemble.features.get_by_id(feat_obj.id)
        for key, value in feat_obj.states.items():
            assert comp_feat.states[key] == value
        
    # Test that all the member states are the same
    for mem_obj in test_batch_ensemble.members:
        comp_mem_obj = test_ensemble.members.get_by_id(mem_obj.id)
        for key, value in mem_obj.states.items():
            assert comp_mem_obj.states[test_ensemble.features.get_by_id(key.id)] == value
        
    # Test that the reactions and feature.base_components are the same objects
    for feat in test_batch_ensemble.features:
        assert feat.base_component == test_batch_ensemble.base_model.reactions.get_by_id(feat.base_component.id)