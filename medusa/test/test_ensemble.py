from cobra.test import create_test_model
from medusa.core.ensemble import Ensemble
from medusa.core.ensemble import update_features_states


def test_concurrent_join():
    # test creation of an ensemble from a list of models with all models loaded
    # in memory
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[2:5])
    model3.id = 'third_textbook'
    textbook_ensemble = Ensemble(model_list=[model1,model2,model3],
                                        base_id='textbook_ensemble')
    # are there three members in the ensemble?
    num_models = 3
    assert len(textbook_ensemble.features.index) == num_models
    # are there 5 reactions in the reaction diff for each model?
    for model in textbook_ensemble.features.index.tolist():
        assert len(textbook_ensemble.states.loc[model]) == 5

def test_iterative_join():
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[2:5])
    model3.id = 'third_textbook'
    textbook_ensemble = Ensemble(model_list=[model1,model2,model3],\
                                base_id='textbook_ensemble',\
                                join_method='iterative')
    # are there three members in the ensemble?
    num_models = 3
    assert len(textbook_ensemble.features.index) == num_models
    # are there 5 reactions in the reaction diff for each model?
    for model in textbook_ensemble.features.index.tolist():
        assert len(textbook_ensemble.states.loc[model]) == 5

def test_iterative_join_size():
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[2:5])
    model3.id = 'third_textbook'
    textbook_ensemble = Ensemble(model_list=[model1,model2,model3],\
                                base_id='textbook_ensemble',\
                                join_method='iterative',join_size=2)

    # are there three members in the ensemble?
    num_models = 3
    assert len(textbook_ensemble.features.index) == num_models
    # are there 5 reactions in the reaction diff for each model?
    for model in textbook_ensemble.features.index.tolist():
        assert len(textbook_ensemble.states.loc[model]) == 5


def test_empty_join():
    # test creation of an ensemble from a list of models when starting with an
    # empty ensemble object
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[2:5])
    model3.id = 'third_textbook'
    # do it the normal way first
    textbook_ensemble = Ensemble(model_list=[model1,model2,model3],
                                        base_id='textbook_ensemble')
    from_empty = Ensemble(base_id='ensemble_from_empty')
    update_features_states(from_empty,model_list=[model1,model2,model3])

    # are there three members in the ensemble?
    assert len(textbook_ensemble.reaction_diffs.keys()) == \
        len(from_empty.reaction_diffs.keys())
    # are there 5 reactions in the reaction diff for each model?
    for model in textbook_ensemble.reaction_diffs.keys():
        assert len(textbook_ensemble.reaction_diffs[model]) == \
            len(from_empty.reaction_diffs[model])
