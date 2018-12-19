
from cobra.test import create_test_model
from medusa.core.ensemble import Ensemble
from medusa.flux_analysis.flux_balance import optimize_ensemble


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

def test_fba_return_dims():
        ensemble = construct_textbook_ensemble()
        fba_fluxes = optimize_ensemble(ensemble)
        rows = fba_fluxes.shape[0]
        columns = fba_fluxes.shape[1]
        assert rows == len(ensemble.members)
        assert columns == len(ensemble.base_model.reactions)

def test_fba_single_return():
        # test return of flux values for a single reaction
        ensemble = construct_textbook_ensemble()
        fba_fluxes = optimize_ensemble(ensemble,
                                return_flux=ensemble.base_model.reactions[2].id)
        rows = fba_fluxes.shape[0]
        columns = fba_fluxes.shape[1]
        assert rows == len(ensemble.members)
        assert columns == 1

def test_fba_random_sample():
        # test return of FBA on randomly sampled ensemble members
        ensemble = construct_mixed_ensemble()
        nmodels = 2
        fba_fluxes = optimize_ensemble(ensemble,
                                num_models=nmodels)
        rows = fba_fluxes.shape[0]
        columns = fba_fluxes.shape[1]
        assert rows == nmodels
        assert columns == len(ensemble.base_model.reactions)

def test_fba_specific_models():
        # test return of FBA on specific ensemble members
        ensemble = construct_mixed_ensemble()
        model1 = ensemble.members[0]
        model2 = ensemble.members[1]
        model_list = [model1,model2]
        fba_fluxes = optimize_ensemble(ensemble,
                                specific_models=model_list)
        rows = fba_fluxes.shape[0]
        columns = fba_fluxes.shape[1]
        rownames = fba_fluxes.index
        assert rows == len(model_list)
        assert columns == len(ensemble.base_model.reactions)

        assert rownames.contains(model1)
        assert rownames.contains(model2)
