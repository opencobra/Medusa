
from cobra.test import create_test_model
from medusa.core.ensemble import Ensemble
from medusa.flux_analysis import flux_balance


def construct_textbook_ensemble():
    # create two identical models and make an ensemble
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    textbook_ensemble = Ensemble(model_list=[model1,model2],
                                        base_id='textbook_ensemble')
    return textbook_ensemble
    
def test_fba_return_dims():
        ensemble = construct_textbook_ensemble()
        fba_fluxes = flux_balance.optimize_ensemble(ensemble)
        rows = fba_fluxes.shape[0]
        columns = fba_fluxes.shape[1]
        assert rows == len(ensemble.reaction_diffs.keys())
        assert columns == len(ensemble.base_model.reactions)
