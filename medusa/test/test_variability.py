from cobra.test import create_test_model
from medusa.core.ensemble import Ensemble
from medusa.flux_analysis.variability import ensemble_fva


def construct_textbook_ensemble():
    # create two identical models and make an ensemble
    model1 = create_test_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = create_test_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    model3 = create_test_model("textbook")
    model3.remove_reactions(model3.reactions[2:5])
    model3.id = 'third_textbook'
    textbook_ensemble = Ensemble(list_of_models=[model1,model2,model3],
                                        identifier='textbook_ensemble')
    return textbook_ensemble

def test_fva_return_dims():
    ensemble = construct_textbook_ensemble()
    # get the exchange reactions to perform FVA
    ex_rxns = [rxn.id for rxn in \
                ensemble.base_model.reactions if rxn.id.startswith('EX')]
    fva_fluxes = ensemble_fva(ensemble,reaction_list=ex_rxns)
    rows = fva_fluxes.shape[0]
    columns = fva_fluxes.shape[1]
    assert rows == (len(ensemble.members)*2) # two rows for each model
    assert columns == (len(ex_rxns) + 1) # one additional column for model name
