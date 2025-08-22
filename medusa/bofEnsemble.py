from cobra.core import DictList
from medusa.core.ensemble import Ensemble
from medusa.core.member import Member
from medusa.core.feature import Feature
import warnings
import numpy as np
import pandas as pd

def bofEnsemble(model, BofDf, BofId=None):
    '''
    Create an ensemble of models where each member has the same reactions but
    different reaction bounds, without the need of first constructing a list of these
    individual models. While all members share the same reactions, the bounds of some
    reactions in some members maybe set to (0,0), i.e., inactivating the reaction.

    Parameters
    ----------
    model : cobra.Model
        The ensemble with which to perform reaction deletions
    boundsDict : A dictionary of pandas.DataFrames
        A dictionary of dataframes in which each row (index) represents a 
        model within the ensemble, and each column represents a reaction for 
        which values of objective when the reaction is deleted are returned.

    Returns
    -------
    Medusa.core.ensemble
        An ensemble where each member has accordingly adjusted reaction bounds
    '''
    # Setup ensemble structure base on single baseline model
    ensemble = Ensemble([model],
                        identifier = "placeholderId",
                        name = "placeholderName")
    ensemble.features = DictList()
    ensemble.members = DictList()

    # Attempt to retreive BofId automatically if not provided
    if BofId is None:
        BofId = _getBofId(model)
    
    hlp = ensemble.base_model.reactions.get_by_id(BofId).metabolites
    states = {}
    for col_idx in range(1, BofDf.shape[1]):  # all columns except the first
        col_name = BofDf.columns[col_idx]
        states[col_name] = _update_dict_from_df_column(hlp, BofDf, col_idx)

    rxn_base = ensemble.base_model.reactions.get_by_id(BofId)
    feature_id = f"{BofId}_metabolites"

    # Create and add the Feature object
    feature = Feature(
        ensemble=ensemble,
        identifier=feature_id,
        name=rxn_base.name,
        base_component=rxn_base,
        component_attribute='metabolites',
        states=states,
    )
    ensemble.features.append(feature)

    if ensemble.base_model.name is not None:
        names = [ensemble.base_model.name] + ['placeholderName' for _ in range(BofDf.shape[1]-2)]
    else:
        names = ['placeholderName' for _ in range(BofDf.shape[1]-1)]

    if ensemble.base_model.id is not None:
        ids = [ensemble.base_model.id] + [f'model_{i}' for i in range(BofDf.shape[1]-2)]
    else:
        ids = ['placeholderId' for _ in range(BofDf.shape[1]-1)]

    for i in range(0,len(ids)):
        model_states = dict()
        for feature in ensemble.features:
            model_states[feature] = feature.get_model_state(ids[i])
        member = Member(ensemble=ensemble,\
                        identifier=ids[i],\
                        name=names[i],\
                        states=model_states)

        ensemble.members += [member]

    return ensemble

def _update_dict_from_df_column(hlp, df, col_idx):
    values = df.iloc[:, col_idx].tolist()
    updated_dict = {}
    for key, value in zip(hlp.keys(), values):
        updated_dict[key] = value
    return updated_dict

def _getBofId(model):

    # Get BofId
    warnings.warn("No 'BofId' provided, attempting to obtain programmatically.")
    BOF_candidates = [var.name for var in model.objective.expression.free_symbols]
    if not BOF_candidates:
        raise ValueError("No objective function reaction found.")
    BofId = BOF_candidates[0]
    print(f"Identified BofId: {BofId}")

    # Assess automatic BofId retrieval
    try:
        reaction = model.reactions.get_by_id(BofId)
    except KeyError as e:
        raise ValueError(
            f"Reaction '{BofId}' not found in model. "
            "Automatic retrieval of BOF failed. "
            "Please provide a valid BOF_id manually."
        ) from e
    return BofId

def _getBofDf(model, BofId=None, n_models=100):
    
    # Attempt to retreive BofId automatically if not provided
    if BofId is None:
        BofId = _getBofId(model)

    # Extract metabolites and coefficients
    reaction = model.reactions.get_by_id(BofId)
    metabolites = reaction.metabolites  # dict: metabolite -> coefficient
    met_ids = [met.id for met in metabolites.keys()]
    template = np.array([coef for coef in metabolites.values()])

    # Prepare the matrix of all model values (template + noise)
    noise = np.random.normal(loc=0, scale=0.1, size=(len(met_ids), n_models))
    model_data = template[:, np.newaxis] + noise  # shape: (n_mets, n_models)

    # Create full DataFrame at once
    col_names = [f"model_{i}" for i in range(n_models)]
    BofDf = pd.DataFrame(model_data, columns=col_names)
    BofDf.insert(0, model.id, template)
    BofDf.insert(0, "Metabolite", met_ids)

    return BofDf