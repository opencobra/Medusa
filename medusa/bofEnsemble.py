from cobra.core import DictList
from medusa.core.ensemble import Ensemble
from medusa.core.member import Member
from medusa.core.feature import Feature
import warnings
import numpy as np
import pandas as pd

def bofEnsemble(model, BofDf, BofId=None):
    '''
    Create an ensemble of models where each member differs only with respect to
    their biomass objective functions (BOF), without the need of first constructing 
    a list of the individual models. Note that it is allowed for members to have
    a coefficient of zero for certain metabolites in the BOF, i.e., removing the
    metabolite from the BOF. However, it is currently not possible for members to
    include additional metabolites compared to the baseline model.

    Parameters
    ----------
    model : cobra.Model
        The ensemble with which to perform reaction deletions
    BofDf : A pandas.DataFrame
        A dataframe in which each row represents a metabolite within
        the BOF. The columns correspond to the different ensemble members, where
        the first column corresponding to the input (baseline) model. The values
        in the dataframe are the coefficients of that metabolite in that ensemble
        member. As a reference, a valid BofDf dataframe can be generated using the
        helper function _getBofDf.
    BofId : str, optional
        Identifier of the biomass objective function reaction.
        If not provided (None, default), the functions attempts to retrieve 
        the BOF automatically. 

    Returns
    -------
    Medusa.core.ensemble
        An ensemble where each member has accordingly adjusted BOF coefficients
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
    
    metDict = ensemble.base_model.reactions.get_by_id(BofId).metabolites
    states = {}
    for col_idx in range(BofDf.shape[1]):
        col_name = BofDf.columns[col_idx]
        states[col_name] = _update_dict_from_df_column(metDict, BofDf, col_idx)

    # TODO if we ever want to allow for new members having additional metabolites in their 
    # BOF compared to the baseline model, we should probably update here.
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
        names = [ensemble.base_model.name] + ['placeholderName' for _ in range(BofDf.shape[1]-1)]
    else:
        names = ["templateName"] + ['placeholderName' for _ in range(BofDf.shape[1])]

    if ensemble.base_model.id is not None:
        ids = [ensemble.base_model.id] + [f'model_{i}' for i in range(BofDf.shape[1]-1)]
    else:
        ids = ["templateId"] + ['placeholderId' for _ in range(BofDf.shape[1])]

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

def _update_dict_from_df_column(metDict, df, col_idx):
    values = df.iloc[:, col_idx].tolist()
    updated_dict = {}
    for key, value in zip(metDict.keys(), values):
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

    '''
    Create a dataframe in which each row represents a metabolite within
    the BOF. The columns correspond to the different ensemble members, where
    the first column corresponding to the input (baseline) model. The values
    in the dataframe are the coefficients of that metabolite in that ensemble
    member. In this (toy) function, coefficients for each model are the sample 
    from a normal distribution with as mean the metabolite coefficient in the 
    original BOF and standard deviation 0.1.

    Parameters
    ----------
    model : cobra.Model
        The ensemble with which to perform reaction deletions
    BofId : str, optional
        Identifier of the biomass objective function reaction.
        If not provided (None, default), the functions attempts to retrieve 
        the BOF automatically.
    n_models : Int
        Define how many ensemble members should be created.

    Returns
    -------
    A pandas Dataframe that serves as input for the bofEnsemble() function.
    '''
    
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
    BofDf = pd.DataFrame(model_data, index=met_ids, columns=col_names)
    BofDf.insert(0, model.id, template)

    return BofDf