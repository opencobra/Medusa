from cobra.core import DictList
from medusa.core.ensemble import Ensemble
from medusa.core.member import Member
from medusa.core.feature import Feature
import warnings
import numpy as np
import pandas as pd
import itertools

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']

def boundsEnsemble(model, boundsDict):
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

    # Set features, similar to _populate_features_base()
    for rxn_id in boundsDict.keys():
        for attr in REACTION_ATTRIBUTES:
            if boundsDict[rxn_id][attr].nunique() > 1:
                rxn_base = ensemble.base_model.reactions.get_by_id(rxn_id)
                feature_id = f"{rxn_id}_{attr}"

                # Create states dict for feature
                states = boundsDict[rxn_id][attr].to_dict()

                # Create and add the Feature object
                feature = Feature(
                    ensemble=ensemble,
                    identifier=feature_id,
                    name=rxn_base.name,
                    base_component=rxn_base,
                    component_attribute=attr,
                    states=states,
                )
                ensemble.features.append(feature)
    
    names = [ensemble.base_model.name] + ['placeholderName' for _ in range(len(boundsDict[next(iter(boundsDict))].index) - 1)]
    ids = list(boundsDict[next(iter(boundsDict))].index)

    # Populate members, similar to _populate_members()
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

def _setBoundsRandom(model, rxn_ids, bound=None, reversibility='respect', n_models=100):
    if bound is None:
        bound = 1000
        warnings.warn("No 'bound' provided for method 'random'. Defaulting to bound = 1000.")

    if not isinstance(bound, (int, float)):
        raise TypeError(f"'bound' must be a float or int, got {type(bound)}")

    ids = [f'model_{i}' for i in range(n_models)]
    
    boundsDict = {
        rxn_id: pd.DataFrame({
            'lower_bound': np.random.uniform(-abs(bound), 0, size=len(ids)),
            'upper_bound': np.random.uniform(0, abs(bound), size=len(ids))
        }, index=ids)
        for rxn_id in rxn_ids
    }

    if reversibility == 'respect':
        for rxn_id in rxn_ids:
            if model.reactions.get_by_id(rxn_id).lower_bound == 0:
                boundsDict[rxn_id].lower_bound = 0
            if model.reactions.get_by_id(rxn_id).upper_bound == 0:
                boundsDict[rxn_id].upper_bound = 0

    return boundsDict


def _setBoundsOnOff(model, rxn_ids):

    combinations = list(itertools.product([0,1], repeat=len(rxn_ids)))
    df = pd.DataFrame(combinations, columns=rxn_ids)
    df.index = [f"model_{i}" for i in range(len(df))]

    boundsDict = {}
    for rxn_id in df.columns:
        values = df[rxn_id]
        boundsDict[rxn_id] = pd.DataFrame({
            'lower_bound': values * model.reactions.get_by_id(rxn_id).lower_bound + 0,
            'upper_bound': values * model.reactions.get_by_id(rxn_id).upper_bound
        }, index=df.index)

    return boundsDict


def _setBoundsFullFactorial(model, rxn_ids, bound=None, reversibility='respect'):
    if bound is None:
        bound = 1000
        warnings.warn("No 'bound' provided for method 'fullFactorial'. Defaulting to bound = 1000.")

    default_options = [(-abs(bound), 0), (-abs(bound), abs(bound)), (0, 0), (0, abs(bound))]

    if reversibility == 'ignore':    
        bound_options_dict = {rxn_id: default_options for rxn_id in rxn_ids}
    else:
        bound_options_dict = {}
        for rxn_id in rxn_ids:
            if model.reactions.get_by_id(rxn_id).reversibility:
                bound_options_dict[rxn_id] = default_options
            else:
                bound_options_dict[rxn_id] = [(0, 0), 
                                              tuple(bound * (x / abs(x)) if x != 0 else 0 for x in model.reactions.get_by_id(rxn_id).bounds)]

    reaction_options = [bound_options_dict[rxn_id] for rxn_id in rxn_ids]
    all_combinations = list(itertools.product(*reaction_options))
    row_labels = [f"model_{i}" for i in range(len(all_combinations))]
    
    boundsDict = {}
    for i, rxn_id in enumerate(rxn_ids):
        lower_bounds = [combo[i][0] for combo in all_combinations]
        upper_bounds = [combo[i][1] for combo in all_combinations]
        boundsDict[rxn_id] = pd.DataFrame({
            'lower_bound': lower_bounds,
            'upper_bound': upper_bounds
        }, index=row_labels)

    return boundsDict

def _setBounds(model, rxn_ids, method='random', reversibility=None, bound=None, n_models=100):
    allowed_methods = ['random', 'onOff', 'fullFactorial']
    if method not in allowed_methods:
        raise ValueError(f"Invalid method '{method}'. Choose one of {allowed_methods}.")

    allowed_reversibility = ['respect', 'ignore']
    if reversibility is not None and reversibility not in allowed_reversibility:
        raise ValueError(f"Invalid reversibility '{reversibility}'. Choose one of {allowed_reversibility}.")

    # Set defaults depending on method
    if method == 'random':
        if reversibility is None:
            reversibility = 'respect'
        boundsDict = _setBoundsRandom(model, rxn_ids, bound=bound, reversibility=reversibility, n_models=n_models)

    elif method == 'onOff':
        if bound is not None:
            warnings.warn("'bound' argument is ignored when method is 'onOff'")
        if reversibility is not None:
            warnings.warn("'reversibility' argument is ignored when method is 'onOff'")
        if n_models != 100:
            warnings.warn("'n_models' argument is ignored when method is 'onOff'")
        boundsDict = _setBoundsOnOff(model, rxn_ids)

    elif method == 'fullFactorial':
        if n_models != 100:
            warnings.warn("'n_models' argument is ignored when method is 'fullFactorial'")
        if bound is None:
            bound = 1000
            warnings.warn("No 'bound' provided for method 'fullFactorial'. Defaulting to bound = 1000.")
        if reversibility is None:
            reversibility = 'respect'
            warnings.warn("No 'reversibility' provided for method 'fullFactorial'. Defaulting to reversibility = 'respect'.")
        boundsDict = _setBoundsFullFactorial(model, rxn_ids, bound=bound, reversibility=reversibility)

    # Add row for base model at the top
    for rxn_id in rxn_ids:
        boundsDict[rxn_id] = pd.concat([
            pd.DataFrame({
                'lower_bound': [model.reactions.get_by_id(rxn_id).lower_bound],
                'upper_bound': [model.reactions.get_by_id(rxn_id).upper_bound]
            }, index=[model.id]),
            boundsDict[rxn_id]
        ])

    return boundsDict