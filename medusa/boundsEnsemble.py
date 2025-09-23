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
    Create an ensemble of models where all members have the same reactions but
    different reaction bounds, without the need of first constructing a list of these
    individual models. While all members share the same reactions, the bounds of some
    reactions in some members maybe set to (0,0), i.e., inactivating the reaction.

    Parameters
    ----------
    model : cobra.Model
        A cobrapy model.
    boundsDict : A dictionary of pandas.DataFrames
        A dictionary with reaction identifiers as keys and dataframes as values.
        Each row in the dataframe is the identifier of a model within the resulting ensemble.
        The dataframe has two columns, one for the lower and one for the upper bound of the reaction.
        As a reference, a valid boundsDict can be generated using the
        helper function _setBounds.

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

def _setBoundsRandom(model, rxn_ids, bound=None, reversibility=True, n_models=100):
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

    if reversibility:
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


def _setBoundsFullFactorial(model, rxn_ids, bound=None, reversibility=True):
    if bound is None:
        bound = 1000
        warnings.warn("No 'bound' provided for method 'fullFactorial'. Defaulting to bound = 1000.")

    default_options = [(-abs(bound), 0), (-abs(bound), abs(bound)), (0, 0), (0, abs(bound))]

    if reversibility:
        bound_options_dict = {}
        for rxn_id in rxn_ids:
            if model.reactions.get_by_id(rxn_id).reversibility:
                bound_options_dict[rxn_id] = default_options
            else:
                bound_options_dict[rxn_id] = [(0, 0), 
                                              tuple(bound * (x / abs(x)) if x != 0 else 0 for x in model.reactions.get_by_id(rxn_id).bounds)]
    else:    
        bound_options_dict = {rxn_id: default_options for rxn_id in rxn_ids}

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

def _setBounds(model, rxn_ids, method='random', reversibility=True, bound=None, n_models=100):

    ''' 
    Helper function for constructing object 'boundsDict', an argument used by 
    the 'boundsEnsemble' function. Users can provide their own boundsDict, this 
    helper function exists only for providing users with reusable code and for 
    internal function testing. 

    Parameters
    ----------
    model : cobra.Model A single cobraPy Model that will be used as a baseline
        to generate an ensemble of models with different reaction bounds for a select
        set of reactions.
    
    rxn_ids : list of str 
        Target reactions for which each ensemble member will have a different 
        combination of bounds.

    method : str, optional 
        Method used to generate reaction bounds.
        Must be one of: 
            - 'random' : The value of the lower and upper bound of each target reaction 
                is chosen as a random integer between 0 and the bound provided in 
                the 'bound' argument.
            - 'onOff' : Sets reactions to either completely off (0, 0) or active (-bound, +bound). 
                The number of models that will be created is two to the power of the number of
                target reactions (length of rxn_ids argument).
            - 'fullFactorial' : Creates all possible combinations of states for all target reactions.
                I.e., each reaction can take on values (0,0), (-bound,0), (0,bound), and (-bound, bound).
                The number of potential combinations is thus four to the power of the number of
                target reactions (length of rxn_ids argument). Note, that if reversibility is 'respect',
                some options will be removed, i.e., the reversibility direction of the reaction in the
                baseline model will be respected.
        Default is 'random'. 

    reversibility : boolean, optional 
        If True, the reaction reversibility of the baseline model will be respected. 
        E.g., if the baseline model only allows for the forward reaction,
        then reversibility or backward reactions will not be allowed for any of the new members.
        If False, previously irreversible reactions will be allowed to be reversible.
        Default is True.

    bound : float or int or None, optional 
        This value is used as the magnitude of the reaction bounds (e.g., ±bound). 
        If None (default), a bound of 1000 is used internally.

    n_models : int
        Number of models that will be created if method is 'random', ignored otherwise.
        Default is 100.

    Returns
    -------
    boundsDict : A dictionary of pandas.DataFrames A dictionary of dataframes in which
        each row (index) represents a model within the ensemble, and each column represents
        a reaction for which values of objective when the reaction is deleted are returned.

    '''

    allowed_methods = ['random', 'onOff', 'fullFactorial']
    if method not in allowed_methods:
        raise ValueError(f"Invalid method '{method}'. Choose one of {allowed_methods}.")

    # Set defaults depending on method
    if method == 'random':
        boundsDict = _setBoundsRandom(model, rxn_ids, bound=bound, reversibility=reversibility, n_models=n_models)

    elif method == 'onOff':
        if bound is not None:
            warnings.warn("'bound' argument is ignored when method is 'onOff'")
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