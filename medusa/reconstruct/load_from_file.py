import medusa
import cobra
import cobra.test
import math
import copy
import random
from copy import deepcopy

from medusa.core.feature import Feature

def parent_attr_of_base_component(base_comp):
    """
    Output a string to indicate the parent attribute of the cobra.core object.
    
    Parameters
    ----------
    base_comp : cobra.core object
        Ensemble base_component of feature. i.e. cobra reaction, metabolite, or gene
    """
    
    if type(base_comp) == cobra.core.Reaction:
        parent_attr = "reactions"
    elif type(base_comp) == cobra.core.Metabolite:
        parent_attr = "metabolites"
    elif type(base_comp) == cobra.core.Gene:
        parent_attr = "genes"
    else:
        raise AttributeError("Only cobra.core.Reaction, cobra.core.Metabolite, and cobra.core.Gene supported for base_component type")
    
    return parent_attr
    
def batch_load_from_files(model_file_names, identifier='ensemble', batchsize=5, verbose = False):
    
    """
    Loads a list of models my file name in batches to and generates an ensemble object.
    This function is meant to be used to limit how much flash memory is required to 
    generate very large ensembles. 
    
    Parameters
    ----------
    model_jsons : List 
        List of json cobra.model file names.
    batchsize : Integer
        Total number of models loaded into memory.
    """
    
    total = len(model_file_names)
    range_lists = []
    iterations = math.ceil(total/batchsize)
    fix_last = 0
    for i in range(iterations):
        start = batchsize*i
        stop = batchsize*(1+i)
        if stop > total:
            stop = total
        if len(range(start,stop)) == 1:
            fix_last = 1
        range_lists.append(list(range(start,stop)))
    if fix_last == 1:
        range_lists[iterations-1] = [range_lists[iterations-2][batchsize-1]] + range_lists[iterations-1]
        del range_lists[iterations-2][batchsize-1]

    for range_list in range_lists:
        model_list = []
        for model_file in [model_file_names[i] for i in range_list]:
            if model_file.endswith(".json"):
                model = cobra.io.load_json_model(model_file)
            elif model_file.endswith(".xml"):
                model = cobra.io.read_sbml_model(model_file)
            else:
                raise AttributeError("Only .json or .xml files supported")
            if not isinstance(model, cobra.core.Model):
                raise AttributeError("Only cobra.core.Model objects supported")
            model_list.append(model)
        if range_list[0] == 0:
            final_ensemble = medusa.core.Ensemble(model_list, identifier = identifier)
            if verbose == True:
                print("Ensemble 1 finished")
        else:
            new_ensemble = medusa.core.Ensemble(model_list)
            final_ensemble = add_ensembles(final_ensemble,new_ensemble)
            del new_ensemble
            if verbose == True:
                print("Next ensemble finished")
        del model_list 
        del model
    return final_ensemble

def add_ensembles(e1,e2,verbose=False):
    
    """
        Adds two ensemble objects together.
        
        Parameters
        ----------
        e1 & e2 : Ensemble Objects
            Generated using medusa.core.Ensemble()
    """
    
    # Deep copy ensembles
    emodel1 = e1
    emodel2 = e2
    emodel3 = copy.deepcopy(e1)

    # Add reactions to new base_model: Base_model1 + Base_model2 = base_model3
    base_model = copy.deepcopy(emodel1.base_model)
    all_reactions = set()
    all_reactions = set([rxn.id for rxn in base_model.reactions])
    new_reactions = set([rxn.id for rxn in emodel2.base_model.reactions]) - all_reactions
    reactions_to_add = [emodel2.base_model.reactions.get_by_id(rxn) for rxn in new_reactions]
    base_model.add_reactions(reactions_to_add)
    emodel3.base_model = base_model

    all_feats = set()
    all_feats = set([feat.id for feat in emodel1.features])

    new_feats = set()
    new_feats = set([feat.id for feat in emodel2.features]) - all_feats

    old_feats = set()
    old_feats = all_feats - new_feats

    # Add new features to base ensemble
    for feat_id in new_feats:
        feat = emodel2.features.get_by_id(feat_id)
        emodel3.features = emodel3.features + [feat]

    # Swap feature objects to make consistent across ensemble
    emodel3.members = emodel1.members + emodel2.members

    em1_feats = set([feat.id for feat in emodel1.features])
    em2_feats = set([feat.id for feat in emodel2.features])
    
    for member in emodel2.members:
        member_obj = emodel3.members.get_by_id(member.id)
        new_feat_dict = dict()
        for feat, val in member_obj.states.items():
            if feat.id in em1_feats:
                base_feat = emodel1.features.get_by_id(feat.id)
                new_feat_dict[base_feat] = val
            else:
                new_feat_dict[feat] = val
        member_obj.states = new_feat_dict
        
    # Make feature.base_components consistent with base_model
        # This may need to be updated to account for metabolite and gene features
    for feat_obj in emodel3.features:
        if isinstance(feat_obj.base_component, cobra.core.Reaction):
            rxn_base_obj = emodel3.base_model.reactions.get_by_id(feat_obj.base_component.id)
            feat_obj.base_component = rxn_base_obj
    
    em1_rxns = set([rxn.id for rxn in emodel1.base_model.reactions])
    em2_rxns = set([rxn.id for rxn in emodel2.base_model.reactions])
    
    # Create features for reactions missing from either base_model without an existing feature
    missing_rxns = (em2_rxns - em1_rxns) | (em1_rxns - em2_rxns)
    exist_feat_ids = set([feat.base_component.id for feat in emodel3.features])
    attr_list = ['lower_bound','upper_bound']
    states1 = emodel1.features[0].states
    states2 = emodel2.features[0].states
    for rxn_id in missing_rxns:
        if not rxn_id in exist_feat_ids:
            for attr_str in attr_list:
                if rxn_id in em1_rxns:
                    attr1 = getattr(getattr(emodel1.base_model, "reactions").get_by_id(rxn_id), attr_str)
                else:
                    attr1 = 0.0
                if rxn_id in em2_rxns:
                    attr2 = getattr(getattr(emodel2.base_model, "reactions").get_by_id(rxn_id), attr_str)
                else:
                    attr2 = 0.0
                rxn_from_base = emodel3.base_model.reactions.get_by_id(rxn_id)
                feature_id = rxn_from_base.id + '_' + attr_str
                states1 = dict.fromkeys(states1, attr1)
                states2 = dict.fromkeys(states2, attr2)
                states = dict(states1, **states2)
                feature = Feature(ensemble=emodel3,\
                        identifier=feature_id,\
                        name=rxn_from_base.name,\
                        base_component=rxn_from_base,\
                        component_attribute=attr_str,\
                        states=states)
                emodel3.features = emodel3.features + [feature]
                if verbose == True:
                    print("New feature added: " + feature_id)

    # Check for new features that need to be added because the base models don't align
        # Needs to be generalized to all feature types beyond reactions
    ovrlp_rxns = (em1_rxns & em2_rxns) - exist_feat_ids

    attr_list = ['lower_bound','upper_bound']
    states1 = emodel1.features[0].states
    states2 = emodel2.features[0].states

    for rxn_id in ovrlp_rxns:
        for attr_str in attr_list:
            attr1 = getattr(getattr(emodel1.base_model, "reactions").get_by_id(rxn_id), attr_str)
            attr2 = getattr(getattr(emodel2.base_model, "reactions").get_by_id(rxn_id), attr_str)
            if attr1 != attr2:
                rxn_from_base = emodel3.base_model.reactions.get_by_id(rxn_id)
                feature_id = rxn_from_base.id + '_' + attr_str
                states1 = dict.fromkeys(states1, attr1)
                states2 = dict.fromkeys(states2, attr2)
                states = dict(states1, **states2)
                feature = Feature(ensemble=emodel3,\
                        identifier=feature_id,\
                        name=rxn_from_base.name,\
                        base_component=rxn_from_base,\
                        component_attribute=attr_str,\
                        states=states)
                emodel3.features = emodel3.features + [feature]
                if verbose == True:
                    print("New feature added: " + feature_id)
    
    # Set feature.states
    for feature_obj in emodel3.features:
        dict1_zeros = dict.fromkeys(emodel1.features[0].states, 0.0)
        dict2_zeros = dict.fromkeys(emodel2.features[0].states, 0.0)
        if feature_obj.id in em1_feats:
            dict1 = emodel1.features.get_by_id(feature_obj.id).states
        elif feature_obj.base_component.id in em1_rxns:
            base_comp = feature_obj.base_component
            comp_attr = feature_obj.component_attribute
            parent_attribute = parent_attr_of_base_component(base_comp)
            attr = getattr(getattr(emodel1.base_model, parent_attribute).get_by_id(base_comp.id), comp_attr)
            for member in emodel1.members:
                dict1_zeros[member.id] = attr
            dict1 = dict1_zeros
        else:
            dict1 = dict1_zeros
        if feature_obj.id in em2_feats:
            dict2 = emodel2.features.get_by_id(feature_obj.id).states
        elif feature_obj.base_component.id in em2_rxns:
            base_comp = feature_obj.base_component
            comp_attr = feature_obj.component_attribute
            attr = getattr(getattr(emodel2.base_model, "reactions").get_by_id(base_comp.id), comp_attr)
            for member in emodel2.members:
                dict2_zeros[member.id] = attr
            dict2 = dict2_zeros
        else:
            dict2 = dict2_zeros

        feature_obj.states = dict(dict1, **dict2)
        
    # Set member.states
    for member in emodel3.members:
        member.ensemble = emodel3
        temp_dict = dict()
        for feat in emodel3.features:
            temp_dict[feat] = feat.states[member.id]
        member.states = temp_dict
    
    return emodel3
