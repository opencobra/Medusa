# ensemble_model class

from __future__ import absolute_import

from cobra.core.object import Object
from cobra.core import Model
from cobra.core import DictList
from cobra.core import Reaction

from medusa.core.member import Member
from medusa.core.feature import Feature

from pickle import dump

import cobra
import random
import pandas as pd

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']
MISSING_ATTRIBUTE_DEFAULT = {'lower_bound':0,'upper_bound':0}

class Ensemble(Object):
    """
    Ensemble of metabolic models

    Parameters
    ----------
    identifier : string
        The identifier to associate with the ensemble as a string.

    list_of_models : list of cobra.core.model.Model
        Either a list of existing Model objects in which case a new Model
        object is instantiated and an ensemble is constructed using the list of
        Models, or None/empty list, in which case an ensemble is created with
        empty attributes.

    name : string
        Human-readable name for the ensemble

    Attributes
    ----------
    base_model : Model
        A cobra.core.Model that contains all variable and invariable components
        of an ensemble.
    members : DictList
        A DictList where the key is the member identifier and the value is a
        medusa.core.member.Member object
    features : DictList
        A DictList where the key is the feature identifier and the value is a
        medusa.core.feature.Feature object
    """
    def __init__(self,list_of_models=[], identifier=None, name=None):
        Object.__init__(self,identifier,name)
        if len(list_of_models) > 1:
            if not all(isinstance(x, Model) for x in list_of_models):
                raise AttributeError("list_of_models may only contain cobra.core.Model objects")
            if len([model.id for model in list_of_models]) > \
                            len(set([model.id for model in list_of_models])):
                raise AssertionError("Ensemble members cannot have duplicate model ids.")
            self.features = DictList()
            self._populate_features_base(list_of_models)

            self.members = DictList()
            self._populate_members(list_of_models)

        else:
            if len(list_of_models) == 0:
                self.base_model = Model(id_or_model=identifier+'_base_model',\
                                        name=name)
            else:
                if not isinstance(list_of_models[0], Model):
                    raise AttributeError("list_of_models may only contain cobra.core.Model objects")
                self.base_model = list_of_models[0]

    def _populate_features_base_new(self, list_of_models):
        # Determine all reactions across all models and construct the base model
        base_model = list_of_models[0].copy()
        all_reactions = set(rxn.id for rxn in base_model.reactions)
        for model in list_of_models:
            model_rxn_ids = set(rxn.id for rxn in model.reactions)
            new_reactions = model_rxn_ids - all_reactions
            if new_reactions:
                reactions_to_add = [model.reactions.get_by_id(rxn_id) for rxn_id in new_reactions]
                base_model.add_reactions(reactions_to_add)
                all_reactions.update(new_reactions)
        all_reactions = list(all_reactions)

        # Pre-cache model reaction attributes in dicts to avoid repeated getattr calls
        model_reaction_attrs = {}
        for model in list_of_models:
            mid = model.id
            model_reaction_attrs[mid] = {}
            for rxn in model.reactions:
                # Direct attribute access instead of getattr
                model_reaction_attrs[mid][rxn.id] = {
                    'lower_bound': rxn.lower_bound,
                    'upper_bound': rxn.upper_bound
                }

        # Iterate over each reaction to detect variable attributes
        for reaction in all_reactions:
            # Collect reaction attributes per model
            rxn_vals = {}
            for model in list_of_models:
                mid = model.id
                rxn_attrs = model_reaction_attrs[mid]
                if reaction in rxn_attrs:
                    rxn_vals[mid] = rxn_attrs[reaction]
                else:
                    # Use default bounds if reaction missing from model
                    rxn_vals[mid] = {'lower_bound': 0, 'upper_bound': 0} # TODO remark: could use MISSING_ATTRIBUTE_DEFAULT instead

            # Create a DataFrame for easier analysis of attribute variability
            rxn_vals = pd.DataFrame.from_dict(rxn_vals, orient='index') # TODO remark: faster than transposing

            for reaction_attribute in REACTION_ATTRIBUTES:
                if rxn_vals[reaction_attribute].nunique() > 1:
                    rxn_from_base = base_model.reactions.get_by_id(reaction)
                    feature_id = f"{reaction}_{reaction_attribute}"
                    feature_id = rxn_from_base.id + '_' + reaction_attribute

                    # Create states dict for feature
                    states = rxn_vals[reaction_attribute].to_dict()

                    # Create and add the Feature object
                    feature = Feature(
                        ensemble=self,
                        identifier=feature_id,
                        name=rxn_from_base.name,
                        base_component=rxn_from_base,
                        component_attribute=reaction_attribute,
                        states=states,
                    )
                    self.features.append(feature)

        self.base_model = base_model

    def _populate_members(self,list_of_models):
        for model in list_of_models:
            model_states = dict()
            for feature in self.features:

                model_states[feature] = feature.get_model_state(model.id)
            member = Member(ensemble=self,\
                            identifier=model.id,\
                            name=model.name,\
                            states=model_states)

            self.members += [member]


    def set_state(self,member):
        """Set the state of the base model to represent a single member.

        Sets all features to the state for the provided member. Only
        reaction states are currently implemented (e.g. GPRs as features
        will not work)

        Parameters
        ----------
        member : str or medusa.Member
            The Member.id, or the Member object itself, to set the state
            of the Ensemble.base_model to represent.

        """
        # if member was passed as an id, get the actual member object
        if isinstance(member, str):
            member = self.members.get_by_id(member)

        for feature in self.features:
            if isinstance(feature.base_component, cobra.core.Reaction):
                setattr(feature.base_component,\
                        feature.component_attribute,\
                        feature.states[member.id])
            else:
                raise AttributeError("Only cobra.core.Reaction supported for base_component type")

    def to_pickle(self, filename):
        """
        Save an ensemble as a pickled object. Pickling is currently the only supported
        method for saving and loading ensembles.

        Parameters
        ----------
        filename : String
            location to save the pickle.
        """

        with open(filename, "wb") as outfile:
            dump(self, outfile, protocol=4)

    def extract_member(self, member):
        """
        Extract an individual member as a cobrapy model (cobra.Model), removing
        any components associated with features that are inactive in member.

        Provided as a more convenient option than medusa.Member.to_model(),
        but is the exact same.

        Parameters
        ----------
        member : str or medusa.Member
            The Member.id, or the Member object itself, to be represented in
            the cobrapy model output.

        Returns
        -------
        model : cobra.Model
            The extracted member as a cobrapy model.
        """
        # if member was passed as an id, get the actual member object
        if isinstance(member, str):
            member = self.members.get_by_id(member)

        model = member.to_model()
        return model
