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
        empty attributes. When `features` is provided, this must contain
        exactly one Model, which is used as the base model.

    name : string
        Human-readable name for the ensemble

    features : list of medusa.core.feature.Feature, optional
        Pre-built Feature objects describing how members vary from the base
        model. When provided, `list_of_models` must contain exactly one Model.
        Member ids are taken from the keys of each feature's `states` dict;
        every supplied feature must define states for exactly the same set of
        member ids. Each feature's base_component is re-resolved against the
        base model so set_state mutates the right reaction object.

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
    def __init__(self, list_of_models=None, identifier=None, name=None,
                 features=None):
        Object.__init__(self,identifier,name)

        # Defaulting to None rather than [] so the default is not a shared
        # mutable object.
        if list_of_models is None:
            list_of_models = []

        if features is not None and len(features) == 0:
            raise ValueError(
                "`features` must be a non-empty list of Feature objects."
            )
        if features:
            if len(list_of_models) != 1:
                raise AttributeError(
                    "When `features` is provided, `list_of_models` must "
                    "contain exactly one cobra.core.Model."
                )
            if not isinstance(list_of_models[0], Model):
                raise AttributeError(
                    "list_of_models may only contain cobra.core.Model objects"
                )
            self.base_model = list_of_models[0]
            self._attach_prebuilt_features(features)
            return

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

    def _populate_features_base(self, list_of_models):
        # Determine all reactions across all models and construct the base model
        base_model = list_of_models[0].copy()
        all_reactions = set(rxn.id for rxn in base_model.reactions)
        for model in list_of_models:
            model_rxn_ids = set(rxn.id for rxn in model.reactions)
            new_reactions = model_rxn_ids - all_reactions
            if new_reactions:
                # Copy before adding. cobra's add_reactions takes ownership of
                # the Reaction objects it is given: it reassigns each
                # reaction._model and remaps its metabolites onto base_model's
                # copies. Handing it live reactions therefore removes them
                # from the caller's model, so set_state would go on to mutate
                # the caller's input and that model would be left structurally
                # inconsistent. Only list_of_models[0] was protected, by the
                # .copy() above.
                reactions_to_add = [
                    model.reactions.get_by_id(rxn_id).copy()
                    for rxn_id in sorted(new_reactions)]
                base_model.add_reactions(reactions_to_add)
                all_reactions.update(new_reactions)
        # Sorted, not just list(). Iterating the set directly makes the order
        # of self.features depend on string hash randomization, so two
        # ensembles built from identical inputs come out with identically
        # valued but differently ordered features between processes.
        all_reactions = sorted(all_reactions)

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
                if rxn_vals[reaction_attribute].nunique() > 1: # TODO remark: used to be len() > 1, seems incorrect
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

    def _attach_prebuilt_features(self, features):
        self.features = DictList()
        self.members = DictList()
        base_reaction_ids = {rxn.id for rxn in self.base_model.reactions}
        for feature in features:
            component = feature.base_component
            if not isinstance(component, cobra.core.Reaction):
                raise AttributeError(
                    "Only cobra.core.Reaction is supported for "
                    "feature.base_component"
                )
            if component.id not in base_reaction_ids:
                raise ValueError(
                    f"Feature '{feature.id}' references reaction "
                    f"'{component.id}', which is not in the base model."
                )
            # Copy rather than rewire in place. Both `base_component` and
            # `ensemble` are rebound below, so attaching the caller's own
            # Feature objects would repoint them at this ensemble; passing the
            # same feature list to two ensembles used to leave the first one
            # holding features bound to the second's base model, and
            # set_state would then mutate the wrong model. The states dict is
            # copied too, shallowly, so that adding a member to one ensemble
            # does not appear in the other.
            attached = Feature(
                identifier=feature.id,
                name=feature.name,
                ensemble=self,
                base_component=self.base_model.reactions.get_by_id(
                    component.id),
                component_attribute=feature.component_attribute,
                states=dict(feature.states),
            )
            self.features.append(attached)

        member_ids = list(self.features[0].states.keys())
        reference_set = set(member_ids)
        for f in self.features[1:]:
            if set(f.states.keys()) != reference_set:
                raise ValueError(
                    "All supplied features must define states for exactly "
                    "the same set of member ids."
                )

        for member_id in member_ids:
            member_states = {f: f.states[member_id] for f in self.features}
            member = Member(
                ensemble=self,
                identifier=member_id,
                name=member_id,
                states=member_states,
            )
            self.members += [member]

    @classmethod
    def from_reaction_states(cls, model, reaction_id, states,
                             component_attribute='metabolites',
                             allow_new_metabolites=False,
                             identifier=None, name=None):
        """Build an ensemble whose members differ in a single reaction attribute.

        Parameters
        ----------
        model : cobra.Model
            The base model; one of its reactions will be varied across members.
        reaction_id : str
            Id of the reaction in `model` to vary. Required — no
            auto-detection from the model objective.
        states : dict
            Mapping of member_id -> attribute value. The keys become the
            ensemble's member ids; no implicit baseline is added. For the
            default component_attribute='metabolites', each value should be a
            dict mapping metabolite (id or cobra.Metabolite) to coefficient,
            suitable for passing to Reaction.add_metabolites(combine=False).
        component_attribute : str, optional
            Reaction attribute that varies across members. Defaults to
            'metabolites' (the alternative-biomass-composition use case).
        allow_new_metabolites : bool, optional
            Only meaningful when component_attribute='metabolites'. When False
            (the default), every metabolite referenced by any state must
            already be in the target reaction; otherwise a ValueError is
            raised at construction. Set True to allow members to introduce
            metabolites not in the baseline reaction (e.g. swapping
            ATP for an alternative energy carrier).
        identifier, name : str, optional
            Passed through to Ensemble.__init__.

        Returns
        -------
        Ensemble
        """
        if not isinstance(model, Model):
            raise AttributeError("`model` must be a cobra.core.Model")
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError as e:
            raise ValueError(
                f"Reaction '{reaction_id}' not found in model."
            ) from e
        if not isinstance(states, dict) or len(states) == 0:
            raise ValueError(
                "`states` must be a non-empty dict of {member_id: value}"
            )

        if component_attribute == 'metabolites' and not allow_new_metabolites:
            existing_met_ids = {met.id for met in reaction.metabolites}
            for member_id, met_dict in states.items():
                if not isinstance(met_dict, dict):
                    raise ValueError(
                        f"State '{member_id}' must be a dict of "
                        "{metabolite: coefficient}."
                    )
                for met_key in met_dict:
                    met_id = met_key.id if hasattr(met_key, 'id') else met_key
                    if met_id not in existing_met_ids:
                        raise ValueError(
                            f"State '{member_id}' references metabolite "
                            f"'{met_id}', which is not in reaction "
                            f"'{reaction_id}'. Pass "
                            "allow_new_metabolites=True to introduce new "
                            "metabolites."
                        )

        feature = Feature(
            identifier=f"{reaction_id}_{component_attribute}",
            name=reaction.name,
            base_component=reaction,
            component_attribute=component_attribute,
            states=dict(states),
        )
        return cls(
            list_of_models=[model],
            features=[feature],
            identifier=identifier,
            name=name,
        )

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

        # Bounds are collected per reaction and applied together, rather than
        # written one setattr at a time. cobra validates each assignment
        # against the bound already in place, so setting lower_bound first
        # raises whenever the new lower bound exceeds the *old* upper bound.
        # That happens for any ordinary on/off ensemble: a reaction switched
        # off holds (0, 0), and switching it back on with a positive minimum
        # flux (ATPM at 8.39, say) fails with "The lower bound must be less
        # than or equal to the upper bound". Assigning reaction.bounds as a
        # pair validates the two together.
        pending_bounds = {}
        pending_metabolites = {}

        for feature in self.features:
            component = feature.base_component
            attr = feature.component_attribute
            value = feature.states[member.id]

            if not isinstance(component, cobra.core.Reaction):
                raise AttributeError("Only cobra.core.Reaction supported for base_component type")

            if attr in ('lower_bound', 'upper_bound'):
                bounds = pending_bounds.setdefault(id(component),
                                                   [component, None, None])
                bounds[1 if attr == 'lower_bound' else 2] = value
                continue

            if attr == 'metabolites':
                pending_metabolites[id(component)] = (component, feature,
                                                      value)
                continue

            try:
                setattr(component, attr, value)
            except AttributeError as e:
                raise AttributeError(f"Cannot set attribute '{attr}' and no handler is defined for it.") from e

        for component, lower, upper in pending_bounds.values():
            current_lower, current_upper = component.bounds
            component.bounds = (current_lower if lower is None else lower,
                                current_upper if upper is None else upper)

        for component, feature, value in pending_metabolites.values():
            self._set_metabolite_state(component, feature, value)

    @staticmethod
    def _set_metabolite_state(reaction, feature, value):
        """Apply a 'metabolites' feature state, clearing the previous member's.

        ``add_metabolites(..., combine=False)`` overwrites the coefficients it
        is given and leaves every other metabolite alone, so a metabolite that
        one member introduces stays on the reaction when the next member's
        state does not mention it. Visiting members in a different order then
        produces different stoichiometry, and re-visiting a member does not
        reproduce its own first result.

        Every metabolite named by any member's state for this feature is
        therefore written on every visit, with a coefficient of zero where the
        current member does not use it. Keys may be Metabolite objects or ids;
        both are resolved against the reaction's model so the two forms are
        interchangeable.
        """
        model = reaction.model

        def _resolve(key):
            if not isinstance(key, str):
                return key
            if model is not None:
                return model.metabolites.get_by_id(key)
            raise KeyError(
                "Cannot resolve metabolite id %r: the feature's reaction is "
                "not attached to a model." % key)

        union_ids = {}
        for state in feature.states.values():
            for key in state:
                metabolite = _resolve(key)
                union_ids[metabolite.id] = metabolite

        payload = {metabolite: 0.0 for metabolite in union_ids.values()}
        for key, coefficient in value.items():
            payload[_resolve(key)] = coefficient

        reaction.add_metabolites(payload, combine=False)

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
