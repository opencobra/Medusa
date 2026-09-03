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

# Number of ids listed before truncating an error message.
_MAX_IDS_IN_ERROR = 10


def _canonicalize_state(value):
    """Return a hashable, order-insensitive stand-in for a feature state.

    Most component_attributes (e.g. 'lower_bound') hold scalars, which are
    already hashable and are returned unchanged. The 'metabolites' attribute
    holds a dict of {metabolite_or_id: coefficient}, which is neither hashable
    nor order-stable, so it is canonicalized to a sorted tuple of
    (str(key), value) pairs. str() is used on the key because cobra Objects
    stringify to their id, so a dict keyed by Metabolite objects and one keyed
    by the equivalent metabolite ids canonicalize identically.

    This exists so that "does this feature vary across members?" can be
    answered the same way for scalar- and dict-valued states.
    """
    if isinstance(value, dict):
        # Sort on the stringified key alone. Sorting whole pairs would fall
        # through to comparing the values whenever two keys tie, which can
        # raise for mixed value types.
        return tuple(sorted(((str(key), value[key]) for key in value),
                            key=lambda item: item[0]))
    return value


def _is_missing_state(value):
    """True if value should be treated as a missing (NaN) cell.

    Guards pandas.isna against container values: a dict-valued state (the
    'metabolites' case) is never missing, and pandas.isna would return an
    array rather than a bool for some containers.
    """
    if isinstance(value, dict):
        return False
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _states_vary(states):
    """True if a feature's states are not identical across all members.

    states : dict of {member_id: value}, i.e. the same structure as
    Feature.states. Values are canonicalized first so that dict-valued states
    ('metabolites') can be compared.
    """
    canonical = [_canonicalize_state(value) for value in states.values()]
    if len(canonical) < 2:
        return False
    try:
        return len(set(canonical)) > 1
    except TypeError:
        # Fall back for values that resist hashing even after
        # canonicalization (e.g. a list-valued state).
        first = canonical[0]
        return any(other != first for other in canonical[1:])


def _find_duplicates(ids):
    """Return the duplicated ids in `ids`, in first-seen order."""
    seen = set()
    duplicates = []
    for identifier in ids:
        if identifier in seen and identifier not in duplicates:
            duplicates.append(identifier)
        seen.add(identifier)
    return duplicates


def _truncated_id_list(ids):
    """Format ids for an error message, capping the number shown."""
    ids = list(ids)
    shown = ids[:_MAX_IDS_IN_ERROR]
    message = ", ".join(str(i) for i in shown)
    if len(ids) > len(shown):
        message += " ... (%i total)" % len(ids)
    return message


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
    def __init__(self, list_of_models=[], identifier=None, name=None,
                 features=None):
        Object.__init__(self,identifier,name)

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
            # An ensemble built from 0 or 1 models has no variable components
            # yet, but features/members must still exist as empty DictLists so
            # that attribute access does not raise AttributeError. Callers that
            # populate them afterwards (e.g. medusa.boundsEnsemble and
            # medusa.reconstruct.expand) simply overwrite them.
            self.features = DictList()
            self.members = DictList()
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
            # Re-resolve to the reaction object that actually lives in
            # base_model so set_state mutates that one.
            feature.base_component = self.base_model.reactions.get_by_id(
                component.id)
            feature.ensemble = self
            self.features.append(feature)

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

    @classmethod
    def from_state_matrix(cls, model, state_matrix,
                          component_attribute='lower_bound',
                          identifier=None, name=None, drop_invariant=True,
                          missing_value=None):
        """Build an ensemble from an explicit member x feature state matrix.

        This is the general-purpose complement to the other two construction
        routes: `Ensemble(list_of_models=[...])` infers features by diffing
        whole models, `Ensemble.from_reaction_states` varies one reaction
        across many members, and this varies many reactions across many
        members without requiring a Model object per member.

        As with `from_reaction_states`, `model` is used directly as the
        ensemble's `base_model` (it is not copied), so `set_state` mutates the
        Model object that was passed in.

        Parameters
        ----------
        model : cobra.Model
            The base model. Every reaction id referenced by `state_matrix`
            must already exist in this model; reactions are never added.
        state_matrix : pandas.DataFrame
            Index is member ids (coerced to str). Columns are either

            (a) reaction ids, in which case every column targets
                `component_attribute`; or
            (b) a 2-level pandas.MultiIndex of
                (reaction_id, component_attribute), which overrides the
                `component_attribute` argument on a per-column basis. This is
                how you build an ensemble in which some reactions vary in
                their lower bound and others in their upper bound.

            Values are the state for that member and feature.
        component_attribute : str, optional
            Reaction attribute targeted by every column, used only when
            `state_matrix.columns` is not a MultiIndex. Defaults to
            'lower_bound'.

            `component_attribute='metabolites'` is supported: states are then
            {metabolite_or_id: coefficient} dicts, because `set_state` routes
            that attribute through `Reaction.add_metabolites(..., combine=False)`
            rather than `setattr`. Invariance testing handles dict-valued cells
            by comparing canonicalized (sorted-tuple) representations.
        identifier, name : str, optional
            Passed through to Ensemble.__init__.
        drop_invariant : bool, optional
            When True (the default), columns whose values do not vary across
            members are skipped, matching the `nunique() > 1` behaviour of
            `_populate_features_base`. With `drop_invariant=False` you get
            Feature objects whose state is the same for every member. That is
            legal — `set_state` will happily set a constant — but it is usually
            a mistake: such a feature adds no ensemble structure, inflates
            `len(ensemble.features)`, and makes every downstream diversity or
            saturation statistic look better than it is.
        missing_value : optional
            Value substituted for NaN cells. When None (the default), the
            fallback is `MISSING_ATTRIBUTE_DEFAULT[component_attribute]` for
            the column's attribute; if that attribute has no entry there (e.g.
            'metabolites'), a NaN cell raises ValueError. The lookup is lazy:
            an attribute with no default is only an error if the matrix
            actually contains a missing cell for it. Note that None therefore
            cannot be used as an explicit fill value.

        Returns
        -------
        Ensemble
        """
        if not isinstance(model, Model):
            raise AttributeError("`model` must be a cobra.core.Model")
        if not isinstance(state_matrix, pd.DataFrame):
            raise AttributeError(
                "`state_matrix` must be a pandas.DataFrame with member ids "
                "as the index and features as the columns."
            )
        if state_matrix.shape[0] == 0:
            raise ValueError(
                "`state_matrix` has no rows; it must contain at least one "
                "member id in its index."
            )
        if state_matrix.shape[1] == 0:
            raise ValueError(
                "`state_matrix` has no columns; it must contain at least one "
                "reaction id in its columns."
            )

        member_ids = [str(member_id) for member_id in state_matrix.index]
        duplicates = _find_duplicates(member_ids)
        if duplicates:
            raise ValueError(
                "`state_matrix` contains duplicate member ids: %s"
                % _truncated_id_list(duplicates)
            )

        # Resolve each column to a (reaction_id, component_attribute) pair.
        if isinstance(state_matrix.columns, pd.MultiIndex):
            if state_matrix.columns.nlevels != 2:
                raise ValueError(
                    "A MultiIndex on `state_matrix.columns` must have exactly "
                    "2 levels, (reaction_id, component_attribute); got %i."
                    % state_matrix.columns.nlevels
                )
            column_targets = [
                (column, (str(column[0]), str(column[1])))
                for column in state_matrix.columns]
        else:
            column_targets = [
                (column, (str(column), component_attribute))
                for column in state_matrix.columns]

        base_reaction_ids = {rxn.id for rxn in model.reactions}
        missing_reactions = []
        seen_missing = set()
        for column, (reaction_id, attribute) in column_targets:
            if reaction_id not in base_reaction_ids and \
                    reaction_id not in seen_missing:
                missing_reactions.append(reaction_id)
                seen_missing.add(reaction_id)
        if missing_reactions:
            raise KeyError(
                "`state_matrix` references reaction ids that are not in "
                "`model`: %s" % _truncated_id_list(missing_reactions)
            )

        duplicate_features = _find_duplicates(
            f"{reaction_id}_{attribute}"
            for column, (reaction_id, attribute) in column_targets)
        if duplicate_features:
            raise ValueError(
                "`state_matrix` columns map to duplicate feature ids: %s. "
                "Each (reaction_id, component_attribute) pair may appear at "
                "most once." % _truncated_id_list(duplicate_features)
            )

        features = []
        for column, (reaction_id, attribute) in column_targets:
            reaction = model.reactions.get_by_id(reaction_id)

            fill = missing_value
            states = {}
            for member_id, value in zip(member_ids, state_matrix[column]):
                if _is_missing_state(value):
                    if fill is None:
                        if attribute not in MISSING_ATTRIBUTE_DEFAULT:
                            raise ValueError(
                                f"`state_matrix` has a missing value for "
                                f"reaction '{reaction_id}', attribute "
                                f"'{attribute}', member '{member_id}', but "
                                f"'{attribute}' has no entry in "
                                "MISSING_ATTRIBUTE_DEFAULT. Pass an explicit "
                                "`missing_value`, or fill the matrix before "
                                "calling from_state_matrix."
                            )
                        fill = MISSING_ATTRIBUTE_DEFAULT[attribute]
                    states[member_id] = fill
                else:
                    states[member_id] = value

            if drop_invariant and not _states_vary(states):
                continue

            features.append(Feature(
                identifier=f"{reaction_id}_{attribute}",
                name=reaction.name,
                base_component=reaction,
                component_attribute=attribute,
                states=states,
            ))

        if not features:
            if drop_invariant and len(member_ids) == 1:
                raise ValueError(
                    "`state_matrix` has only one member, so no column can "
                    "vary across members and every column was dropped. Pass "
                    "drop_invariant=False to build a single-member ensemble."
                )
            if drop_invariant:
                raise ValueError(
                    "No columns of `state_matrix` vary across members, so the "
                    "ensemble would have no features. Check the matrix, or "
                    "pass drop_invariant=False if constant features are "
                    "intended."
                )
            raise ValueError(
                "`state_matrix` produced no features."
            )

        return cls(
            list_of_models=[model],
            features=features,
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

    def _resolve_from_dictlist(self, requested, dictlist, kind):
        """Resolve ids or objects to the objects held by one of our DictLists.

        Accepts a single id/object or an iterable of them. Everything is
        resolved by id, so passing an equivalent object from a copied ensemble
        still returns this ensemble's object rather than the caller's.
        """
        if requested is None:
            return list(dictlist)
        if isinstance(requested, str) or not hasattr(requested, '__iter__'):
            requested = [requested]
        resolved = []
        seen = set()
        for item in requested:
            item_id = item if isinstance(item, str) else item.id
            # A repeated id is treated as a set-style subset rather than as a
            # request for a duplicated row/column, which would produce a
            # DataFrame with non-unique labels.
            if item_id in seen:
                continue
            try:
                resolved.append(dictlist.get_by_id(item_id))
            except KeyError:
                raise KeyError(
                    f"'{item_id}' is not a {kind} of this ensemble."
                ) from None
            seen.add(item_id)
        return resolved

    def feature_state_matrix(self, *, features=None, members=None):
        """Return the ensemble's states as a member x feature DataFrame.

        This is the flat view of the ensemble's structure: one row per member,
        one column per feature, each cell the value that
        `set_state` would assign for that member/feature pair. It is the same
        information held redundantly in `Feature.states` (keyed by member id)
        and `Member.states` (keyed by Feature object), in the orientation most
        analyses want.

        Both axes are sorted. This matters for reproducibility rather than
        cosmetics: `Ensemble._populate_features_base` iterates a Python `set`
        of reaction ids when creating features, so `ensemble.features` comes
        out in an order that is not stable across runs (or across processes,
        given string hash randomization). Two ensembles built from identical
        inputs can therefore hold identically-valued but differently-ordered
        `features` DictLists. Sorting the columns here means callers can
        compare, hash, or concatenate matrices without tripping over that.

        Parameters
        ----------
        features : iterable of str or medusa.core.feature.Feature, optional
            Restrict the columns to these features, given as feature ids or
            Feature objects. A single id/object may be passed directly. When
            None (the default), all features are included.
        members : iterable of str or medusa.core.member.Member, optional
            Restrict the rows to these members, given as member ids or Member
            objects. A single id/object may be passed directly. When None (the
            default), all members are included.

        Returns
        -------
        pandas.DataFrame
            Index is member ids (sorted), columns are feature ids (sorted).
            The dtype follows from the states themselves: numeric for bound
            attributes, object for dict-valued attributes such as
            'metabolites'.
        """
        feature_list = self._resolve_from_dictlist(
            features, self.features, 'feature')
        member_list = self._resolve_from_dictlist(
            members, self.members, 'member')

        member_ids = [member.id for member in member_list]
        columns = {}
        for feature in feature_list:
            columns[feature.id] = pd.Series(
                [feature.states[member_id] for member_id in member_ids],
                index=member_ids, dtype=object)

        matrix = pd.DataFrame(
            columns, index=member_ids,
            columns=[feature.id for feature in feature_list])
        # infer_objects recovers numeric dtypes for scalar-valued attributes
        # while leaving dict-valued ('metabolites') columns as object.
        matrix = matrix.infer_objects()
        return matrix.sort_index(axis=0).sort_index(axis=1)

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
            component = feature.base_component
            attr = feature.component_attribute
            value = feature.states[member.id]

            if not isinstance(component, cobra.core.Reaction):
                raise AttributeError("Only cobra.core.Reaction supported for base_component type")

            try:
                # Try direct assignment first
                setattr(component, attr, value)
            except AttributeError as e:
                # Handle known read-only attributes 
                # TODO only metabolites for now , could add to this
                if attr == "metabolites":
                    component.add_metabolites(value, combine=False)
                else:
                    raise AttributeError(f"Cannot set attribute '{attr}' and no handler is defined for it.") from e

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
