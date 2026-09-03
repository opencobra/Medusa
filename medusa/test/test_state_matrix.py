from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.core.ensemble import MISSING_ATTRIBUTE_DEFAULT

import pandas as pd
import pytest

REACTION_ATTRIBUTES = ['lower_bound', 'upper_bound']


def construct_textbook_ensemble():
    # Two textbook models with different reactions removed, so that features
    # are inferred for both bounds of the removed reactions.
    model1 = load_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = load_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    return Ensemble(list_of_models=[model1, model2],
                    identifier='textbook_ensemble')


def construct_state_matrix(n_members=4):
    # A member x reaction matrix of lower bounds for three real textbook
    # reactions. PGI varies, PFK varies, and PGK is deliberately constant so
    # drop_invariant has something to drop.
    member_ids = ['member_%i' % i for i in range(n_members)]
    return pd.DataFrame(
        {
            'PGI': [-1000.0 * (i + 1) / n_members for i in range(n_members)],
            'PFK': [0.0 if i % 2 else -10.0 for i in range(n_members)],
            'PGK': [-1000.0] * n_members,
        },
        index=member_ids,
    )


def test_feature_state_matrix_shape_and_content():
    ensemble = construct_textbook_ensemble()
    matrix = ensemble.feature_state_matrix()

    assert isinstance(matrix, pd.DataFrame)
    assert matrix.shape == (len(ensemble.members), len(ensemble.features))
    assert set(matrix.index) == {member.id for member in ensemble.members}
    assert set(matrix.columns) == {feature.id for feature in ensemble.features}

    # Every cell must equal the corresponding Feature.states entry.
    for feature in ensemble.features:
        for member in ensemble.members:
            assert matrix.loc[member.id, feature.id] == \
                feature.states[member.id]


def test_feature_state_matrix_is_sorted_and_deterministic():
    ensemble = construct_textbook_ensemble()
    matrix = ensemble.feature_state_matrix()

    assert list(matrix.index) == sorted(matrix.index)
    assert list(matrix.columns) == sorted(matrix.columns)

    # Two calls must agree, and the sorted axes must not depend on the
    # (nondeterministic) order of ensemble.features.
    again = ensemble.feature_state_matrix()
    pd.testing.assert_frame_equal(matrix, again)

    reordered = ensemble.feature_state_matrix(
        features=[feature.id for feature in
                  sorted(ensemble.features, key=lambda f: f.id,
                         reverse=True)],
        members=[member.id for member in
                 sorted(ensemble.members, key=lambda m: m.id, reverse=True)],
    )
    pd.testing.assert_frame_equal(matrix, reordered)


def test_feature_state_matrix_subsets():
    ensemble = construct_textbook_ensemble()
    feature = ensemble.features[0]
    member = ensemble.members[0]

    # Subset by object.
    subset = ensemble.feature_state_matrix(
        features=[feature], members=[member])
    assert subset.shape == (1, 1)
    assert subset.loc[member.id, feature.id] == feature.states[member.id]

    # Subset by id, passed as a bare string rather than a list.
    single = ensemble.feature_state_matrix(
        features=feature.id, members=member.id)
    pd.testing.assert_frame_equal(subset, single)

    # Unknown ids are an error on both axes.
    with pytest.raises(KeyError):
        ensemble.feature_state_matrix(features=['not_a_feature'])
    with pytest.raises(KeyError):
        ensemble.feature_state_matrix(members=['not_a_member'])

    # A repeated id is a set-style subset, not a request for a duplicated
    # label, which would yield a DataFrame with non-unique columns.
    repeated = ensemble.feature_state_matrix(
        features=[feature.id, feature.id], members=[member.id, member.id])
    pd.testing.assert_frame_equal(subset, repeated)
    assert not repeated.columns.has_duplicates
    assert not repeated.index.has_duplicates


def test_feature_state_matrix_roundtrips_through_from_state_matrix():
    # feature_state_matrix and from_state_matrix should be inverses when the
    # columns carry their attributes as a MultiIndex.
    ensemble = construct_textbook_ensemble()
    matrix = ensemble.feature_state_matrix()

    # Rebuild the (reaction_id, attribute) pairs from the features so the
    # matrix can be fed back in.
    pairs = {}
    for feature in ensemble.features:
        pairs[feature.id] = (feature.base_component.id,
                             feature.component_attribute)
    multi = matrix.copy()
    multi.columns = pd.MultiIndex.from_tuples(
        [pairs[column] for column in matrix.columns])

    rebuilt = Ensemble.from_state_matrix(
        load_model("textbook"), multi, identifier='rebuilt')

    assert len(rebuilt.features) == len(ensemble.features)
    assert {member.id for member in rebuilt.members} == \
        {member.id for member in ensemble.members}
    pd.testing.assert_frame_equal(
        rebuilt.feature_state_matrix(), matrix, check_dtype=False)


def test_from_state_matrix_single_attribute():
    textbook = load_model("textbook")
    matrix = construct_state_matrix()

    ensemble = Ensemble.from_state_matrix(
        textbook, matrix, component_attribute='lower_bound',
        identifier='from_matrix')

    # PGK is invariant and should have been dropped by default.
    assert {feature.id for feature in ensemble.features} == \
        {'PGI_lower_bound', 'PFK_lower_bound'}
    assert {member.id for member in ensemble.members} == set(matrix.index)

    for feature in ensemble.features:
        assert feature.component_attribute == 'lower_bound'
        assert feature.ensemble is ensemble
        # base_component must be re-resolved against the ensemble's base model.
        assert feature.base_component in ensemble.base_model.reactions
        assert feature.base_component is \
            ensemble.base_model.reactions.get_by_id(
                feature.base_component.id)

    for member in ensemble.members:
        assert member.ensemble is ensemble
        assert len(member.states) == len(ensemble.features)

    # set_state must push the matrix value onto the base model reaction.
    pgi = ensemble.base_model.reactions.get_by_id('PGI')
    for member_id in matrix.index:
        ensemble.set_state(member_id)
        assert pgi.lower_bound == pytest.approx(matrix.loc[member_id, 'PGI'])


def test_from_state_matrix_multiindex_columns():
    textbook = load_model("textbook")
    matrix = pd.DataFrame(
        {
            ('PGI', 'lower_bound'): [-1000.0, -500.0],
            ('PGI', 'upper_bound'): [1000.0, 250.0],
            ('PFK', 'upper_bound'): [1000.0, 0.0],
        },
        index=['a', 'b'],
    )
    matrix.columns = pd.MultiIndex.from_tuples(matrix.columns)

    # component_attribute is overridden per column, so passing a bogus one
    # must have no effect.
    ensemble = Ensemble.from_state_matrix(
        textbook, matrix, component_attribute='ignored_because_multiindex')

    assert {feature.id for feature in ensemble.features} == \
        {'PGI_lower_bound', 'PGI_upper_bound', 'PFK_upper_bound'}
    attributes = {feature.id: feature.component_attribute
                  for feature in ensemble.features}
    assert attributes['PGI_lower_bound'] == 'lower_bound'
    assert attributes['PGI_upper_bound'] == 'upper_bound'
    assert attributes['PFK_upper_bound'] == 'upper_bound'

    ensemble.set_state('b')
    pgi = ensemble.base_model.reactions.get_by_id('PGI')
    assert pgi.lower_bound == pytest.approx(-500.0)
    assert pgi.upper_bound == pytest.approx(250.0)
    assert ensemble.base_model.reactions.get_by_id('PFK').upper_bound == \
        pytest.approx(0.0)


def test_from_state_matrix_multiindex_rejects_wrong_nlevels():
    textbook = load_model("textbook")
    matrix = pd.DataFrame({('PGI', 'lower_bound', 'extra'): [-1.0, -2.0]},
                          index=['a', 'b'])
    matrix.columns = pd.MultiIndex.from_tuples(matrix.columns)
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(textbook, matrix)


def test_from_state_matrix_drop_invariant_both_ways():
    textbook = load_model("textbook")
    matrix = construct_state_matrix()

    dropped = Ensemble.from_state_matrix(textbook, matrix,
                                         drop_invariant=True)
    assert len(dropped.features) == 2
    assert 'PGK_lower_bound' not in \
        {feature.id for feature in dropped.features}

    kept = Ensemble.from_state_matrix(load_model("textbook"), matrix,
                                      drop_invariant=False)
    assert len(kept.features) == 3
    constant = kept.features.get_by_id('PGK_lower_bound')
    assert len(set(constant.states.values())) == 1

    # A matrix in which nothing varies is an error with the default, and
    # constructible without it.
    invariant = matrix[['PGK']]
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(load_model("textbook"), invariant)
    forced = Ensemble.from_state_matrix(
        load_model("textbook"), invariant, drop_invariant=False)
    assert len(forced.features) == 1


def test_from_state_matrix_metabolites_attribute():
    textbook = load_model("textbook")
    biomass_id = 'Biomass_Ecoli_core'
    biomass = textbook.reactions.get_by_id(biomass_id)
    target_met = sorted(biomass.metabolites, key=lambda m: m.id)[0]

    matrix = pd.DataFrame(
        {biomass_id: [{target_met.id: -0.5}, {target_met.id: -2.0}]},
        index=['alt_a', 'alt_b'],
    )
    ensemble = Ensemble.from_state_matrix(
        textbook, matrix, component_attribute='metabolites')

    assert len(ensemble.features) == 1
    feature = ensemble.features[0]
    assert feature.id == f"{biomass_id}_metabolites"
    assert feature.component_attribute == 'metabolites'

    base_biomass = ensemble.base_model.reactions.get_by_id(biomass_id)
    ensemble.set_state('alt_a')
    assert base_biomass.metabolites[target_met] == pytest.approx(-0.5)
    ensemble.set_state('alt_b')
    assert base_biomass.metabolites[target_met] == pytest.approx(-2.0)

    # The state matrix accessor must survive dict-valued cells.
    states = ensemble.feature_state_matrix()
    assert states.shape == (2, 1)
    assert states.loc['alt_a', feature.id] == {target_met.id: -0.5}


def test_from_state_matrix_drop_invariant_handles_dict_cells():
    # Dict-valued states must be compared by content, not identity, and key
    # order must not make two equal dicts look different.
    textbook = load_model("textbook")
    biomass_id = 'Biomass_Ecoli_core'
    biomass = textbook.reactions.get_by_id(biomass_id)
    mets = sorted(biomass.metabolites, key=lambda m: m.id)[:2]

    same = pd.DataFrame(
        {biomass_id: [
            {mets[0].id: -0.5, mets[1].id: -1.5},
            {mets[1].id: -1.5, mets[0].id: -0.5},
        ]},
        index=['a', 'b'],
    )
    with pytest.raises(ValueError):
        # Identical content in a different key order is invariant.
        Ensemble.from_state_matrix(
            textbook, same, component_attribute='metabolites')

    different = pd.DataFrame(
        {biomass_id: [
            {mets[0].id: -0.5, mets[1].id: -1.5},
            {mets[0].id: -0.5, mets[1].id: -9.0},
        ]},
        index=['a', 'b'],
    )
    ensemble = Ensemble.from_state_matrix(
        load_model("textbook"), different,
        component_attribute='metabolites')
    assert len(ensemble.features) == 1


def test_from_state_matrix_missing_values():
    textbook = load_model("textbook")
    matrix = pd.DataFrame(
        {'PGI': [-1000.0, float('nan')], 'PFK': [-10.0, 0.0]},
        index=['a', 'b'],
    )

    # Default fill comes from MISSING_ATTRIBUTE_DEFAULT.
    default_fill = Ensemble.from_state_matrix(textbook, matrix)
    feature = default_fill.features.get_by_id('PGI_lower_bound')
    assert feature.states['b'] == MISSING_ATTRIBUTE_DEFAULT['lower_bound']

    # Explicit missing_value wins.
    explicit = Ensemble.from_state_matrix(
        load_model("textbook"), matrix, missing_value=-7.5)
    assert explicit.features.get_by_id('PGI_lower_bound').states['b'] == -7.5


def test_from_state_matrix_missing_value_without_default_raises():
    # 'metabolites' has no MISSING_ATTRIBUTE_DEFAULT entry, so a NaN cell has
    # no defensible fill and must raise.
    textbook = load_model("textbook")
    biomass_id = 'Biomass_Ecoli_core'
    matrix = pd.DataFrame(
        {biomass_id: [{'atp_c': -1.0}, float('nan')]},
        index=['a', 'b'],
    )
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(
            textbook, matrix, component_attribute='metabolites')


def test_from_state_matrix_error_paths():
    textbook = load_model("textbook")
    matrix = construct_state_matrix()

    # Not a Model.
    with pytest.raises(AttributeError):
        Ensemble.from_state_matrix('not_a_model', matrix)

    # Not a DataFrame.
    with pytest.raises(AttributeError):
        Ensemble.from_state_matrix(textbook, {'PGI': [-1.0, -2.0]})

    # Empty on either axis.
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(textbook, matrix.iloc[0:0])
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(textbook, matrix.iloc[:, 0:0])

    # Duplicate member ids.
    duplicated = pd.concat([matrix, matrix.iloc[[0]]])
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(textbook, duplicated)

    # Reaction id not in the model. Must be a KeyError naming the offenders.
    unknown = matrix.rename(columns={'PGI': 'NOT_A_REACTION',
                                     'PFK': 'ALSO_MISSING'})
    with pytest.raises(KeyError) as excinfo:
        Ensemble.from_state_matrix(textbook, unknown)
    assert 'NOT_A_REACTION' in str(excinfo.value)
    assert 'ALSO_MISSING' in str(excinfo.value)

    # Duplicate (reaction, attribute) pairs.
    duplicate_pairs = pd.DataFrame(
        {('PGI', 'lower_bound'): [-1.0, -2.0],
         ('PGI', 'lower_bound '): [-3.0, -4.0]},
        index=['a', 'b'],
    )
    duplicate_pairs.columns = pd.MultiIndex.from_tuples(
        [('PGI', 'lower_bound'), ('PGI', 'lower_bound')])
    with pytest.raises(ValueError):
        Ensemble.from_state_matrix(textbook, duplicate_pairs)


def test_from_state_matrix_single_member():
    # One member cannot vary, so the default drops everything and says why.
    textbook = load_model("textbook")
    matrix = pd.DataFrame({'PGI': [-1000.0]}, index=['solo'])

    with pytest.raises(ValueError) as excinfo:
        Ensemble.from_state_matrix(textbook, matrix)
    assert 'one member' in str(excinfo.value)

    single = Ensemble.from_state_matrix(
        load_model("textbook"), matrix, drop_invariant=False)
    assert len(single.members) == 1
    assert len(single.features) == 1
    single.set_state('solo')
    assert single.base_model.reactions.get_by_id('PGI').lower_bound == \
        pytest.approx(-1000.0)


def test_from_state_matrix_coerces_member_ids_to_str():
    textbook = load_model("textbook")
    matrix = pd.DataFrame({'PGI': [-1000.0, -500.0]}, index=[0, 1])
    ensemble = Ensemble.from_state_matrix(textbook, matrix)
    assert {member.id for member in ensemble.members} == {'0', '1'}
    assert set(ensemble.features[0].states.keys()) == {'0', '1'}


def test_empty_ensemble_has_empty_features_and_members():
    # Regression: Ensemble.__init__ used to leave features/members unassigned
    # when list_of_models had 0 or 1 entries, so attribute access raised
    # AttributeError.
    empty = Ensemble(identifier='empty')
    assert len(empty.features) == 0
    assert len(empty.members) == 0

    single = Ensemble(list_of_models=[load_model("textbook")],
                      identifier='single')
    assert len(single.features) == 0
    assert len(single.members) == 0

    # And the basic accessor works on them rather than blowing up.
    matrix = single.feature_state_matrix()
    assert matrix.shape == (0, 0)
