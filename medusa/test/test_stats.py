from cobra.io import load_model

from medusa.core.ensemble import Ensemble
from medusa.stats import (
    variable_features, n_variable_features, diversity_curve,
    prediction_saturation_curve, consensus_fraction)

import numpy as np
import pandas as pd
import pytest

CURVE_COLUMNS = ['size', 'mean', 'std', 'n_draws']


def construct_sparse_ensemble(n_members=20, n_features=40):
    """An ensemble whose features each vary in exactly one member.

    Feature i differs only in member i (i < n_features is not required; the
    varying member is i % n_members). Constructed this way, a subsample of
    size k drawn without replacement contains the odd member for feature i
    with probability k/n_members, so the expected number of variable features
    in a subsample is exactly n_features * k / n_members. That makes the
    diversity curve a known straight line and gives the monotonicity check
    something with real signal in it.
    """
    textbook = load_model("textbook")
    reaction_ids = [rxn.id for rxn in textbook.reactions][:n_features]
    member_ids = ['member_%02i' % i for i in range(n_members)]

    columns = {}
    for i, reaction_id in enumerate(reaction_ids):
        odd_member = i % n_members
        columns[reaction_id] = [
            -999.0 if j == odd_member else -1000.0
            for j in range(n_members)]
    matrix = pd.DataFrame(columns, index=member_ids)
    matrix = matrix[reaction_ids]
    return Ensemble.from_state_matrix(
        textbook, matrix, component_attribute='lower_bound',
        identifier='sparse_ensemble')


def construct_degenerate_ensemble(n_members=20):
    """An ensemble with exactly one variable feature.

    This is the pathological case diversity_curve is meant to expose: 20
    "members" that differ in a single bound.
    """
    textbook = load_model("textbook")
    matrix = pd.DataFrame(
        {'PGI': [-1000.0 + i for i in range(n_members)]},
        index=['member_%02i' % i for i in range(n_members)],
    )
    return Ensemble.from_state_matrix(
        textbook, matrix, identifier='degenerate_ensemble')


def construct_textbook_ensemble():
    model1 = load_model("textbook")
    model1.remove_reactions(model1.reactions[1:3])
    model1.id = 'first_textbook'
    model2 = load_model("textbook")
    model2.remove_reactions(model2.reactions[4:6])
    model2.id = 'second_textbook'
    return Ensemble(list_of_models=[model1, model2],
                    identifier='textbook_ensemble')


def test_variable_features_on_inferred_ensemble():
    ensemble = construct_textbook_ensemble()
    # Every feature inferred by _populate_features_base varies by
    # construction, so all of them must come back.
    variable = variable_features(ensemble)
    assert len(variable) == len(ensemble.features)
    assert n_variable_features(ensemble) == len(ensemble.features)
    assert all(feature in list(ensemble.features) for feature in variable)


def test_variable_features_ignores_constant_features():
    # from_state_matrix(drop_invariant=False) is the supported way to end up
    # with a constant feature, and variable_features must not count it.
    textbook = load_model("textbook")
    matrix = pd.DataFrame(
        {'PGI': [-1000.0, -500.0], 'PFK': [-10.0, -10.0]},
        index=['a', 'b'],
    )
    ensemble = Ensemble.from_state_matrix(
        textbook, matrix, drop_invariant=False)
    assert len(ensemble.features) == 2
    assert n_variable_features(ensemble) == 1
    assert variable_features(ensemble)[0].id == 'PGI_lower_bound'


def test_variable_features_handles_dict_states():
    textbook = load_model("textbook")
    biomass_id = 'Biomass_Ecoli_core'
    biomass = textbook.reactions.get_by_id(biomass_id)
    met = sorted(biomass.metabolites, key=lambda m: m.id)[0]
    matrix = pd.DataFrame(
        {biomass_id: [{met.id: -0.5}, {met.id: -2.0}]},
        index=['a', 'b'],
    )
    ensemble = Ensemble.from_state_matrix(
        textbook, matrix, component_attribute='metabolites')
    assert n_variable_features(ensemble) == 1


def test_diversity_curve_shape_and_reproducibility():
    ensemble = construct_sparse_ensemble()
    curve = diversity_curve(ensemble, n_draws=50, random_state=0)

    assert list(curve.columns) == CURVE_COLUMNS
    assert curve.shape[0] > 1
    assert (curve['n_draws'] == 50).all()
    assert curve['size'].min() == 2
    assert curve['size'].max() == len(ensemble.members)
    assert list(curve['size']) == sorted(curve['size'])

    # Same seed, same curve. Different seed, different curve (the point of
    # taking a seed at all).
    same = diversity_curve(ensemble, n_draws=50, random_state=0)
    pd.testing.assert_frame_equal(curve, same)
    other = diversity_curve(ensemble, n_draws=50, random_state=1)
    assert not curve['mean'].equals(other['mean'])

    # The global numpy RNG must not be touched.
    np.random.seed(1234)
    before = np.random.random()
    np.random.seed(1234)
    diversity_curve(ensemble, n_draws=5, random_state=None)
    assert np.random.random() == before


def test_diversity_curve_is_monotone_and_saturates():
    ensemble = construct_sparse_ensemble(n_members=20, n_features=40)
    total_variable = n_variable_features(ensemble)
    assert total_variable == 40

    curve = diversity_curve(ensemble, n_draws=200, random_state=7,
                            replace=False)

    # Monotone non-decreasing in expectation. Allow a small amount of Monte
    # Carlo slack rather than asserting strict monotonicity of the estimates.
    means = list(curve['mean'])
    for earlier, later in zip(means, means[1:]):
        assert later >= earlier - 1.0

    # It must actually climb, not merely fail to fall.
    assert means[-1] > means[0] + 5.0

    # Without replacement, the largest size is the whole ensemble, so every
    # draw recovers every variable feature exactly.
    largest = curve.iloc[-1]
    assert largest['size'] == len(ensemble.members)
    assert largest['mean'] == pytest.approx(float(total_variable))
    assert largest['std'] == pytest.approx(0.0)

    # And the curve tracks the analytic expectation n_features * k / n.
    n_members = len(ensemble.members)
    for _, row in curve.iterrows():
        expected = total_variable * row['size'] / n_members
        assert row['mean'] == pytest.approx(expected, abs=2.0)


def test_diversity_curve_degenerate_ensemble_is_flat_and_low():
    # The warning case from the docstring: lots of members, almost no
    # structure. The curve must be flat at 1 rather than climbing.
    ensemble = construct_degenerate_ensemble(n_members=20)
    assert n_variable_features(ensemble) == 1

    curve = diversity_curve(ensemble, n_draws=50, random_state=3,
                            replace=False)
    assert curve['mean'].max() <= 1.0
    assert curve['mean'].min() == pytest.approx(1.0)
    assert curve['std'].max() == pytest.approx(0.0)

    # Contrast with the structurally rich ensemble, which does climb.
    rich = diversity_curve(construct_sparse_ensemble(), n_draws=50,
                           random_state=3, replace=False)
    assert rich['mean'].iloc[-1] > curve['mean'].iloc[-1] * 10


def test_diversity_curve_replace_semantics():
    ensemble = construct_sparse_ensemble(n_members=20, n_features=40)

    # With replacement, a size-n draw usually misses some members, so it
    # cannot recover every variable feature.
    with_replacement = diversity_curve(
        ensemble, sizes=[20], n_draws=100, random_state=11, replace=True)
    without = diversity_curve(
        ensemble, sizes=[20], n_draws=100, random_state=11, replace=False)
    assert with_replacement['mean'].iloc[0] < without['mean'].iloc[0]
    assert without['mean'].iloc[0] == pytest.approx(40.0)

    # Sizes above the member count are legal with replacement, an error
    # without it.
    diversity_curve(ensemble, sizes=[50], n_draws=5, random_state=0,
                    replace=True)
    with pytest.raises(ValueError):
        diversity_curve(ensemble, sizes=[50], n_draws=5, random_state=0,
                        replace=False)


def test_diversity_curve_explicit_sizes_and_error_paths():
    ensemble = construct_sparse_ensemble()
    curve = diversity_curve(ensemble, sizes=[2, 5, 10], n_draws=10,
                            random_state=0)
    assert list(curve['size']) == [2, 5, 10]

    with pytest.raises(ValueError):
        diversity_curve(ensemble, sizes=[], n_draws=10)
    with pytest.raises(ValueError):
        diversity_curve(ensemble, sizes=[0, 2], n_draws=10)
    with pytest.raises(ValueError):
        diversity_curve(ensemble, n_draws=0)

    # std is undefined for a single draw.
    single = diversity_curve(ensemble, sizes=[5], n_draws=1, random_state=0)
    assert pd.isna(single['std'].iloc[0])

    # Fewer than 2 members means there is nothing to compare.
    tiny = Ensemble(list_of_models=[load_model("textbook")], identifier='one')
    with pytest.raises(ValueError):
        diversity_curve(tiny)


def test_diversity_curve_no_features_returns_zeros():
    # An ensemble with no features at all is the extreme degenerate case. It
    # should report zero diversity rather than raise, so the caller sees it.
    textbook = load_model("textbook")
    matrix = pd.DataFrame(
        {'PGI': [-1000.0, -500.0]}, index=['a', 'b'])
    ensemble = Ensemble.from_state_matrix(textbook, matrix)
    ensemble.features = type(ensemble.features)()
    curve = diversity_curve(ensemble, n_draws=5, random_state=0)
    assert (curve['mean'] == 0.0).all()


def construct_results_frame(n_members=20, n_disagreeing=3, n_agreeing=5):
    """Binary results with a known number of always-disagreeing columns."""
    member_ids = ['member_%02i' % i for i in range(n_members)]
    columns = {}
    for i in range(n_disagreeing):
        # Alternating, so any subsample of size >= 2 drawn without
        # replacement from consecutive-ish positions will usually disagree;
        # made certain below by splitting the members in half.
        columns['disagree_%i' % i] = [
            1 if j < n_members // 2 else 0 for j in range(n_members)]
    for i in range(n_agreeing):
        columns['agree_%i' % i] = [1] * n_members
    return pd.DataFrame(columns, index=member_ids)


def test_prediction_saturation_curve_known_disagreement():
    results = construct_results_frame(n_members=20, n_disagreeing=3,
                                      n_agreeing=5)
    curve = prediction_saturation_curve(results, n_draws=100, random_state=0,
                                        replace=False)

    assert list(curve.columns) == CURVE_COLUMNS
    assert (curve['n_draws'] == 100).all()
    assert curve['size'].max() == len(results.index)

    # Exactly 3 columns disagree across the full set of members, and no
    # subsample can ever exceed that.
    assert curve['mean'].max() <= 3.0
    largest = curve.iloc[-1]
    assert largest['mean'] == pytest.approx(3.0)
    assert largest['std'] == pytest.approx(0.0)

    # The 5 agreeing columns are never counted, at any size.
    small = prediction_saturation_curve(
        results, sizes=[2], n_draws=200, random_state=0, replace=False)
    assert small['mean'].iloc[0] <= 3.0

    # Reproducible.
    pd.testing.assert_frame_equal(
        curve,
        prediction_saturation_curve(results, n_draws=100, random_state=0,
                                    replace=False))


def test_prediction_saturation_curve_rises_with_sparse_disagreement():
    # One odd member per disagreeing column, so disagreement is only visible
    # in subsamples that happen to include that member.
    n_members = 20
    n_disagreeing = 10
    columns = {}
    for i in range(n_disagreeing):
        columns['disagree_%i' % i] = [
            1 if j == i else 0 for j in range(n_members)]
    results = pd.DataFrame(
        columns, index=['member_%02i' % i for i in range(n_members)])

    curve = prediction_saturation_curve(
        results, n_draws=200, random_state=5, replace=False)
    means = list(curve['mean'])
    for earlier, later in zip(means, means[1:]):
        assert later >= earlier - 1.0
    assert means[-1] > means[0] + 2.0
    assert curve['mean'].iloc[-1] == pytest.approx(float(n_disagreeing))


def test_prediction_saturation_curve_threshold_consensus():
    # Continuous predictions: two columns spread well beyond tolerance, two
    # that differ only at float noise level.
    n_members = 10
    results = pd.DataFrame(
        {
            'wide_0': [float(i) for i in range(n_members)],
            'wide_1': [-float(i) for i in range(n_members)],
            'noise_0': [1.0 + i * 1e-12 for i in range(n_members)],
            'noise_1': [2.0 - i * 1e-12 for i in range(n_members)],
        },
        index=['member_%i' % i for i in range(n_members)],
    )

    # 'exact' flags all four, because the noise columns are not bit-identical.
    exact = prediction_saturation_curve(
        results, sizes=[n_members], n_draws=20, random_state=0, replace=False)
    assert exact['mean'].iloc[0] == pytest.approx(4.0)

    # 'threshold' flags only the two real spreads.
    threshold = prediction_saturation_curve(
        results, sizes=[n_members], n_draws=20, random_state=0,
        consensus='threshold', tolerance=1e-6, replace=False)
    assert threshold['mean'].iloc[0] == pytest.approx(2.0)

    # A tolerance wide enough to swallow everything flags nothing.
    permissive = prediction_saturation_curve(
        results, sizes=[n_members], n_draws=20, random_state=0,
        consensus='threshold', tolerance=1e6, replace=False)
    assert permissive['mean'].iloc[0] == pytest.approx(0.0)


def test_prediction_saturation_curve_error_paths():
    results = construct_results_frame()

    with pytest.raises(AttributeError):
        prediction_saturation_curve({'a': [1, 0]})
    with pytest.raises(ValueError):
        prediction_saturation_curve(results.iloc[0:0])
    with pytest.raises(ValueError):
        prediction_saturation_curve(results, consensus='not_a_mode')

    duplicated = pd.concat([results, results[['agree_0']]], axis=1)
    with pytest.raises(ValueError):
        prediction_saturation_curve(duplicated)

    # threshold consensus needs numbers.
    non_numeric = results.copy()
    non_numeric['label'] = ['up'] * results.shape[0]
    with pytest.raises(ValueError) as excinfo:
        prediction_saturation_curve(non_numeric, consensus='threshold')
    assert 'label' in str(excinfo.value)

    with pytest.raises(ValueError):
        prediction_saturation_curve(results, consensus='threshold',
                                    tolerance=-1.0)


def test_consensus_fraction_arithmetic():
    results = pd.DataFrame(
        {
            'unanimous': [1, 1, 1, 1],
            'three_quarters': [1, 1, 1, 0],
            'even_split': [1, 1, 0, 0],
            'three_way': [1, 1, 2, 3],
        },
        index=['a', 'b', 'c', 'd'],
    )
    fractions = consensus_fraction(results)

    assert isinstance(fractions, pd.Series)
    assert list(fractions.index) == list(results.columns)
    assert fractions['unanimous'] == pytest.approx(1.0)
    assert fractions['three_quarters'] == pytest.approx(0.75)
    assert fractions['even_split'] == pytest.approx(0.5)
    assert fractions['three_way'] == pytest.approx(0.5)


def test_consensus_fraction_booleans_and_nan():
    results = pd.DataFrame(
        {
            'boolean': [True, True, False, True, True],
            'all_nan': [np.nan] * 5,
            'some_nan': [1.0, 1.0, np.nan, np.nan, np.nan],
        },
        index=list('abcde'),
    )
    fractions = consensus_fraction(results)
    assert fractions['boolean'] == pytest.approx(0.8)
    # NaN is a value in its own right here, not dropped.
    assert fractions['all_nan'] == pytest.approx(1.0)
    assert fractions['some_nan'] == pytest.approx(0.6)


def test_consensus_fraction_error_paths():
    with pytest.raises(AttributeError):
        consensus_fraction([1, 0, 1])
    with pytest.raises(ValueError):
        consensus_fraction(pd.DataFrame({'a': []}))
