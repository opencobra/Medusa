
from __future__ import absolute_import

import numpy as np
import pandas as pd

from medusa.core.ensemble import _canonicalize_state
from medusa.core.ensemble import _states_vary

# Structural and prediction diagnostics for an ensemble.
#
# Everything here is method-agnostic: it asks only "how much do the members of
# this ensemble differ from one another?", either structurally (in their
# Feature states) or in their predictions (in a member x prediction results
# table). Nothing here knows or cares how the ensemble was generated.

CURVE_COLUMNS = ['size', 'mean', 'std', 'n_draws']

DEFAULT_CURVE_POINTS = 20


def variable_features(ensemble):
    """Return the Features whose states are not identical across all members.

    A Feature is only informative if at least two members disagree about it.
    Features built by `Ensemble._populate_features_base` and
    `Ensemble.from_state_matrix` (with the default `drop_invariant=True`) are
    variable by construction, but features can be attached directly, and
    `from_state_matrix(..., drop_invariant=False)` deliberately allows constant
    ones, so this cannot be assumed.

    Dict-valued states (a 'metabolites' component_attribute) are compared by
    their canonicalized representation, so key order does not matter and a
    dict keyed by Metabolite objects compares equal to the same dict keyed by
    metabolite ids.

    Parameters
    ----------
    ensemble : medusa.core.Ensemble

    Returns
    -------
    list of medusa.core.feature.Feature
        In `ensemble.features` order, which is not deterministic across runs;
        sort by `.id` if you need a stable order.
    """
    return [feature for feature in ensemble.features
            if _states_vary(feature.states)]


def n_variable_features(ensemble):
    """Number of Features whose states vary across members.

    Parameters
    ----------
    ensemble : medusa.core.Ensemble

    Returns
    -------
    int
    """
    return len(variable_features(ensemble))


def diversity_curve(ensemble, sizes=None, n_draws=100, random_state=None,
                    replace=True):
    """Structural diversity as a function of the number of members sampled.

    For each subsample size, `n_draws` random subsets of members are drawn and
    the number of Features that vary *within that subset* is counted. This is
    the ensemble analogue of a rarefaction curve: the x axis is sampling
    effort (members), the y axis is the structural richness that effort
    recovers (variable features).

    How to read the result
    ----------------------
    - **Rising and then plateauing** is the expected, healthy shape. The
      plateau means additional members stop contributing new structural
      variation, i.e. the ensemble has saturated the structural space it can
      represent *given its inputs*. That is a statement about the ensemble and
      its generating procedure, not about the biology: a plateau says nothing
      about whether the represented space is the right one.
    - **Still climbing steeply at the largest size** means the ensemble is
      undersampled — more members would still add structure.
    - **Flat and low from the very start** is a warning, not a result. It means
      there was very little variation to begin with, so the "ensemble" is close
      to a single model wearing many labels. Any downstream statistic computed
      over it (consensus fractions, prediction spread, uncertainty intervals)
      will look reassuringly tight for the trivial reason that the members
      barely differ. Check `n_variable_features` and the ensemble's
      construction before interpreting anything else.

    Parameters
    ----------
    ensemble : medusa.core.Ensemble
        Must contain at least 2 members. An ensemble with no variable features
        returns an all-zero curve rather than raising — that is the degenerate
        case described above and is worth seeing rather than hiding.
    sizes : iterable of int, optional
        Subsample sizes to evaluate. When None (the default), roughly
        `DEFAULT_CURVE_POINTS` (20) evenly spaced sizes from 2 to
        `len(ensemble.members)` are used, deduplicated after rounding, so small
        ensembles yield fewer than 20 points.
    n_draws : int, optional
        Number of random subsets drawn per size. Defaults to 100.
    random_state : int, array_like, numpy.random.SeedSequence or Generator, optional
        Seed for `numpy.random.default_rng`. Pass an int for a reproducible
        curve. The global numpy RNG is never used or disturbed.
    replace : bool, optional
        Whether members are drawn with replacement within a single subsample.
        Defaults to True (bootstrap-style resampling), which is why sizes may
        exceed `len(ensemble.members)`. Note the consequence: a size-n draw
        from an n-member ensemble will usually miss some members, so the curve
        approaches but does not reach `n_variable_features(ensemble)`. Pass
        `replace=False` for strict rarefaction, in which case every size must
        be <= `len(ensemble.members)` and the largest size reproduces
        `n_variable_features(ensemble)` exactly.

    Returns
    -------
    pandas.DataFrame
        Columns ['size','mean','std','n_draws']: the subsample size, the mean
        and sample standard deviation (ddof=1, hence NaN when `n_draws` is 1)
        of the number of variable features across draws, and the number of
        draws.
    """
    member_ids = [member.id for member in ensemble.members]
    feature_ids = [feature.id for feature in ensemble.features]

    columns = {}
    for feature in ensemble.features:
        columns[feature.id] = pd.Series(
            [_canonicalize_state(feature.states[member_id])
             for member_id in member_ids],
            index=member_ids, dtype=object)
    canonical = pd.DataFrame(columns, index=member_ids, columns=feature_ids)

    def count_variable(drawn):
        if not feature_ids:
            return 0
        subsample = canonical.iloc[drawn]
        return int((subsample.nunique(dropna=False) > 1).sum())

    return _saturation_curve(
        n_units=len(member_ids), count_fn=count_variable, sizes=sizes,
        n_draws=n_draws, random_state=random_state, replace=replace,
        unit_name='member')


def prediction_saturation_curve(results, sizes=None, n_draws=100,
                                random_state=None, consensus='exact',
                                tolerance=1e-6, replace=True):
    """Prediction disagreement as a function of the number of members sampled.

    The prediction-space counterpart to `diversity_curve`. Where
    `diversity_curve` asks how much *structure* a subsample of members
    recovers, this asks how many *predictions* they disagree about. Structural
    diversity that produces no prediction disagreement is diversity the
    downstream analysis cannot see, so the two curves are worth reading
    together.

    Parameters
    ----------
    results : pandas.DataFrame
        Index is member ids, columns are predictions (e.g. reaction ids, gene
        ids, or condition names), values are the per-member prediction. Binary,
        boolean, or continuous. Column labels must be unique.
    sizes : iterable of int, optional
        Subsample sizes to evaluate. When None (the default), roughly
        `DEFAULT_CURVE_POINTS` (20) evenly spaced sizes from 2 to
        `len(results.index)`.
    n_draws : int, optional
        Number of random subsets drawn per size. Defaults to 100.
    random_state : int, array_like, numpy.random.SeedSequence or Generator, optional
        Seed for `numpy.random.default_rng`. The global numpy RNG is never
        used or disturbed.
    consensus : {'exact', 'threshold'}, optional
        How disagreement is decided for a column within a subsample.

        - 'exact' (the default): the column is non-consensus if
          `nunique(dropna=False) > 1`. Appropriate for binary, boolean, or
          categorical predictions. Applied to floats it will flag columns that
          differ only in the last bits, which is usually not what you want.
        - 'threshold': the column is non-consensus if its range
          (`max - min`) exceeds `tolerance`. For continuous predictions such
          as fluxes or growth rates. Requires numeric values.
    tolerance : float, optional
        Range above which a column counts as non-consensus when
        `consensus='threshold'`. Ignored for `consensus='exact'`. Defaults to
        1e-6.
    replace : bool, optional
        Whether members are drawn with replacement within a subsample.
        Defaults to True; see `diversity_curve` for the consequences.

    Returns
    -------
    pandas.DataFrame
        Columns ['size','mean','std','n_draws'], with the same meaning as
        `diversity_curve` except that the counted quantity is non-consensus
        prediction columns rather than variable features. A curve that
        plateaus means additional members stop revealing new disagreement; a
        curve that is flat and low from the start means the members already
        agree almost everywhere, which — as in `diversity_curve` — may reflect
        a degenerate ensemble rather than a confident one.
    """
    _validate_results(results)
    if consensus not in ('exact', 'threshold'):
        raise ValueError(
            "`consensus` must be one of 'exact' or 'threshold'; got "
            "%r." % (consensus,)
        )
    if consensus == 'threshold':
        numeric = results.select_dtypes(include=[np.number])
        if numeric.shape[1] != results.shape[1]:
            non_numeric = [column for column in results.columns
                           if column not in set(numeric.columns)]
            raise ValueError(
                "consensus='threshold' requires numeric values, but these "
                "columns are not numeric: %s"
                % ", ".join(str(column) for column in non_numeric[:10])
            )
        if tolerance < 0:
            raise ValueError("`tolerance` must be non-negative.")

    def count_disagreement(drawn):
        if results.shape[1] == 0:
            return 0
        subsample = results.iloc[drawn]
        if consensus == 'exact':
            return int((subsample.nunique(dropna=False) > 1).sum())
        spread = subsample.max() - subsample.min()
        return int((spread > tolerance).sum())

    return _saturation_curve(
        n_units=results.shape[0], count_fn=count_disagreement, sizes=sizes,
        n_draws=n_draws, random_state=random_state, replace=replace,
        unit_name='member')


def consensus_fraction(results):
    """Fraction of members taking the modal value, per prediction column.

    For each column of `results`, the count of the most common value divided
    by the number of members. 1.0 means unanimity; for a binary prediction the
    minimum is ~0.5 (an even split). This is a raw agreement statistic, not a
    probability and not a calibrated confidence: whether 0.8 agreement means
    an 80% chance of being right is an empirical question that has to be
    settled against held-out observations. Calibrating it is the consumer's
    job, deliberately not done here.

    Parameters
    ----------
    results : pandas.DataFrame
        Index is member ids, columns are predictions. Column labels must be
        unique. NaN is treated as a value in its own right rather than being
        dropped, so a column that is entirely NaN has a consensus fraction of
        1.0.

    Returns
    -------
    pandas.Series
        Indexed by the columns of `results`, in the order they appear there.
    """
    _validate_results(results)
    n_members = results.shape[0]
    fractions = []
    for column in results.columns:
        counts = results[column].value_counts(dropna=False)
        fractions.append(float(counts.max()) / n_members)
    return pd.Series(fractions, index=results.columns, dtype=float)


def _validate_results(results):
    """Shared shape checks for a member x prediction results table."""
    if not isinstance(results, pd.DataFrame):
        raise AttributeError(
            "`results` must be a pandas.DataFrame with member ids as the "
            "index and predictions as the columns."
        )
    if results.shape[0] == 0:
        raise ValueError(
            "`results` has no rows; it must contain at least one member id "
            "in its index."
        )
    if results.columns.has_duplicates:
        raise ValueError(
            "`results` has duplicate column labels; prediction labels must "
            "be unique."
        )


def _default_sizes(n_units, unit_name='member',
                   n_points=DEFAULT_CURVE_POINTS):
    """Roughly n_points evenly spaced sizes from 2 to n_units.

    Rounded to ints and deduplicated, so an ensemble with fewer than
    n_points+1 units simply yields every size from 2 upward.
    """
    if n_units < 2:
        raise ValueError(
            "At least 2 %ss are required to measure variation; got %i."
            % (unit_name, n_units)
        )
    n_points = min(n_points, n_units - 1)
    sizes = np.linspace(2, n_units, num=n_points)
    sizes = np.unique(np.rint(sizes).astype(int))
    return [int(size) for size in sizes]


def _saturation_curve(n_units, count_fn, sizes, n_draws, random_state,
                      replace, unit_name='member'):
    """Shared subsampling loop behind the two curve functions.

    count_fn takes an array of positional indices into the units and returns
    the count of interest for that subsample.
    """
    if sizes is None:
        sizes = _default_sizes(n_units, unit_name=unit_name)
    else:
        sizes = [int(size) for size in sizes]
        if not sizes:
            raise ValueError("`sizes` must not be empty.")
        if min(sizes) < 1:
            raise ValueError(
                "`sizes` must all be >= 1; got %i." % min(sizes)
            )
        if not replace and max(sizes) > n_units:
            raise ValueError(
                "With replace=False, every size must be <= the number of "
                "%ss (%i); got %i." % (unit_name, n_units, max(sizes))
            )
    if n_draws < 1:
        raise ValueError("`n_draws` must be >= 1; got %i." % n_draws)

    rng = np.random.default_rng(random_state)
    positions = np.arange(n_units)

    rows = []
    for size in sizes:
        counts = []
        for _ in range(n_draws):
            drawn = rng.choice(positions, size=size, replace=replace)
            counts.append(count_fn(drawn))
        counts = pd.Series(counts, dtype=float)
        rows.append({
            'size': int(size),
            'mean': counts.mean(),
            'std': counts.std(),
            'n_draws': int(n_draws),
        })
    return pd.DataFrame(rows, columns=CURVE_COLUMNS)
