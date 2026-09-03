# Changelog

All notable changes to `medusa-cobra` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `Ensemble.from_state_matrix`, building an ensemble from an explicit
  member-by-feature DataFrame without requiring a `cobra.Model` per member. A
  2-level MultiIndex of `(reaction_id, component_attribute)` on the columns
  lets a single ensemble vary different attributes for different reactions,
  which is the case the "composable ensemble constructors" backlog entry
  identified as impossible: bounds and biomass coefficients varying together.
- `Ensemble.feature_state_matrix`, the inverse view, returning the ensemble's
  states as a member-by-feature DataFrame with both axes sorted so that
  matrices from separately built ensembles can be compared or concatenated.
- `medusa.stats`, a method-agnostic diagnostics module: `variable_features`,
  `n_variable_features`, `diversity_curve`, `prediction_saturation_curve` and
  `consensus_fraction`. These ask how much the members of an ensemble differ
  from one another, structurally or in their predictions, without knowing how
  the ensemble was generated.

### Fixed
- `ensemble_fva` no longer transposes each member's minimum and maximum.
  cobrapy returns the columns ordered `["minimum", "maximum"]` and the result
  was relabelled positionally as `["maximum_<id>", "minimum_<id>"]`, so every
  row labelled as a maximum in fact held that member's minima. On the E. coli
  core model at `fraction_of_optimum=0.1`, PGI's maximum was reported as
  -46.033 against a reported minimum of 9.982; a reported maximum could be
  smaller than its own reported minimum. Columns are now selected by name.
  Present in 0.3.0.
- `Ensemble` built from zero or one model now exposes empty `features` and
  `members` DictLists instead of leaving the attributes unassigned, which made
  ordinary attribute access raise AttributeError. Both shapes are reachable:
  zero models is the documented way to create an empty ensemble, and the
  single-model form is what `boundsEnsemble` and `medusa.reconstruct.expand`
  build before populating features themselves.

## [0.3.0] - 2026-05-20

### Added
- `Ensemble.from_reaction_states` classmethod for building an ensemble
  whose members differ in a single reaction attribute (e.g. alternative
  biomass compositions) (#133).
- `Ensemble.__init__` now accepts a `features=` list of pre-built
  `Feature` objects alongside a single base model (#133).
- `boundsEnsemble` restored alongside `bofEnsemble` as a supported
  ensemble constructor (#114, #133).

### Changed
- `BofDf` rows are now looked up by metabolite id rather than positional
  index, making `bofEnsemble` robust to row reordering (#133).
- `medusa.__version__` is now resolved from installed-package metadata
  via `importlib.metadata`. Previously it was a hardcoded `"0.1"` that
  drifted from `setup.py`'s declared version.

### Fixed
- `test_load_from_file` no longer writes `model[1-4].json` into the repo
  root; tests are routed through pytest's `tmp_path` (#135).

### Infrastructure
- README build-status badge switched from Travis CI to GitHub Actions
  (#134).
- Added `.github/workflows/publish.yml` for tag-triggered PyPI releases
  via trusted publishing (OIDC).
- Removed `medusa_cobra.egg-info/` from version control and added
  `*.egg-info/` to `.gitignore`.
- Set `python_requires='>=3.8'` in `setup.py` to match the runtime
  requirements of `importlib.metadata` and the project's actual Python
  support.

## [0.2.1] - 2020-07-28

Last release before the 2020–2026 hiatus. See git history at tag
`v0.2.1` for details.
