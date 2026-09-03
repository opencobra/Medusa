# Changelog

All notable changes to `medusa-cobra` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `medusa.quality.mass_balance.energy_generating_cycles`, which closes every
  boundary reaction and maximizes ATP hydrolysis to detect networks that
  create energy from nothing.
- `medusa.quality.mass_balance.mass_charge_balance`, reporting internal
  reactions that are not elementally and charge balanced, per member.
- `optimize_ensemble` gained an `infeasible` argument (`'warn'`, `'nan'` or
  `'raise'`) controlling how members that fail to solve are reported, and now
  attaches per-member solver statuses to the result as
  `.attrs['member_status']`.
- Test coverage for running two simulations in sequence on one ensemble,
  which previously had none.

### Fixed
- `optimize_ensemble` no longer leaves `ensemble.base_model` holding the
  last-solved member's state. Its two sibling functions already wrapped their
  loop in the model's context manager; this one did not, so every subsequent
  use of the ensemble silently answered a different question, and when the
  final member was infeasible the base model was left unsolvable.
- `optimize_ensemble` no longer reports the solver's stale primal values as
  fluxes for members that did not solve to optimality. Because all members
  share one `base_model`, those values were typically the previously-solved
  member's solution, so an infeasible member silently inherited a neighbour's
  flux distribution. Such members are now NaN.
- `ensemble_single_gene_deletion` no longer defaults `specific_genes` to an
  empty list. cobrapy treats only `None` as "all genes" and an empty list as
  "no genes", so the default call deleted nothing and returned an empty
  result for every member while the docstring promised the opposite.
  `specific_genes` is now required and every empty form raises.
- The deletion functions accept member ids as well as `Member` objects,
  matching `optimize_ensemble`. Previously they required `Member` objects
  while `ensemble_fva` required ids, so no single value worked across all
  three entry points.
- `medusa.quality.mass_balance.leak_test` now runs at all. It previously
  referenced `self` inside a module-level function, called `cobra.Reaction`
  with no import, called `ensemble.optimize_ensemble` which is not a method
  of `Ensemble`, and indexed its results by reaction id where member ids were
  intended. It also accepted an `exchange_prefix` argument and never used it,
  so no boundary reaction was closed and no leak could have been detected.
- Row order of the `optimize_ensemble` result no longer depends on
  `num_processes` or on solver timing.

### Changed
- `optimize_ensemble` emits one aggregated warning naming every member that
  failed, rather than relying on cobrapy's one-warning-per-solve.

## [0.3.0] - 2026-05-20

### Added
- `Ensemble.from_reaction_states` classmethod for building an ensemble
  whose members differ in a single reaction attribute (e.g. alternative
  biomass compositions) (#133).
- `Ensemble.__init__` now accepts a `features=` list of pre-built
  `Feature` objects alongside a single base model (#133).
- `boundsEnsemble` restored as a supported ensemble constructor
  (#114, #133).

### Removed
- `medusa/bofEnsemble.py`, replaced by `Ensemble.from_reaction_states`
  (#133). It was added and removed between 0.2.1 and 0.3.0, so it was
  never part of a published release and no released API is affected.
  (Corrected after the fact: the 0.3.0 entry above originally described
  `bofEnsemble` as a supported constructor, which was true of an
  intermediate commit within #133 but not of the release. The 0.3.0 wheel
  contains `boundsEnsemble.py` and no `bofEnsemble.py`.)

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
