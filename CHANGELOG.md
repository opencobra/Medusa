# Changelog

All notable changes to `medusa-cobra` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
