# Feature backlog

Proposed enhancements that have been intentionally deferred. Each entry
should record what is proposed, the motivation, and why it is not being
done now. Add new entries at the top.

## Standalone feature-construction helper suite

**Source:** Review of PR #115 (Fast BOF) — 2026-05-19

The accepted redesign of PR #115 builds the BOF-style ensemble path
directly into `Ensemble` (extended `__init__` accepting pre-built
features, plus `Ensemble.from_reaction_states` for the
single-reaction/single-attribute case).

**Alternative considered:** keep `Ensemble` minimal and instead start a
separate suite of feature-construction helpers — small functions whose
job is to build `Feature` objects (or lists of them) from common inputs
(a model + reaction id + states dict; a list of models diffed on
bounds; etc.). Callers then pass the resulting features into
`Ensemble(...)`. This keeps `Ensemble`'s surface area small and lets the
helper layer grow independently as new feature-construction recipes are
added.

**Why deferred:** The classmethod approach is enough for the BOF use
case today and matches how `Ensemble` already wraps its
diff-from-models path. Splitting helpers out is worth revisiting once
there are 2+ recipes that would share the helper layer.

## Composable ensemble constructors and shared scaffold

**Source:** Review of PR #115 (Fast BOF) — 2026-05-19

`boundsEnsemble` (varies reaction bounds) and `bofEnsemble` (varies BOF
metabolite coefficients) each return a fresh `Ensemble` built from a
single `cobra.Model`. There is no way today to build a single ensemble
whose members vary in both bounds and BOF coefficients.

Both functions also duplicate the same scaffold: seed an ensemble from
one model, reset `features`/`members` to empty `DictList`s, append
`Feature` objects, then iterate to build `Member` objects. The same
member-construction loop also lives in `Ensemble._populate_members`,
giving three copies overall.

**Proposed direction:**
- Factor a private helper (e.g.
  `_attach_features_to_ensemble(ensemble, features, member_ids, member_names)`)
  shared by `Ensemble._populate_members`, `boundsEnsemble`, and
  `bofEnsemble`.
- Change the two ensemble constructors to accept either a `cobra.Model`
  *or* an existing `Ensemble`, so calls can be chained to layer feature
  types onto the same ensemble.

**Why deferred:** Each constructor is correct on its own. This is
structural cleanup best done alongside the next ensemble-constructor
addition, not in the PR that surfaced it.

**Update — 2026-09-03.** The capability gap is closed:
`Ensemble.from_state_matrix` with a `(reaction_id, component_attribute)`
MultiIndex builds an ensemble varying bounds and biomass coefficients
together. `bofEnsemble` no longer exists, having been replaced by
`Ensemble.from_reaction_states` in #133.

The *structural* half of this entry is now more pressing rather than
less. There are three ways to build a bounds-varying ensemble:
`Ensemble(list_of_models=...)`, `boundsEnsemble`, and
`from_state_matrix`. The last two take essentially the same information
in different shapes — a dict of per-reaction DataFrames versus one
member-by-feature DataFrame — and differ only in that `boundsEnsemble`
is bounds-only while `from_state_matrix` handles any
component_attribute. The open question is whether `from_state_matrix`
should absorb `boundsEnsemble` outright, with `_setBounds` retargeted to
emit a state matrix, rather than the two continuing side by side.
