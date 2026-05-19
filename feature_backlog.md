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
