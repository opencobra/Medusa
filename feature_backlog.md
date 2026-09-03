# Feature backlog

Proposed enhancements that have been intentionally deferred. Each entry
should record what is proposed, the motivation, and why it is not being
done now. Add new entries at the top.

## Further GENRE quality assessments for medusa.quality

**Source:** Repair of `medusa/quality/mass_balance.py` — 2026-09-03

`medusa.quality` now holds three checks: `leak_test` (can a metabolite be
produced from nothing), `energy_generating_cycles` (can ATP be produced
from nothing), and `mass_charge_balance` (is every internal reaction
elementally and charge balanced). Those cover the defects that make a
reconstruction produce free mass or free energy, which are the ones that
corrupt every downstream flux prediction. Several further assessments are
worth adding.

**Framing.** MEMOTE already implements most single-model GENRE tests well,
and medusa should not reimplement it wholesale. What MEMOTE cannot do is
report per member and separate defects *inherited from the base
reconstruction* from defects *introduced by gapfilling in some members*.
That distinction is the reason to have these in medusa at all: an
ensemble in which 3 of 500 members leak is a gapfilling problem, whereas
one in which all 500 leak is a problem with the input reconstruction. Any
check added here should return a member-indexed result so that split is
readable, and the highest-value future work may simply be a summary layer
that runs a check across members and reports the all/some/none split.

**Tier 1 — structural defects, cheap and high value**

- *Stoichiometric consistency* (Gevorgyan, Poolman & Fell 2008). Tests
  whether a strictly positive molecular-mass vector exists that balances
  every internal reaction. Unlike `leak_test` it is a property of the
  stoichiometric matrix alone, so it detects mass creation that no
  particular flux distribution happens to expose, and unlike
  `mass_charge_balance` it does not depend on metabolite formulas being
  annotated.
- *Blocked reactions.* Reactions that cannot carry flux in any steady
  state, via FVA. Ensemble-aware reporting is the interesting part:
  blocked in every member, versus unblocked only in members that received
  a particular gapfill.
- *Dead-end and orphan metabolites.* Metabolites that are only ever
  produced or only ever consumed. These are the usual cause of a blocked
  reaction and point directly at the missing reaction.

**Tier 2 — growth and biomass sanity**

- *Biomass precursor producibility.* For each biomass component, open a
  demand reaction and test whether it can be produced on the given medium.
  This is the standard first step in debugging a model that will not grow,
  and per member it identifies which gapfill solutions actually closed the
  gap.
- *Biomass consistency.* Whether the biomass reaction's coefficients sum
  to approximately 1 g/gDW, and whether it can carry flux at all. A
  biomass reaction that is off by an order of magnitude silently rescales
  every growth rate the ensemble reports.
- *Unbounded flux.* Any reaction reaching the default ±1000 bound with
  exchanges closed indicates an unconstrained internal cycle.

**Tier 3 — thermodynamics and annotation**

- *Thermodynamically infeasible loops.* Compare FVA against loopless FVA
  and report reactions whose feasible range shrinks. More expensive than
  `energy_generating_cycles` but localizes the loop rather than only
  detecting that one exists.
- *Duplicate reactions.* Distinct ids with identical stoichiometry, which
  inflate apparent network size and create artificial loops.
- *Annotation coverage.* Fraction of reactions carrying gene associations,
  and of metabolites carrying formula, charge, and database
  cross-references. Not a correctness check, but it bounds how far the
  other checks can be trusted: `mass_charge_balance` cannot evaluate a
  reaction whose metabolites have no formulas.

**Tier 4 — checks that only make sense for an ensemble**

- *Invalid feature states.* A member whose state sets `lower_bound` above
  `upper_bound`. This cannot arise in a single model but is easy to
  construct in an ensemble, and `set_state` currently raises partway
  through when it happens, leaving the base model half-updated.
- *Degenerate members.* Members whose feature states are identical to
  another member's. They inflate `len(ensemble.members)` and every
  diversity statistic without adding information.
- *Medium consistency.* Whether every member agrees on which exchange
  reactions are open, so that a growth-rate comparison across members is
  actually a comparison of networks rather than of media.

**Why deferred:** the three checks that landed are the ones that make
downstream flux predictions wrong, and they were the ones needed to make
the module run at all. The rest are additive and each wants its own
design discussion, particularly around whether medusa should wrap MEMOTE
rather than reimplement it.

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
