# MICA Integration Plan for `hubmeta`

## Goal

Add a publishable MICA workflow to `hubmeta` as a Stage-1.5 bridge between
pairwise meta-analytic effect pooling and downstream Stage-2 SEM or path
modeling.

The package story should become:

1. `meta_analysis()` and `morris_weight_analysis()` pool bivariate
   correlations.
2. `mica()` accepts a partial correlation matrix plus optional `k` and `tau`
   matrices.
3. `mica()` returns a positive-definite completed matrix, per-cell
   diagnostics, and a reportability/triage layer for downstream use.

We are **not** bundling real sample databases in v1. If examples are needed,
they should use tiny simulated matrices. If we later want one empirical example,
international entrepreneurship is acceptable; procrastination is out of scope.

## Package Scope for V1

The first package release should expose the mature, reviewable core:

- Input validation and matrix helpers
- Phase 1 diagnostics
- Phase 2 PD-projected deterministic fill
- Triage / reportability rules
- A single front-door `mica()` function
- S3 methods: `print()`, `summary()`, `as.matrix()`

The first package release should **not** expose these as stable public
features yet:

- POET / factor-structured completion
- MNAR sensitivity
- Large benchmark runners
- Stage-2 wrappers for `lavaan` / `metaSEM`
- Real bundled datasets

## Architecture

### Public API

- `mica()`
- `mica_diagnostics()`
- `mica_fill_pd()`
- `mica_triage()`

### Internal layout

- `R/mica.R`
- `R/mica_utils.R`
- `R/mica_diagnostics.R`
- `R/mica_fill.R`
- `R/mica_triage.R`
- `R/mica_methods.R`

### Return object

`mica()` should return an S3 object with:

- `input_matrix`
- `k_matrix`
- `tau_matrix`
- `diagnostics`
- `fill`
- `triage`
- `completed_matrix`
- `method`
- `call`

## Mapping from Current MICA Pipeline

These source files are the right starting point:

- `mica/R/00_utils.R`
- `mica/R/01_diagnostics.R`
- `mica/R/01b_recommender.R`
- `mica/R/02_chordal_baseline.R`

These should inform v2+ package work, but not be required in the first
shipping package surface:

- `mica/R/03b_stan_sampler.R`
- `mica/R/04_factor_poet.R`
- `mica/R/05_mnar_sensitivity.R`
- `mica/R/06_reporting.R`

## Packaging Strategy

### Phase 0: hygiene and baseline

- Fix `DESCRIPTION` metadata
- Remove tracked junk files like `.Rhistory` and `.DS_Store`
- Normalize README install/use instructions
- Get `R CMD check` to a clean baseline

### Phase 1: deterministic MICA core

- Package the utilities, diagnostics, deterministic fill, and triage logic
- Add examples using simulated matrices only
- Add basic tests for:
  - symmetric / diagonal validation
  - deterministic fill returns a PD correlation matrix
  - diagnostics return expected columns
  - triage labels are present for missing cells

### Phase 2: Bayesian core

Development target:

- First integrate the current Stan workflow behind a soft interface so we can
  test packaging and object structure.
- Before CRAN submission, migrate the package-facing implementation to an
  `rstan` / `rstantools` layout if we decide the Stan model must ship in v1.

Important honesty constraints:

- Do not describe Phase 2 deterministic fill as formal GJSW max-det.
- Describe it as a regression-based or PD-projected deterministic fill.
- Do not claim calibrated intervals for all cells.
- Public docs should emphasize reportable cells vs point-only cells.

### Phase 3: release prep

- Vignettes with simulated examples
- GitHub Actions / package checks
- pkgdown site
- CRAN readiness pass

## CRAN-Facing Positioning

Recommended package claim:

`hubmeta` supports Stage-1 meta-analysis and partial correlation-matrix
completion for downstream multivariate modeling. MICA provides a
positive-definite deterministic completion now, with diagnostically guided
reportability and a clear path to Bayesian uncertainty support in later
versions.

That is narrow, honest, and publishable.

## Immediate Build Order

1. Clean package metadata and tracked junk.
2. Add deterministic MICA files and exports.
3. Add lightweight tests.
4. Create a feature branch and push.
5. Decide whether Bayesian MICA lands in the same branch or a follow-up branch.
