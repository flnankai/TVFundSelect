# Verification scope (2026-10-25)

The checks below are implementation checks, not a claim that every full Monte
Carlo experiment was re-executed while building the package. Reproduction APIs
are tested through fresh calculations, not by substituting saved result files.

## Core and portable execution

- All five method scores, p-values and selections agree with the source
  algorithms, with independent OLS, local weighted-regression, BIC and joint
  spatial-equation checks.
- Constancy statistics and bootstrap draws agree with independent full
  restricted/unrestricted regression refits; residual correlations are updated
  in every draw. The base-R implementation needs no native compiler.
- Installed-package tests verify serial versus 25 actual local PSOCK workers
  for simulations, bandwidth comparisons, rolling evaluation, data calibration
  and constancy bootstrap. Task seeds/draws and results agree; caller RNG and
  numerical-library thread settings are restored.
- The complete small article workflow regenerates panels, synthetic raw input,
  fits and bootstrap draws, then writes new Table 1-9 CSVs and plots. Output
  overwrite refusal and explicit same-run simulation checkpoint resumption are
  tested. No original author result is loaded by the workflows.

## Fresh empirical recomputation

The authorized local original raw 726-by-368 return matrix was used to freshly
fit all 74 formation windows for both TV methods and all four fixed-knot
alternatives. Only AFTER recomputation were the reference calculations read for
verification. The Table 6/8 alpha and HAC values and all sixteen Table 9
selection/participation/overlap rows agree within 1e-10. Table 7 paired HAC
values also agree, using the single unadjusted two-sided contrast.

Both original constancy samples and all three panels were recalculated with
1999 NEW bootstrap draws, using 25 workers. The longest-history sample contains
725 estimable funds, with one unavailable 13-month history; the common sample
contains 499. Individual normal/bootstrap p-values and all panel hybrid
p-values agree with the reference after recomputation. In particular, the raw
5%/10%/BH10 rejection counts are respectively 695/698/697 (longest normal),
695/697/697 (longest bootstrap), 460/469/469 (common normal), and 458/470/469
(common bootstrap). These values document the check; they are not inputs to R.

## Fresh Monte Carlo smoke checks

- Eight NEW paper-dimension panels (N=200, T=200/400, all four errors) produced
  forty successful method fits. All twelve example/distribution combinations
  also regenerate identical DGP panels for matching source seeds.
- Six empirical TV calibration cases, two NEW full-size N=504 replications
  each, produced sixty method fits. Fitted generators were reconstructed from
  raw data. Seeds, selections, FDP and TDP agree with the matching reference
  replications only checked after computation.
- No new full 216-setting, 1000-replication suite was claimed as validation.
  The package functions default to those full settings for users to rerun.
  A new portable main seed schedule may give different Monte Carlo digits;
  uncertainty is recorded with MCSEs, not hidden by recycling old summaries.

## Source acquisition

The original public author MAT archive and historical French RF ZIP were
downloaded from their primary links and processed offline with the new base-R
reader. The last 368 fund/factor observations and all RF entries exactly match
the article inputs, verified using canonical numeric hashes. Compressed/plain
MAT fixtures and failure cases are tested. Direct live R downloads from
Dropbox timed out in this environment even with external network permission;
the documented browser-download/offline-ZIP path is fully verified. No private
fund observations, original result tables or unpublished PDF are hosted here.

## Package checks

`R CMD build` and `R CMD check --no-manual` are used for source installation,
namespace/dependency consistency, static R analysis, Rd documentation, examples
and all tests. GitHub Actions repeats checks on Windows and Linux. Current
remote run status is displayed by GitHub; do not infer success from dispatch.
