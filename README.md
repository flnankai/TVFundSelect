# TVFundSelect

An installable R package for **Robust Mutual Fund Selection with Time-Varying
Alphas and Factor Loadings**, manuscript revised October 27, 2026. Version 0.2.0.

**Reproduction means new computation.** The package generates new simulated
panels, refits the selection methods, reconstructs empirical portfolios from
raw fund returns, and recomputes bootstrap and HAC inference. It does not ship
the authors' fitted models, result tables, Monte Carlo checkpoints or frozen
fund selections as substitutes for these calculations.

## Install

```r
install.packages("remotes")
remotes::install_github("flnankai/TVFundSelect")
library(TVFundSelect)
```

R >= 4.1 is required. All computation uses base/recommended R packages; no
compiler, Rcpp, MATLAB or Python is required. Local validation uses R 4.5.2.

## Call each method on your own data

`Y` is a finite **fund-by-time** matrix of excess returns, and `factors` is a
**time-by-factor** matrix without an intercept. Returns and factors must be in
consistent units. Do not mix net returns with excess returns.

```r
panel <- simulate_fund_panel(N = 200, TT = 200, example = 3,
                            distribution = "t3", signal_shift = .10, seed = 100)
Y <- panel$Y
F <- panel$factors

tv_robust <- tv_ss_bh(Y, F)  # proposed time-varying spatial-sign procedure
tv_ls    <- tv_ls_bh(Y, F)  # matched spline least squares, residual-df variance
cc_robust <- cc_ss_bh(Y, F) # constant-coefficient spatial-sign procedure
cc_ls    <- cc_ls_bh(Y, F)  # constant least squares, residual-df variance
local_ls <- ccy_bh(Y, F, B = 9999, seed = 20261027, cores = 25)
# Ordinary residual-bootstrap p-values followed by one-sided BH

all_fits <- fit_fund_methods(Y, F, gamma = .1)
all_fits[["TV-SS-BH"]]$selected  # integer fund row indices
all_fits[["TV-SS-BH"]]$p_values
```

All methods return one-sided p-values and BH selections. The four non-CCY
methods retain their existing formulas and normal reference tails. Both TV
spline procedures share the restricted-residual BIC over 0:3 interior knots
and quadratic splines. TV-SS uses its residual-df multiplier once; TV-LS uses
RSS/(T-q-1), and CC-LS uses RSS/(T-p-1).

**CCY-BH now uses the ordinary residual bootstrap.** Local-linear estimation
retains the normal-reference bandwidth. The fitted null subtracts only the
estimated average alpha, retaining its centered time variation and fitted
factor loadings. The centered residual pool is scaled by sqrt(T/nu), where
nu is the squared Frobenius norm of C(I-S). Each draw resamples whole date
vectors with a shared date index across funds; it uses no random signs and
no within-draw recentering. The p-value is (1 + number of bootstrap average
alphas at least as large as the observed average alpha)/(B + 1). Estimated
smoothing bias is retained. Normal-reference fields are diagnostic only.

`ccy_bh()` always selects this ordinary, df-corrected implementation. The
lower-level `ccy_bootstrap_bh()` additionally exposes explicitly labeled signed
and centered diagnostics. No bootstrap procedure is asserted to guarantee BH
FDR control under arbitrary dependence. With B=9999 the minimum p-value is
0.0001; compare this with gamma/N for your fund universe.

Input validation rejects unidentified designs, nonfinite data, zero residual
radii and nonconvergent spatial equations. It does not silently ridge-fit,
discard failed replications or turn failed methods into empty discoveries.

## One command to recalculate

First verify the complete workflow with fresh synthetic inputs:

```r
reproduce_paper(output_dir = "new_demo", quick = TRUE, cores = 2)
```

Quick mode is a software demonstration, **not evidence reproducing the paper's
empirical findings**. For the full article design, supply the raw empirical
panel and a NEW output directory:

```r
raw <- read_fund_panel("data/fund.csv", "data/factor.csv", "data/rf.csv")
reproduce_paper(data = raw, output_dir = "new_full_run", cores = 25)
```

Defaults are 1,000 new replications per simulation/calibration setting and
1,999 new draws for constancy bootstraps and 9,999 for each CCY bootstrap. Main simulations contain 216 settings
and five fits per panel: over one million new method fits. The removed bandwidth study is not run
by the default article workflow. Full runs are intentionally expensive. Choose smaller
replication counts or individual components for a pilot, not for final claims.

```r
paper_settings()
paper_targets()
reproduce_paper(raw, output_dir = "new_empirical_only",
                components = c("constancy", "calibration", "empirical"),
                calibration_reps = 1000, bootstrap_reps = 1999, cores = 25)
```

Existing output directories are protected. Intermediate arrays and raw
replication metrics stay in memory: no RDS files, checkpoint directories or
raw replication CSVs are written. `reproduce_paper()` writes just
`paper_summary.txt`; separate module functions write their final summary and,
where requested, a final figure. `resume=TRUE` is unsupported because no
checkpoints are stored. Returned R objects still contain metrics and seeds.

The 25-worker outer loops use one inner CCY worker to avoid nested parallelism.
For a single `ccy_bh()` fit, `cores=25` parallelizes bootstrap batches. Explicit
seeds make serial and parallel p-values/selections identical. BLAS threads are
capped at one per worker and caller RNG/thread settings are restored.

## Separate reproduction functions

| Article result | Function | What is recalculated |
|---|---|---|
| Table 1 | `paper_settings()` | Fixed factor-process parameters |
| Figures 1-6 | `reproduce_simulations()` | New DGP panels, all five methods, FDP/TDP and MCSE |
| Tables 2-3, Figure 7 | `reproduce_constancy()` | Raw-history samples, individual/panel fits and null draws |
| Table 4 | `reproduce_calibration()` | New fitted TV means/residuals, independent signed month draws, five fits |
| Table 5, Figure 8 | `reproduce_empirical()` | New rolling selections/holdings, paired and common-model HAC, knot sensitivity |

```r
sim <- reproduce_simulations(reps = 1000, cores = 25,
                             output_dir = "new_simulations")
cc <- reproduce_constancy(raw, B = 1999, cores = 25,
                          out_dir = "new_constancy")
cal <- reproduce_calibration(raw, replications = 1000, cores = 25,
                             out_dir = "new_calibration")
emp <- reproduce_empirical(raw, cores = 25, out_dir = "new_rolling")
```

The full orchestrator computes every requested component anew and collects
the article tables in one summary document. `ccy_B` controls bootstrap draws
in simulation, calibration, rolling-selection and paper entry points. The
optional `reproduce_bandwidth()` is retained only for explicit diagnostics;
it is not part of the current article targets. All functions have `help()`
pages and runnable examples.

## Empirical inputs and economic conventions

Raw CRSP-derived fund data are **not hosted in this repository**. An explicit
function fetches them directly from the original authors' public archive and
uses the verified historical risk-free vintage:

```r
raw <- fetch_empirical_data("data-private")
reproduce_paper(raw, output_dir = "new_full_run", cores = 25)
```

This downloads raw inputs, NOT result tables. It executes no downloaded code.
If R cannot reach the source host, download the two ZIP files in a browser from
the documented original links and process them locally:

```r
raw <- fetch_empirical_data("data-private",
  author_archive = "GPY_code.zip",
  rf_archive = "F-F_Research_Data_Factors_CSV.zip")
```

The real-source offline transformation is numerically verified. R network
connectivity to Dropbox was unavailable in the validation environment; it is
not a requirement for computation after the authorized raw inputs are obtained.
Check your permission to use the external data; the package's code license
does not grant data rights. Offline users can provide the documented three-CSV
format through `read_fund_panel()`. See `inst/doc/data-sources.md` for source links,
historical vintage requirements and reference hashes. `example_fund_data()`
is explicitly synthetic and cannot reproduce empirical results.

The article panel has 726 series over February 1987-September 2017 (368 months).
Fund/factor returns are decimals; the archived RF CSV is percentages. Zero-coded
missing fund returns are identified BEFORE RF subtraction. Row labels are not
CRSP identifiers or consolidated share-class identifiers. Inherited investment
screens cannot be reconstructed from those row labels.

Rolling eligibility uses prior training observations only. Portfolios start
equal-weighted and are bought and held for six months; weights drift. Empty
selections earn RF. A sleeve with missing holding returns earns RF from its
first missing month through block end; this is a valuation convention, not an
observed liquidation return. Final short holding blocks are retained. Alpha
and SE are annualized by multiplying monthly values by 12. Economic evaluation
uses Bartlett HAC with lag 6 and regression-df correction; paired contrasts are
unadjusted, two-sided tests, not one-sided selection tests.

The calibration function defaults to the six TV-generator cases in Table 4:
120/180-month windows ending July 2017, 0/.002/.004 imposed signals, independent
whole-month residual draws and common random signs. Optional CC generators or
longer residual blocks are explicit extensions, not the Table 4 default.
Real observed-fund skill is unknown; no real-data FDP or TDP is reported.

## Reproducibility and uncertainty

All four simulation error families share covariance Sigma; both t3 families
are standardized. Positive signals, including their random baseline, scale
with sqrt(200/T), and alpha paths are centered to preserve their exact average.
The main and bandwidth experiments use independent seed bases. The portable
main task schedule is `seed + 1009*task_id`, not a promise of identical digits
to every historical six-batch run. New Monte Carlo means should be compared
with their MCSEs; different seeds do not imply a different method.

Tests verify independent OLS/kernel/spatial calculations, fixed-design QR
constancy refits, fresh DGP/calibration draws, no-lookahead portfolio handling,
HAC algebra and serial/25-worker equality. Raw empirical reruns and small
paper-dimension simulations are checked separately from full Monte Carlo.
See `inst/doc/verification.md` for the actual validation scope. No claim is
made that packaging re-executed every 1,000-replication experiment.

## References and scope

Methodological antecedents include Benjamini and Hochberg (1995), Ang and
Kristensen (2012), Cheng, Cai and Yan (2018), Ma, Lan, Su and Tsai (2020),
Wang, Zhao, Feng and Wang (2025), and Zhao and Wang (2026). The individual
Ma constancy diagnostic is distinct from the unpublished panel combination
implemented in Appendix B. FDR claims require the manuscript's assumptions;
residual-df correction alone is not a general proof of FDR control.

Code is licensed under GPL-3. External source data retain their own permissions
and access conditions; the code license grants no rights to CRSP-derived data.
