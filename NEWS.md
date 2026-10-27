# TVFundSelect 0.2.0

- `ccy_bh()` and the `CCY-BH` entry in every simulation and empirical workflow
  now use ordinary whole-date residual bootstrap p-values with 9,999 draws,
  centered smoother residual degrees-of-freedom correction, and one-sided BH.
- The fitted null retains average-zero time-varying alpha and fitted loadings.
  Common date indices preserve cross-fund residual dependence. Bootstrap p-values
  retain fitted-null smoothing bias and include the plus-one correction.
- Added public `ccy_bootstrap_bh()` and explicit bootstrap seed/draw controls.
  Parallel reproducibility and independent local-regression refits are tested.
- The four other fund-selection methods and their degrees-of-freedom conventions
  are unchanged. The old analytic CCY calculation remains internal/diagnostic.
- Intermediate results stay in memory. Reproduction writes final summary text
  and optional final figures, never RDS checkpoints or raw replication dumps.
  `resume=TRUE` is no longer supported. `reproduce_paper()` saves one summary.
- Updated the article target mapping and removed the bandwidth study from the
  default paper workflow. The optional diagnostic function remains callable.

# TVFundSelect 0.1.0

- First installable package for the October 24, 2026 manuscript.
- Separate callable implementations of all five fund-screening procedures.
- Fresh Monte Carlo, bandwidth, raw-data rolling and calibrated experiments.
- Portable base-R constancy diagnostics and bootstrap coefficient refits.
- Article target mapping, protected output directories, deterministic parallel
  seeds and documented data/interpretation boundaries.
