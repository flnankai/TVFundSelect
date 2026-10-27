# Article-level orchestration. Every component computes new data/fits/draws.
paper_settings <- function() {
  list(manuscript = "TV_SS_BH_FDR_bootstrap_CCY_20261027.pdf",
    manuscript_sha256 = "63dff7244f0d9a5f5c98aab8528dfab7d55d61313b56a2deaa37c6c0a03710ee",
    methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"),
    N = 200L, times = c(200L, 400L), examples = 1:3,
    distributions = c("normal", "t3", "mixture", "icm"),
    signals = c(.02, .04, .06, .08, .10, .13, .175, .23, .30),
    replications = 1000L, cores = 25L, gamma = .1, rho = .5,
    spline_degree = 2L, knot_candidates = 0:3,
    ccy_bootstrap_replications = 9999L, ccy_inference = "ordinary residual bootstrap",
    lookbacks = c(120L, 180L), holding = 6L, hac_lag = 6L,
    bootstrap_replications = 1999L,
    calibration_shifts = c(0, .002, .004), calibration_generator = "TV",
    calibration_train_end = 366L, calibration_block_length = 1L,
    factor_parameters = data.frame(factor = 1:3, mu = c(.34,.04,.06),
      phi = c(.05,.07,.04), omega = c(.32,.33,.26), eta = c(.67,.51,.72),
      kappa = c(.13,.03,.05)))
}

paper_targets <- function() {
  data.frame(target = c("Table 1", "Figures 1-6", "Figure 7", "Table 2",
    "Table 3", "Table 4", "Table 5", "Figure 8"),
    function_name = c("paper_settings", "reproduce_simulations", "reproduce_constancy",
      "reproduce_constancy", "reproduce_constancy", "reproduce_calibration",
      "reproduce_empirical", "reproduce_empirical"),
    fresh_computation = c("fixed DGP parameters", "new panels and five methods",
      "new individual regressions", "new individual regressions and null bootstrap",
      "new panel regressions and null bootstrap", "new fitted generators and bootstrap CCY fits",
      "new rolling selection and portfolio evaluation", "new wealth and selection paths"),
    stringsAsFactors = FALSE)
}

# Synthetic raw data are clearly marked; never represented as article inputs.
example_fund_data <- function(N = 30L, TT = 240L, seed = 20261025L) {
  N <- .tf_integer(N, "N", 2L); TT <- .tf_integer(TT, "TT", 60L)
  .tf_with_seed(seed, function() {
    u <- seq_len(TT)/TT
    F <- matrix(stats::rnorm(TT*4L, sd = .02), TT, 4L,
      dimnames = list(NULL, c("Mkt-RF", "SMB", "HML", "MOM")))
    rf <- rep(.002, TT)
    beta <- matrix(c(.8,.2,-.1,.15), 4L, TT) + .05*rep(sin(2*pi*u), each = 4L)
    delta <- seq(-.002, .004, length.out = N)
    returns <- matrix(rep(rf + colSums(t(F)*beta), each = N), N, TT) +
      delta + matrix(stats::rnorm(N*TT, sd = .015), N, TT)
    observed <- matrix(TRUE, N, TT)
    observed[1L, 2L] <- FALSE
    returns[!observed] <- NA_real_
    list(returns = returns, observed = observed, factors = F, rf = rf,
      dates = seq(as.Date("1990-01-01"), by = "month", length.out = TT),
      fund_id = sprintf("synthetic_%04d", seq_len(N)),
      input_md5 = character(), synthetic = TRUE,
      input_convention = list(zero_is_missing = FALSE, rf_in_percent = FALSE,
                              description = "Synthetic API demonstration, not empirical article data"))
  })
}

.paper_export_tables <- function(result, output_dir) {
  tables <- list(Table_1 = paper_settings()$factor_parameters)
  if (!is.null(result$simulations)) tables$Simulation_summary <- result$simulations$summary
  if (!is.null(result$bandwidth)) tables$Optional_bandwidth_diagnostic <- result$bandwidth$summary
  if (!is.null(result$constancy)) {
    tables$Table_2 <- result$constancy$individual_summary
    tables$Table_3 <- result$constancy$panel
  }
  if (!is.null(result$calibration)) tables$Table_4 <- result$calibration$calibration_summary
  if (!is.null(result$empirical)) {
    empirical <- result$empirical
    evaluate <- empirical$strategy_alpha
    tables$Table_5 <- merge(empirical$performance_summary,
      evaluate[evaluate$benchmark == "FF4" & evaluate$hac_lag == 6L, ],
      by = c("L", "method"), all.x = TRUE, sort = FALSE)
  }
  tables
}

reproduce_paper <- function(data = NULL, output_dir = "paper_recomputed",
    components = c("simulations", "constancy", "calibration", "empirical"),
    cores = 25L, simulation_reps = 1000L, calibration_reps = 1000L,
    bootstrap_reps = 1999L, seed = 420261010L, quick = FALSE, ccy_B = 9999L) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  allowed <- c("simulations", "bandwidth", "constancy", "calibration", "empirical")
  if (!length(components) || anyNA(components) || anyDuplicated(components) ||
      any(!components %in% allowed)) stop("Invalid components.")
  for (key in c("cores", "simulation_reps", "calibration_reps", "bootstrap_reps")) {
    .tf_integer(get(key), key)
  }
  .tf_integer(seed, "seed", 0L)
  ccy_B <- .tf_integer(ccy_B, "ccy_B")
  if (!is.logical(quick) || length(quick) != 1L || is.na(quick)) stop("Invalid quick flag.")
  needs_data <- any(components %in% c("constancy", "calibration", "empirical"))
  if (needs_data && is.null(data)) {
    if (quick) data <- example_fund_data(seed = seed)
    else stop("Supply raw empirical data via read_fund_panel(); no archived result data are substituted.")
  }
  if (!is.null(data)) .emp_validate_data(data)
  .tf_output_dir(output_dir)
  if (is.null(output_dir)) stop("reproduce_paper requires a new output_dir.")
  subdir <- function(name) NULL
  effective_ccy_B <- if (quick) min(99L, ccy_B) else ccy_B
  result <- list()
  if ("simulations" %in% components) result$simulations <- if (quick) {
    reproduce_simulations(reps = 2L, cores = cores, N = 40L, examples = 1L,
      times = 80L, distributions = c("normal", "t3"), signals = c(.04,.10),
      seed = seed, output_dir = subdir("simulations"), chunk_size = 1L, ccy_B = effective_ccy_B)
  } else reproduce_simulations(reps = simulation_reps, cores = cores, seed = seed,
                               output_dir = subdir("simulations"), ccy_B = effective_ccy_B)
  if ("bandwidth" %in% components) result$bandwidth <- if (quick) {
    reproduce_bandwidth(reps = 2L, cores = cores, N = 40L, examples = 1L,
      times = 80L, distributions = "normal", signals = c(.04,.10),
      output_dir = subdir("bandwidth"), chunk_size = 1L, ccy_B = effective_ccy_B)
  } else reproduce_bandwidth(reps = simulation_reps, cores = cores,
                              output_dir = subdir("bandwidth"), ccy_B = effective_ccy_B)
  if ("constancy" %in% components) result$constancy <- if (quick) {
    reproduce_constancy(data, B = 9L, cores = cores,
      panel_lengths = min(120L, ncol(data$returns)), common_length = min(120L,ncol(data$returns)),
      fund_rows = seq_len(min(8L,nrow(data$returns))), out_dir = subdir("constancy"))
  } else reproduce_constancy(data, B = bootstrap_reps, cores = cores, out_dir = subdir("constancy"))
  if ("calibration" %in% components) result$calibration <- if (quick) {
    reproduce_calibration(data, replications = 2L, cores = cores, lookbacks = 120L,
      train_end = ncol(data$returns)-2L, shifts = c(0,.002), out_dir = subdir("calibration"), ccy_B = effective_ccy_B)
  } else reproduce_calibration(data, replications = calibration_reps, cores = cores,
                                out_dir = subdir("calibration"), ccy_B = effective_ccy_B)
  if ("empirical" %in% components) result$empirical <- reproduce_empirical(data,
    lookbacks = if (quick) 120L else c(120L,180L), cores = cores,
    out_dir = subdir("empirical"))
  result$tables <- .paper_export_tables(result, output_dir)
  result$manifest <- list(version = as.character(utils::packageVersion("TVFundSelect")),
    manuscript_sha256 = paper_settings()$manuscript_sha256, components = components,
    quick = quick, synthetic_data = !is.null(data$synthetic) && isTRUE(data$synthetic),
    fresh_computation = TRUE, cores = cores, simulation_reps = if (quick) 2L else simulation_reps,
    calibration_reps = if (quick) 2L else calibration_reps,
    bootstrap_reps = if (quick) 9L else bootstrap_reps,
    seed = seed, ccy_B = effective_ccy_B, R_version = R.version.string,
    caveat = "Quick mode is a software demonstration, never article evidence. Monte Carlo estimates vary with seeds.")
  .tf_write_summary(c(list(run_settings = result$manifest), result$tables), output_dir, "paper_summary")
  invisible(result)
}
