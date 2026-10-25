# Article-level orchestration. Every component computes new data/fits/draws.
paper_settings <- function() {
  list(manuscript = "TV_SS_BH_FDR_appendix_revised_20261024.pdf",
    manuscript_sha256 = "FDB45238AAB70CBAA08D4CF7486987F540596FA9CA9804E23098352A4B4171FB",
    methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"),
    N = 200L, times = c(200L, 400L), examples = 1:3,
    distributions = c("normal", "t3", "mixture", "icm"),
    signals = c(.02, .04, .06, .08, .10, .13, .175, .23, .30),
    replications = 1000L, cores = 25L, gamma = .1, rho = .5,
    spline_degree = 2L, knot_candidates = 0:3,
    bandwidth_multipliers = c(1, 2, 4, 8, 16),
    lookbacks = c(120L, 180L), holding = 6L, hac_lag = 6L,
    bootstrap_replications = 1999L,
    calibration_shifts = c(0, .002, .004), calibration_generator = "TV",
    calibration_train_end = 366L, calibration_block_length = 1L,
    factor_parameters = data.frame(factor = 1:3, mu = c(.34,.04,.06),
      phi = c(.05,.07,.04), omega = c(.32,.33,.26), eta = c(.67,.51,.72),
      kappa = c(.13,.03,.05)))
}

paper_targets <- function() {
  data.frame(target = c("Table 1", "Table 2", "Figures 1-6", "Figure 7",
    "Table 3", "Table 4", "Table 5", "Table 6", "Table 7", "Table 8",
    "Table 9", "Figure 8", "Figures D.1-D.6"),
    function_name = c("paper_settings", "reproduce_bandwidth", "reproduce_simulations",
      "reproduce_constancy", "reproduce_constancy", "reproduce_constancy",
      "reproduce_calibration", rep("reproduce_empirical",5), "reproduce_bandwidth"),
    fresh_computation = c("fixed DGP parameters", "new panels and all bandwidth fits",
      "new panels and five methods", "new individual regressions",
      "new individual regressions and null bootstrap", "new panel regressions and null bootstrap",
      "new fitted means/residuals, residual draws and all five methods",
      "new rolling selections, holdings and FF4 HAC evaluation",
      "new paired return-difference HAC regression", "new common TV-FF4 HAC evaluation",
      "new rolling fits at every fixed knot count", "new wealth and selection paths",
      "new panels and all bandwidth fits"), stringsAsFactors = FALSE)
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
  save_table <- function(tab, name) if (!is.null(tab) && is.data.frame(tab)) {
    utils::write.csv(tab, file.path(output_dir, paste0(name,".csv")), row.names = FALSE)
  }
  save_table(paper_settings()$factor_parameters, "Table_1")
  if (!is.null(result$bandwidth)) {
    x <- result$bandwidth$summary
    table2 <- do.call(rbind, lapply(split(x, x$bandwidth_multiplier), function(d) {
      data.frame(c = d$bandwidth_multiplier[1L],
        h_T200 = 1.06*stats::sd((1:200)/200)*200^(-1/5)*d$bandwidth_multiplier[1L],
        h_T400 = 1.06*stats::sd((1:400)/400)*400^(-1/5)*d$bandwidth_multiplier[1L],
        mean_FDP_min = min(d$mean_FDP), mean_FDP_max = max(d$mean_FDP),
        settings_above_nominal = sum(d$mean_FDP > result$bandwidth$config$gamma), settings = nrow(d))
    }))
    save_table(table2, "Table_2")
  }
  if (!is.null(result$constancy)) {
    save_table(result$constancy$individual_summary, "Table_3")
    save_table(result$constancy$panel, "Table_4")
  }
  if (!is.null(result$calibration)) save_table(result$calibration$calibration_summary, "Table_5")
  if (!is.null(result$empirical)) {
    empirical <- result$empirical
    evaluate <- empirical$strategy_alpha
    common <- merge(empirical$performance_summary,
      evaluate[evaluate$benchmark == "FF4" & evaluate$hac_lag == 6L, ],
      by = c("L", "method"), all.x = TRUE, sort = FALSE)
    save_table(common, "Table_6")
    save_table(empirical$paired_alpha[empirical$paired_alpha$benchmark == "FF4" &
      empirical$paired_alpha$hac_lag == 6L, ], "Table_7")
    save_table(evaluate[evaluate$benchmark == "TV-FF4" & evaluate$hac_lag == 6L, ], "Table_8")
    save_table(empirical$sensitivity_summary, "Table_9")
  }
}

reproduce_paper <- function(data = NULL, output_dir = "paper_recomputed",
    components = c("simulations", "bandwidth", "constancy", "calibration", "empirical"),
    cores = 25L, simulation_reps = 1000L, calibration_reps = 1000L,
    bootstrap_reps = 1999L, seed = 420261010L, quick = FALSE) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  allowed <- c("simulations", "bandwidth", "constancy", "calibration", "empirical")
  if (!length(components) || anyNA(components) || anyDuplicated(components) ||
      any(!components %in% allowed)) stop("Invalid components.")
  for (key in c("cores", "simulation_reps", "calibration_reps", "bootstrap_reps")) {
    .tf_integer(get(key), key)
  }
  .tf_integer(seed, "seed", 0L)
  if (!is.logical(quick) || length(quick) != 1L || is.na(quick)) stop("Invalid quick flag.")
  needs_data <- any(components %in% c("constancy", "calibration", "empirical"))
  if (needs_data && is.null(data)) {
    if (quick) data <- example_fund_data(seed = seed)
    else stop("Supply raw empirical data via read_fund_panel(); no archived result data are substituted.")
  }
  if (!is.null(data)) .emp_validate_data(data)
  .tf_output_dir(output_dir)
  if (is.null(output_dir)) stop("reproduce_paper requires a new output_dir.")
  utils::write.csv(paper_targets(), file.path(output_dir, "targets.csv"), row.names = FALSE)
  subdir <- function(name) file.path(output_dir, name)
  result <- list()
  if ("simulations" %in% components) result$simulations <- if (quick) {
    reproduce_simulations(reps = 2L, cores = cores, N = 40L, examples = 1L,
      times = 80L, distributions = c("normal", "t3"), signals = c(.04,.10),
      seed = seed, output_dir = subdir("simulations"), chunk_size = 1L)
  } else reproduce_simulations(reps = simulation_reps, cores = cores, seed = seed,
                               output_dir = subdir("simulations"))
  if ("bandwidth" %in% components) result$bandwidth <- if (quick) {
    reproduce_bandwidth(reps = 2L, cores = cores, N = 40L, examples = 1L,
      times = 80L, distributions = "normal", signals = c(.04,.10),
      output_dir = subdir("bandwidth"), chunk_size = 1L)
  } else reproduce_bandwidth(reps = simulation_reps, cores = cores,
                              output_dir = subdir("bandwidth"))
  if ("constancy" %in% components) result$constancy <- if (quick) {
    reproduce_constancy(data, B = 9L, cores = cores,
      panel_lengths = min(120L, ncol(data$returns)), common_length = min(120L,ncol(data$returns)),
      fund_rows = seq_len(min(8L,nrow(data$returns))), out_dir = subdir("constancy"))
  } else reproduce_constancy(data, B = bootstrap_reps, cores = cores, out_dir = subdir("constancy"))
  if ("calibration" %in% components) result$calibration <- if (quick) {
    reproduce_calibration(data, replications = 2L, cores = cores, lookbacks = 120L,
      train_end = ncol(data$returns)-2L, shifts = c(0,.002), out_dir = subdir("calibration"))
  } else reproduce_calibration(data, replications = calibration_reps, cores = cores,
                                out_dir = subdir("calibration"))
  if ("empirical" %in% components) result$empirical <- reproduce_empirical(data,
    lookbacks = if (quick) 120L else c(120L,180L), cores = cores,
    out_dir = subdir("empirical"))
  .paper_export_tables(result, output_dir)
  result$manifest <- list(version = as.character(utils::packageVersion("TVFundSelect")),
    manuscript_sha256 = paper_settings()$manuscript_sha256, components = components,
    quick = quick, synthetic_data = !is.null(data$synthetic) && isTRUE(data$synthetic),
    fresh_computation = TRUE, cores = cores, simulation_reps = if (quick) 2L else simulation_reps,
    calibration_reps = if (quick) 2L else calibration_reps,
    bootstrap_reps = if (quick) 9L else bootstrap_reps,
    seed = seed, R_version = R.version.string,
    caveat = "Quick mode is a software demonstration, never article evidence. Monte Carlo estimates vary with seeds.")
  saveRDS(result$manifest, file.path(output_dir, "run_manifest.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(result)
}
