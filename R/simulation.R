# Fresh Monte Carlo replication of the October 2026 manuscript. No archived
# summaries, checkpoints or original workspace paths are needed by this module.

.tfsim_integer <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be one finite integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

.tfsim_with_seed <- function(seed, fun) {
  seed <- .tfsim_integer(seed, "seed", 0L)
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  fun()
}

#' Generate a fresh panel from the article's simulation model
#' @export
simulate_fund_panel <- function(N = 200L, TT = 200L, example = 3L,
    distribution = c("normal", "t3", "mixture", "icm"), signal_shift = 0.02,
    rho = 0.5, nonnull_fraction = 0.1, seed = NULL, burn_in = 50L,
    positive_range = c(0.02, 0.34), alpha_deviation_amplitude = 0.3) {
  N <- .tfsim_integer(N, "N", 10L)
  TT <- .tfsim_integer(TT, "TT", 10L)
  example <- .tfsim_integer(example, "example")
  burn_in <- .tfsim_integer(burn_in, "burn_in", 0L)
  if (!example %in% 1:3) stop("example must be 1, 2 or 3.", call. = FALSE)
  distribution <- match.arg(distribution)
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("rho must lie strictly between -1 and 1.", call. = FALSE)
  }
  if (!is.numeric(nonnull_fraction) || length(nonnull_fraction) != 1L ||
      !is.finite(nonnull_fraction) || nonnull_fraction <= 0 || nonnull_fraction > 0.5 ||
      floor(N * nonnull_fraction) < 1) {
    stop("nonnull_fraction must produce nonempty disjoint positive/negative sets.", call. = FALSE)
  }
  if (!is.numeric(signal_shift) || length(signal_shift) != 1L || !is.finite(signal_shift)) {
    stop("signal_shift must be one finite number.", call. = FALSE)
  }
  if (!is.numeric(positive_range) || length(positive_range) != 2L ||
      any(!is.finite(positive_range)) || positive_range[1L] >= positive_range[2L] ||
      positive_range[1L] + signal_shift <= 0) {
    stop("positive_range and signal_shift must define strictly positive signals.", call. = FALSE)
  }
  if (!is.numeric(alpha_deviation_amplitude) || length(alpha_deviation_amplitude) != 1L ||
      !is.finite(alpha_deviation_amplitude) || alpha_deviation_amplitude < 0) {
    stop("alpha_deviation_amplitude must be nonnegative.", call. = FALSE)
  }
  generate <- function() simulate_tv_fund_panel(N = N, TT = TT, example = example,
    distribution = distribution, signal_shift = signal_shift, rho = rho,
    nonnull_fraction = nonnull_fraction, burn_in = burn_in,
    positive_range = positive_range, positive_signal_scale = sqrt(200 / TT),
    alpha_deviation = "centered_logistic",
    alpha_deviation_amplitude = alpha_deviation_amplitude)
  if (is.null(seed)) generate() else .tfsim_with_seed(seed, generate)
}

.tfsim_metrics <- function(selected, delta) {
  positive <- which(delta > 0)
  true <- sum(selected %in% positive)
  false <- length(selected) - true
  c(FDP = false / max(length(selected), 1L), TDP = true / length(positive),
    selected = length(selected), true_discoveries = true,
    false_discoveries = false, alternatives = length(positive))
}

.tfsim_one_task <- function(task, config) {
  panel <- simulate_fund_panel(N = config$N, TT = task$T, example = task$example,
    distribution = task$distribution, signal_shift = task$signal_shift,
    rho = config$rho, nonnull_fraction = config$nonnull_fraction,
    seed = task$task_seed)
  score_fit <- function(fit, method, multiplier = NA_real_) {
    if (any(!is.finite(fit$statistics)) || any(!is.finite(fit$p_values))) {
      stop("A fit returned non-finite statistics or p-values.", call. = FALSE)
    }
    metric <- .tfsim_metrics(fit$selected, panel$delta)
    out <- data.frame(task, N = config$N, method = method,
      bandwidth_multiplier = multiplier,
      bandwidth = if (is.null(fit$bandwidth)) NA_real_ else stats::median(fit$bandwidth),
      success = TRUE, stringsAsFactors = FALSE)
    for (name in names(metric)) out[[name]] <- unname(metric[[name]])
    out
  }
  tryCatch({
    if (config$experiment == "main") {
      fits <- fit_fund_methods(panel$Y, panel$factors, gamma = config$gamma,
        knot_candidates = config$knot_candidates, degree = config$degree,
        tolerance = config$tolerance, max_iter = config$max_iter,
        methods = config$methods)
      rows <- lapply(names(fits), function(name) score_fit(fits[[name]], name))
    } else {
      h0 <- 1.06 * stats::sd(seq_len(task$T) / task$T) * task$T^(-1/5)
      rows <- lapply(config$multipliers, function(multiplier) {
        fit <- ccy_bh(panel$Y, panel$factors, gamma = config$gamma,
                      bandwidth = multiplier * h0)
        score_fit(fit, "CCY-BH", multiplier)
      })
    }
    do.call(rbind, rows)
  }, error = function(e) stop("Simulation task ", task$task_id, " (example ", task$example,
    ", T=", task$T, ", ", task$distribution, ", s=", task$signal_shift,
    ", replication ", task$replication, ", seed ", task$task_seed, "): ",
    conditionMessage(e), call. = FALSE))
}

.tfsim_chunk <- function(indices, grid, config, output_dir, resume) {
  chunk_id <- as.integer(indices[1L])
  path <- if (is.null(output_dir)) NULL else file.path(output_dir, "checkpoints",
    sprintf("chunk_%09d.rds", chunk_id))
  if (!is.null(path) && file.exists(path)) {
    if (!resume) stop("Refusing to overwrite a completed simulation chunk.", call. = FALSE)
    saved <- readRDS(path)
    if (!identical(saved$task_ids, grid$task_id[indices])) {
      stop("Checkpoint task indices do not match this run.", call. = FALSE)
    }
    return(saved$raw)
  }
  raw <- do.call(rbind, lapply(indices, function(i) {
    .tfsim_one_task(grid[i, , drop = FALSE], config)
  }))
  if (!is.null(path)) {
    temporary <- paste0(path, ".", Sys.getpid(), ".tmp")
    saveRDS(list(task_ids = grid$task_id[indices], raw = raw), temporary)
    if (!file.rename(temporary, path)) stop("Could not commit simulation checkpoint.")
  }
  raw
}

.tfsim_summary <- function(raw) {
  keys <- c("example", "T", "N", "distribution", "signal_shift", "method")
  use_multiplier <- any(!is.na(raw$bandwidth_multiplier))
  if (use_multiplier) keys <- c(keys, "bandwidth_multiplier", "bandwidth")
  index <- interaction(raw[keys], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(seq_len(nrow(raw)), index), function(at) {
    out <- raw[at[1L], keys, drop = FALSE]
    out$replications <- length(at)
    for (metric in c("FDP", "TDP", "selected", "true_discoveries", "false_discoveries")) {
      values <- raw[[metric]][at]
      out[[paste0("mean_", metric)]] <- mean(values)
      out[[paste0("mcse_", metric)]] <- if (length(values) > 1L) {
        stats::sd(values) / sqrt(length(values))
      } else NA_real_
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[do.call(order, unname(out[keys])), , drop = FALSE]
}

.tfsim_prepare_output <- function(output_dir, manifest, resume) {
  if (is.null(output_dir)) {
    if (resume) stop("resume=TRUE requires output_dir.", call. = FALSE)
    return(NULL)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) stop("output_dir must be NULL or one nonempty path.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output_dir.", call. = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  manifest_path <- file.path(output_dir, "manifest.rds")
  if (file.exists(manifest_path)) {
    if (!resume) stop("Output already contains a run; choose a fresh output_dir or explicitly resume it.")
    if (!identical(readRDS(manifest_path), manifest)) {
      stop("Stored configuration, R version or package implementation differs; choose a fresh output_dir.")
    }
  } else {
    if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
      stop("Refusing to write into a nonempty directory without a matching manifest.")
    }
    saveRDS(manifest, manifest_path)
  }
  dir.create(file.path(output_dir, "checkpoints"), showWarnings = FALSE)
  output_dir
}

.tfsim_implementation <- function() {
  names <- c("simulate_tv_fund_panel", "simulate_reference_errors", "simulate_reference_factors",
    "simulate_reference_beta", "make_alpha_paths", "generate_average_alphas", "fit_fund_methods",
    "tv_ss_bh", "tv_ols_bh", "ss_bh_constant", "cc_ls_bh", "ccy_ll_bh", ".tfsim_one_task")
  # Keep literal algorithm text in the manifest: no external hash package or
  # machine-specific source path is needed to reject stale checkpoints.
  stats::setNames(lapply(names, function(name) {
    paste(deparse(get(name, envir = environment(.tfsim_implementation))), collapse = "\n")
  }), names)
}

.tfsim_run <- function(reps, cores, examples, times, distributions, signals, N, gamma,
    seed, rho, nonnull_fraction, knot_candidates, degree, tolerance, max_iter,
    methods, multipliers, experiment, output_dir, plots, resume, chunk_size) {
  restore_locale <- .tf_utf8_scope()
  on.exit(restore_locale(), add = TRUE)
  reps <- .tfsim_integer(reps, "reps")
  cores <- .tfsim_integer(cores, "cores")
  chunk_size <- .tfsim_integer(chunk_size, "chunk_size")
  N <- .tfsim_integer(N, "N", 10L)
  seed <- .tfsim_integer(seed, "seed", 0L)
  degree <- .tfsim_integer(degree, "degree", 1L)
  max_iter <- .tfsim_integer(max_iter, "max_iter")
  if (!is.numeric(examples) || !length(examples) || anyNA(examples) ||
      any(!examples %in% 1:3) || anyDuplicated(examples)) stop("Invalid or duplicated examples.")
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
      any(times < 10 | times != floor(times)) || anyDuplicated(times)) stop("Invalid or duplicated times.")
  if (!is.character(distributions) || !length(distributions) ||
      any(!distributions %in% c("normal", "t3", "mixture", "icm")) || anyDuplicated(distributions)) {
    stop("Invalid or duplicated error distributions.")
  }
  if (!is.numeric(signals) || !length(signals) || any(!is.finite(signals)) ||
      any(signals <= -0.02) || anyDuplicated(signals)) stop("Invalid or duplicated signal shifts.")
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 || gamma >= 0.5) {
    stop("gamma must lie in (0, 0.5).")
  }
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) stop("Invalid rho.")
  if (!is.numeric(nonnull_fraction) || length(nonnull_fraction) != 1L || !is.finite(nonnull_fraction) ||
      nonnull_fraction <= 0 || nonnull_fraction > 0.5 || floor(N * nonnull_fraction) < 1L) {
    stop("Invalid nonnull_fraction.")
  }
  if (!is.numeric(knot_candidates) || !length(knot_candidates) || any(!is.finite(knot_candidates)) ||
      any(knot_candidates < 0 | knot_candidates != floor(knot_candidates)) || anyDuplicated(knot_candidates)) {
    stop("Invalid knot_candidates.")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be positive.")
  }
  if (!is.character(methods) || !length(methods) || anyDuplicated(methods) ||
      any(!methods %in% c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"))) {
    stop("Invalid or duplicated methods.")
  }
  if (!is.numeric(multipliers) || !length(multipliers) || any(!is.finite(multipliers)) ||
      any(multipliers <= 0) || anyDuplicated(multipliers)) stop("Invalid bandwidth multipliers.")
  for (flag in list(plots, resume)) if (!is.logical(flag) || length(flag) != 1L || is.na(flag)) {
    stop("plots and resume must be single nonmissing logical values.")
  }
  grid <- expand.grid(distribution = distributions, signal_shift = signals,
    replication = seq_len(reps), T = as.integer(times), example = as.integer(examples),
    stringsAsFactors = FALSE)
  if (as.double(seed) + 1009 * nrow(grid) > .Machine$integer.max) {
    stop("The task seed schedule exceeds R's integer range; use a smaller seed or task grid.")
  }
  grid$task_id <- seq_len(nrow(grid))
  grid$task_seed <- as.integer(seed + 1009 * grid$task_id)
  config <- list(experiment = experiment, reps = reps, examples = as.integer(examples),
    times = as.integer(times), distributions = distributions, signals = signals,
    N = N, gamma = gamma, seed = seed, rho = rho, nonnull_fraction = nonnull_fraction,
    knot_candidates = as.integer(knot_candidates), degree = degree, tolerance = tolerance,
    max_iter = max_iter, methods = methods, multipliers = multipliers,
    chunk_size = chunk_size, seed_rule = "Mersenne-Twister/Inversion/Rejection; seed+1009*task_id",
    source = "TV_SS_BH_FDR_appendix_revised_20261024, Section 4 and Appendix D",
    fresh_computation = TRUE)
  version <- tryCatch(as.character(utils::packageVersion("TVFundSelect")), error = function(e) "development")
  manifest <- list(config = config, R_version = R.version.string,
    package_version = version, implementation = .tfsim_implementation())
  output_dir <- .tfsim_prepare_output(output_dir, manifest, resume)
  chunks <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / chunk_size))
  caps <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  old_caps <- Sys.getenv(caps, unset = NA_character_)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  cl <- NULL
  on.exit({
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    absent <- is.na(old_caps)
    if (any(absent)) Sys.unsetenv(caps[absent])
    if (any(!absent)) do.call(Sys.setenv, as.list(old_caps[!absent]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(caps)), caps)))
  message("Fresh ", experiment, " experiment: ", nrow(grid), " generated panels; ",
    if (experiment == "main") length(methods) else length(multipliers),
    " fits per panel; ", cores, " worker(s).")
  if (cores == 1L) {
    blocks <- lapply(chunks, .tfsim_chunk, grid, config, output_dir, resume)
  } else {
    cl <- parallel::makePSOCKcluster(cores, master = "127.0.0.1", setup_timeout = 30)
    initialize <- function(paths, vars, ctype) {
      if (.Platform$OS.type == "windows") {
        changed <- suppressWarnings(Sys.setlocale("LC_CTYPE", ctype))
        if (!nzchar(changed)) stop("Cannot initialize a UTF-8 worker locale.")
      }
      .libPaths(paths)
      do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(vars)), vars)))
      loadNamespace("TVFundSelect")
      NULL
    }
    environment(initialize) <- baseenv()
    parallel::clusterCall(cl, initialize, .libPaths(), caps, Sys.getlocale("LC_CTYPE"))
    worker <- function(indices, task_grid, cfg, output, restart) {
      get(".tfsim_chunk", envir = asNamespace("TVFundSelect"))(
        indices, task_grid, cfg, output, restart)
    }
    environment(worker) <- baseenv()
    blocks <- parallel::parLapplyLB(cl, chunks, worker, grid, config, output_dir, resume)
  }
  raw <- do.call(rbind, blocks)
  rownames(raw) <- NULL
  summary <- .tfsim_summary(raw)
  result <- structure(list(raw = raw, summary = summary, tasks = grid, config = config,
    cores = cores, manifest = manifest, output_dir = output_dir), class = "fund_simulation")
  if (!is.null(output_dir)) {
    utils::write.csv(raw, file.path(output_dir, "replications.csv"), row.names = FALSE)
    utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
    saveRDS(result, file.path(output_dir, "simulation.rds"))
    if (plots) plot_simulation_results(result, file = file.path(output_dir, "fdp_tdp.pdf"),
                                       overwrite = resume)
  }
  result
}

#' Recompute the five-method simulation study
#' @export
reproduce_simulations <- function(reps = 1000L, cores = 25L, examples = 1:3,
    times = c(200L, 400L), distributions = c("normal", "t3", "mixture", "icm"),
    signals = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.175, 0.23, 0.30),
    N = 200L, gamma = 0.1, seed = 420261010L, rho = 0.5,
    nonnull_fraction = 0.1, knot_candidates = 0:3, degree = 2L,
    tolerance = 1e-6, max_iter = 10000L,
    methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"),
    output_dir = NULL, plots = TRUE, resume = FALSE, chunk_size = 20L) {
  .tfsim_run(reps, cores, examples, times, distributions, signals, N, gamma,
    seed, rho, nonnull_fraction, knot_candidates, degree, tolerance, max_iter,
    methods, 1, "main", output_dir, plots, resume, chunk_size)
}

#' Recompute the CCY-BH bandwidth sensitivity study
#' @export
reproduce_bandwidth <- function(reps = 1000L, cores = 25L, examples = 1:3,
    times = c(200L, 400L), distributions = c("normal", "t3", "mixture", "icm"),
    signals = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.175, 0.23, 0.30),
    multipliers = c(1, 2, 4, 8, 16), N = 200L, gamma = 0.1,
    seed = 1120261010L, rho = 0.5, nonnull_fraction = 0.1,
    output_dir = NULL, plots = TRUE, resume = FALSE, chunk_size = 20L) {
  .tfsim_run(reps, cores, examples, times, distributions, signals, N, gamma,
    seed, rho, nonnull_fraction, 0:3, 2L, 1e-6, 10000L,
    "CCY-BH", multipliers, "bandwidth", output_dir, plots, resume, chunk_size)
}

print.fund_simulation <- function(x, ...) {
  cat("Fresh fund-selection Monte Carlo experiment\n")
  cat("  Experiment:", x$config$experiment, "\n")
  cat("  Newly generated panels:", nrow(x$tasks), "\n")
  cat("  Method/bandwidth fits:", nrow(x$raw), "\n")
  cat("  Settings:", nrow(x$summary), "method-level cells;",
      x$config$reps, "replications each\n")
  cat("  Workers:", x$cores, "\n")
  cat("  Mean FDP range:", paste(format(range(x$summary$mean_FDP), digits = 4), collapse = " -- "), "\n")
  cat("  Mean TDP range:", paste(format(range(x$summary$mean_TDP), digits = 4), collapse = " -- "), "\n")
  if (!is.null(x$output_dir)) cat("  Results:", x$output_dir, "\n")
  invisible(x)
}

#' Draw FDP/TDP curves calculated by a fresh replication run
#' @export
plot_simulation_results <- function(x, file = NULL, overwrite = FALSE) {
  restore_locale <- .tf_utf8_scope()
  on.exit(restore_locale(), add = TRUE)
  if (!inherits(x, "fund_simulation") || !is.data.frame(x$summary) || !nrow(x$summary)) {
    stop("x must be a nonempty result from reproduce_simulations or reproduce_bandwidth.")
  }
  opened <- FALSE
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) stop("Invalid file.")
    if (file.exists(file) && !isTRUE(overwrite)) stop("Refusing to overwrite figure: ", file)
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(file, width = 10, height = 11, onefile = TRUE)
    opened <- TRUE
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit({graphics::par(old); if (opened) grDevices::dev.off()}, add = TRUE)
  summary <- x$summary
  panels <- unique(summary[c("example", "T")])
  distributions <- x$config$distributions
  palette <- c("#1B6CA8", "#D1495B", "#3B7A57", "#7A5195", "#D18F00")
  for (k in seq_len(nrow(panels))) {
    graphics::par(mfrow = c(length(distributions), 2L), mar = c(3.4, 3.8, 2.2, 0.8),
                  oma = c(1.5, 0, 2.4, 0))
    for (distribution in distributions) for (metric in c("FDP", "TDP")) {
      z <- summary[summary$example == panels$example[k] & summary$T == panels$T[k] &
                     summary$distribution == distribution, , drop = FALSE]
      label <- if (x$config$experiment == "bandwidth") {
        paste0(z$bandwidth_multiplier, " h0")
      } else z$method
      levels <- if (x$config$experiment == "bandwidth") paste0(x$config$multipliers, " h0") else x$config$methods
      range_x <- range(x$config$signals)
      if (diff(range_x) == 0) range_x <- range_x + c(-0.01, 0.01)
      values <- z[[paste0("mean_", metric)]]
      graphics::plot(NA_real_, NA_real_, xlim = range_x,
        ylim = if (metric == "TDP") c(0, 1) else c(0, max(0.15, values) * 1.06),
        xlab = "Signal shift s", ylab = paste("Mean", metric),
        main = paste(distribution, metric))
      if (metric == "FDP") graphics::abline(h = x$config$gamma, lty = 2, col = "gray45")
      for (j in seq_along(levels)) {
        at <- which(label == levels[j])
        at <- at[order(z$signal_shift[at])]
        graphics::lines(z$signal_shift[at], values[at], type = "b",
          col = palette[(j - 1L) %% length(palette) + 1L], pch = j, lty = j)
      }
      if (identical(distribution, distributions[1L]) && metric == "FDP") {
        graphics::legend("topleft", legend = levels,
          col = rep(palette, length.out = length(levels)),
          pch = seq_along(levels), lty = seq_along(levels), cex = 0.68, bty = "n")
      }
    }
    graphics::mtext(paste0("Example ", panels$example[k], ", T = ", panels$T[k],
      "; ", x$config$reps, " newly calculated replications per setting"),
      outer = TRUE, side = 3, line = 0.5, cex = 1.1)
  }
  invisible(file)
}
