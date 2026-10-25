# Conditional calibration from newly fitted raw-data means and residual pools.
# No archived fitted generator, residuals, selections or summary table is read.

.emp_calibration_designs <- function(data, config) {
  designs <- list(); metadata <- truth <- list()
  for (TT in config$lookbacks) {
    train <- seq.int(config$train_end - TT + 1L, config$train_end)
    rows <- which(rowSums(!data$observed[, train, drop = FALSE]) == 0L)
    if (length(rows) < 2L) stop("Calibration requires at least two training-complete funds.")
    Y <- sweep(data$returns[rows, train, drop = FALSE], 2L, data$rf[train], "-")
    F <- data$factors[train, , drop = FALSE]
    set.seed(config$truth_seed + TT)
    assignment <- sample.int(length(rows)); m <- floor(config$signal_fraction * length(rows))
    signs <- integer(length(rows))
    if (m > 0L) { signs[assignment[seq_len(m)]] <- 1L; signs[assignment[m + seq_len(m)]] <- -1L }
    for (dgp in config$generators) {
      if (dgp == "CC") {
        X <- cbind(1, F); nk <- NA_integer_; alpha_centered <- matrix(0, nrow(Y), TT)
      } else {
        selected <- select_tv_design_bic(Y, F, config$knot_candidates, config$degree)
        D <- selected$design; X <- cbind(1, D$Z); nk <- D$n_internal
      }
      qrx <- qr(X, tol = 1e-10, LAPACK = FALSE)
      if (qrx$rank != ncol(X) || TT <= qrx$rank) stop("Calibrated generating regression is unidentified.")
      coefficients <- t(qr.coef(qrx, t(Y))); fitted <- coefficients %*% t(X)
      if (dgp == "TV") {
        ca <- ncol(D$C)
        alpha_centered <- coefficients[, 1L + seq_len(ca), drop = FALSE] %*% t(D$C)
      }
      null_mean <- sweep(fitted, 1L, coefficients[, 1L], "-")
      errors <- Y - fitted
      errors <- sweep(errors, 1L, rowMeans(errors), "-") * sqrt(TT / (TT - qrx$rank))
      if (max(abs(rowMeans(errors))) > 1e-10 || max(abs(rowMeans(alpha_centered))) > 1e-10)
        stop("Calibrated residual/alpha centering check failed.")
      id <- paste0(dgp, "_", TT)
      designs[[id]] <- list(id = id, L = TT, dgp = dgp, N = nrow(Y), rows = rows,
        factors = F, null_mean = null_mean, errors = errors, alpha_centered = alpha_centered,
        signal_sign = signs, rank = qrx$rank, n_internal = nk,
        coefficients = coefficients, original_alpha = coefficients[, 1L], train = train)
      metadata[[id]] <- data.frame(design = id, dgp = dgp, L = TT, N = nrow(Y),
        train_start = as.character(data$dates[min(train)]), train_end = as.character(data$dates[max(train)]),
        factors = ncol(F), rank = qrx$rank, n_internal = nk,
        residual_scale = sqrt(TT / (TT - qrx$rank)), positive_funds = m, negative_funds = m)
      truth[[id]] <- data.frame(design = id, fund_row = rows, signal_sign = signs,
        alpha_original = coefficients[, 1L], centered_alpha_mean = rowMeans(alpha_centered),
        residual_mean = rowMeans(errors), residual_sd = sqrt(rowMeans(errors^2)))
    }
  }
  list(designs = designs, metadata = do.call(rbind, metadata), truth = do.call(rbind, truth))
}

.emp_calibration_draw <- function(design, shift, replication, config) {
  seed <- config$seed + design$L * 10000L + replication
  if (seed > .Machine$integer.max) stop("Calibration task seed exceeds integer range.")
  set.seed(seed)
  # block_length=1 is the paper's independent-month experiment. Independent
  # signs are shared across funds; the original pool is centered only once.
  nblocks <- ceiling(design$L / config$block_length)
  starts <- sample.int(design$L, nblocks, replace = TRUE)
  block_signs <- sample(c(-1, 1), nblocks, replace = TRUE)
  index <- unlist(lapply(starts, function(s) ((s - 1L + seq_len(config$block_length) - 1L) %% design$L) + 1L),
                  use.names = FALSE)[seq_len(design$L)]
  signs <- rep(block_signs, each = config$block_length)[seq_len(design$L)]
  delta <- design$signal_sign * shift
  errors <- sweep(design$errors[, index, drop = FALSE], 2L, signs, "*")
  list(Y = sweep(design$null_mean + errors, 1L, delta, "+"), delta = delta,
       seed = seed, index = index, signs = signs)
}

.emp_calibration_worker <- function(i, jobs, designs, config) {
  job <- jobs[i, , drop = FALSE]; design <- designs[[job$design]]
  panel <- .emp_calibration_draw(design, job$shift, job$replication, config)
  fits <- fit_fund_methods(panel$Y, design$factors, gamma = config$gamma,
    knot_candidates = config$knot_candidates, degree = config$degree,
    tolerance = config$tolerance, max_iter = config$max_iter, methods = config$methods)
  rows <- lapply(names(fits), function(method) {
    fit <- fits[[method]]; selected <- fit$selected; R <- length(selected)
    V <- sum(panel$delta[selected] <= 0); S <- sum(panel$delta[selected] > 0); m1 <- sum(panel$delta > 0)
    if (!identical(as.integer(selected), as.integer(which(stats::p.adjust(fit$p_values, "BH") <= config$gamma))))
      stop("Calibration BH identity failed.")
    data.frame(case_id = job$case_id, design = job$design, dgp = design$dgp,
      L = design$L, N = design$N, shift = job$shift, replication = job$replication,
      seed = panel$seed, method = method, R = R, V = V, S = S, m1 = m1,
      FDP = V / max(R, 1L), TDP = if (m1) S / m1 else NA_real_,
      any_rejection = as.integer(R > 0L), zero_selection = as.integer(R == 0L))
  })
  list(metrics = do.call(rbind, rows),
    audit = if (job$replication <= 3L) list(index = panel$index, signs = panel$signs,
      sum = sum(panel$Y), sumsq = sum(panel$Y^2)) else NULL)
}

#' Recompute the paper's data-calibrated false-discovery experiment
#' @export
reproduce_calibration <- function(data, replications = 1000L, cores = 25L,
    lookbacks = c(120L, 180L), train_end = 366L, shifts = c(0, .002, .004),
    generators = "TV", methods = five_method_names, signal_fraction = .1,
    gamma = .1, knot_candidates = 0:3, degree = 2L, block_length = 1L,
    seed = 31019000L, truth_seed = 20261019L, tolerance = 1e-6,
    max_iter = 10000L, out_dir = NULL) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  .emp_validate_data(data)
  replications <- .emp_integer(replications, "replications"); cores <- .emp_integer(cores, "cores")
  train_end <- .emp_integer(train_end, "train_end"); block_length <- .emp_integer(block_length, "block_length")
  seed <- .emp_integer(seed, "seed", 0L); truth_seed <- .emp_integer(truth_seed, "truth_seed", 0L)
  if (train_end > ncol(data$returns) || !length(lookbacks) || any(!is.finite(lookbacks)) ||
      any(lookbacks < 5) || any(lookbacks != floor(lookbacks)) || any(lookbacks > train_end) || anyDuplicated(lookbacks))
    stop("Invalid calibration lookbacks or train_end.")
  if (block_length > min(lookbacks) || !length(shifts) || any(!is.finite(shifts)) || any(shifts < 0) || anyDuplicated(shifts))
    stop("Invalid shifts or block_length.")
  if (!length(generators) || any(!generators %in% c("TV", "CC")) || anyDuplicated(generators) ||
      !length(methods) || any(!methods %in% five_method_names) || anyDuplicated(methods)) stop("Invalid generators or methods.")
  if (!is.numeric(signal_fraction) || length(signal_fraction) != 1L || !is.finite(signal_fraction) ||
      signal_fraction < 0 || signal_fraction > .5) stop("signal_fraction must lie in [0,.5].")
  if (max(lookbacks) + truth_seed > .Machine$integer.max ||
      seed + max(lookbacks) * 10000 + replications > .Machine$integer.max) stop("Seed range overflow.")
  out_dir <- .emp_output_dir(out_dir)
  state <- .emp_seed_state(); on.exit(.emp_restore_seed(state), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  config <- list(replications = replications, cores = cores, lookbacks = as.integer(lookbacks),
    train_end = train_end, shifts = shifts, generators = generators, methods = methods,
    signal_fraction = signal_fraction, gamma = gamma, knot_candidates = knot_candidates, degree = degree,
    block_length = block_length, seed = seed, truth_seed = truth_seed, tolerance = tolerance, max_iter = max_iter,
    convention = "excess", residual_scaling = "sqrt(T/(T-rank(full generating design)))")
  built <- .emp_calibration_designs(data, config)
  cases <- expand.grid(design = names(built$designs), shift = shifts, stringsAsFactors = FALSE)
  cases$case_id <- seq_len(nrow(cases))
  jobs <- do.call(rbind, lapply(seq_len(nrow(cases)), function(i)
    cbind(cases[rep(i, replications), , drop = FALSE], replication = seq_len(replications))))
  rownames(jobs) <- NULL
  runs <- .emp_map(seq_len(nrow(jobs)), .emp_calibration_worker, cores,
    jobs = jobs, designs = built$designs, config = config)
  metrics <- do.call(rbind, lapply(runs, `[[`, "metrics")); rownames(metrics) <- NULL
  groups <- split(seq_len(nrow(metrics)), interaction(metrics$case_id, metrics$method, drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(ix) {
    x <- metrics[ix, ]; cbind(x[1L, c("case_id", "design", "dgp", "L", "N", "shift", "method", "m1")],
      replications = nrow(x), FDR = mean(x$FDP), TDP = if (x$m1[1L]) mean(x$TDP) else NA_real_,
      mean_selected = mean(x$R), mean_false = mean(x$V), mean_true = mean(x$S),
      any_rejection = mean(x$any_rejection), zero_selection = mean(x$zero_selection),
      mcse_FDR = stats::sd(x$FDP) / sqrt(nrow(x)),
      mcse_TDP = if (x$m1[1L]) stats::sd(x$TDP) / sqrt(nrow(x)) else NA_real_)
  }))
  if (any(metrics$FDP[metrics$m1 == 0] != metrics$any_rejection[metrics$m1 == 0]))
    stop("All-null FDP identity failed.")
  out <- list(calibration_summary = summary, replication_metrics = metrics,
    design_metadata = built$metadata, simulation_truth = built$truth,
    calibrated_designs = built$designs, cases = cases, job_schedule = jobs,
    audit_draws = lapply(runs, `[[`, "audit"), config = config, input_md5 = data$input_md5,
    R_version = R.version.string, fresh_computation = TRUE,
    interpretation = paste("Conditional known-truth experiment; generating means/residuals freshly fitted from raw data.",
      "Whole residual vectors are resampled with a common Rademacher sign; no per-panel recentering.",
      "block_length=1 removes temporal dependence; this is not an estimate of real-data FDR."))
  class(out) <- c("fund_calibration_reproduction", "list")
  if (!is.null(out_dir)) {
    for (nm in names(out)) if (is.data.frame(out[[nm]]))
      utils::write.csv(out[[nm]], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE, na = "NA")
    saveRDS(out, file.path(out_dir, "calibration_recomputed.rds"))
    writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  }
  out
}
