# Independent numerical verification of the CCY mean-alpha bootstrap extension.
# Run from the project root with Rscript tests/test_ccy_bootstrap.R.
# TVFUNDSELECT_LIB may point to an installed TVFundSelect library.
# This is a numerical test suite, not a Monte Carlo FDR experiment.

run_ccy_bootstrap_tests <- function() {
  old_locale <- Sys.getlocale("LC_CTYPE")
  old_libs <- .libPaths()
  original_kind <- RNGkind()
  original_had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  original_seed <- if (original_had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit({
    suppressWarnings(Sys.setlocale("LC_CTYPE", old_locale))
    .libPaths(old_libs)
    do.call(RNGkind, as.list(original_kind))
    if (original_had_seed) assign(".Random.seed", original_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (.Platform$OS.type == "windows") {
    changed <- suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
    if (!nzchar(changed)) stop("Cannot initialize UTF-8 locale.")
  }
  ccy_bootstrap_bh <- TVFundSelect::ccy_bootstrap_bh
  .ccy_boot_batch <- getFromNamespace(".ccy_boot_batch", "TVFundSelect")
  near <- function(x, y, tolerance = 2e-9) {
    identical(dim(x), dim(y)) && length(x) == length(y) &&
      all(is.finite(x)) && all(is.finite(y)) &&
      max(abs(as.numeric(x) - as.numeric(y)), 0) <=
        tolerance * max(1, abs(as.numeric(x)), abs(as.numeric(y)))
  }
  assert_near <- function(x, y, label, tolerance = 2e-9) {
    if (!near(x, y, tolerance)) {
      stop(label, "; max absolute difference = ",
           if (length(x) == length(y)) max(abs(as.numeric(x) - as.numeric(y))) else Inf)
    }
  }
  expect_error <- function(expr) {
    result <- tryCatch(force(expr), error = identity)
    if (!inherits(result, "error")) stop("Expected a validation error.")
    invisible(result)
  }
  section <- function(message) cat("PASS: ", message, "\n", sep = "")

  # Independent QR estimator: unscaled time slopes, no production smoother,
  # no production design cache, and no production bootstrap weight shortcut.
  reference_fit <- function(Y, factors, bandwidth) {
    n <- ncol(Y)
    u <- seq_len(n) / n
    X <- cbind(1, factors)
    p <- ncol(X)
    alpha <- fitted <- matrix(NA_real_, nrow(Y), n)
    kernel_rows <- matrix(0, n, n)
    for (j in seq_len(n)) {
      distance <- u - u[j]
      kernel <- 0.75 * pmax(1 - (distance / bandwidth)^2, 0)
      D <- cbind(X, X * distance)
      qr_fit <- stats::lm.wfit(D, t(Y), w = kernel)
      if (qr_fit$rank != ncol(D) || anyNA(qr_fit$coefficients)) {
        stop("Independent QR reference is not identifiable.")
      }
      coef <- matrix(qr_fit$coefficients, ncol(D), nrow(Y))
      alpha[, j] <- coef[1L, ]
      fitted[, j] <- drop(X[j, , drop = FALSE] %*% coef[seq_len(p), , drop = FALSE])
      kernel_rows[j, ] <- kernel / sum(kernel)
    }
    list(alpha = alpha, fitted = fitted, delta = rowMeans(alpha),
         residuals = Y - fitted, kernel_rows = kernel_rows)
  }
  bh_reference <- function(p, gamma) which(stats::p.adjust(p, "BH") <= gamma)
  upper_tail <- function(draws, observed) {
    (1 + vapply(seq_len(nrow(draws)), function(i) {
      sum(draws[i, ] >= observed[i])
    }, numeric(1L))) / (ncol(draws) + 1)
  }
  fit_args <- function(Y, factors, resampling = "residual", B = 47L,
                       batch_size = 13L, cores = 1L, ...) {
    ccy_bootstrap_bh(Y, factors, bandwidth = h, gamma = gamma, B = B,
      seed = 725019L, resampling = resampling, batch_size = batch_size,
      cores = cores, keep_bootstrap = TRUE, keep_resamples = TRUE, ...)
  }

  set.seed(2701026)
  TT <- 64L
  N <- 6L
  u <- seq_len(TT) / TT
  factors <- cbind(0.7 + rnorm(TT), -0.25 + rnorm(TT))
  alpha_deviation <- sin(2 * pi * u) + 0.5 * cos(4 * pi * u)
  alpha_deviation <- alpha_deviation - mean(alpha_deviation)
  Y <- matrix(rnorm(N * TT), N, TT) +
    c(-0.25, 0, 0.1, 0.45, 0.9, 1.5) %o% rep(1, TT) +
    seq(0.5, 1.5, length.out = N) %o% alpha_deviation +
    outer(seq_len(N) / N, factors[, 1L] * (1 + sin(3 * pi * u))) +
    outer(seq_len(N) / (2 * N), factors[, 2L] * cos(2 * pi * u))
  h <- 0.28
  gamma <- 0.1
  reference <- reference_fit(Y, factors, h)
  impulse <- reference_fit(diag(TT), factors, h)
  S <- t(impulse$fitted)
  A <- t(impulse$alpha)
  a <- colMeans(A)
  C <- diag(TT) - matrix(1 / TT, TT, TT)
  M <- diag(TT) - S
  centered_df <- sum((C %*% M)^2)
  uncentered_df <- sum(M^2)
  expected_pool <- sweep(reference$residuals, 1L,
                          rowMeans(reference$residuals), "-") * sqrt(TT / centered_df)
  expected_null <- sweep(reference$fitted, 1L, reference$delta, "-")
  expected_alpha_null <- sweep(reference$alpha, 1L, reference$delta, "-")
  expected_bias <- reference_fit(expected_null, factors, h)$delta
  # Nonzero curves and bootstrap bias are deliberate: a zero-alpha null or
  # accidental removal of estimated smoothing bias must fail this fixture.
  stopifnot(max(abs(expected_alpha_null)) > 0.5,
            max(abs(expected_bias)) > 1e-3,
            max(abs(rowMeans(expected_alpha_null))) < 1e-12)

  checked_fits <- list()
  for (resampling in c("residual", "signed_residual")) {
    fit <- fit_args(Y, factors, resampling)
    checked_fits[[resampling]] <- fit
    assert_near(fit$delta_hat, reference$delta, "Average alpha against QR")
    assert_near(fit$average_alpha_weights, a, "Average-alpha impulse weights")
    assert_near(fit$original_residuals, reference$residuals, "Residuals against QR")
    assert_near(fit$fitted_null, expected_null, "Mean-alpha fitted null")
    assert_near(fit$alpha_null, expected_alpha_null, "Retained time-varying null alpha")
    assert_near(fit$residual_df, centered_df, "Centered smoother residual df")
    assert_near(fit$residual_df_uncentered, uncentered_df, "Uncentered smoother residual df")
    assert_near(fit$residual_pool, expected_pool, "Centered, df-calibrated residual pool")
    assert_near(fit$bootstrap_bias, expected_bias, "Retained bootstrap smoothing bias")
    stopifnot(max(abs(rowMeans(fit$residual_pool))) < 1e-12,
              identical(dim(fit$bootstrap_indices), c(TT, 47L)),
              identical(dim(fit$bootstrap_signs), c(TT, 47L)),
              all(fit$bootstrap_indices >= 1L & fit$bootstrap_indices <= TT),
              all(fit$bootstrap_signs %in% c(-1L, 1L)))
    if (resampling == "residual") stopifnot(all(fit$bootstrap_signs == 1L))
    else stopifnot(any(fit$bootstrap_signs == -1L), any(fit$bootstrap_signs == 1L))

    # Reconstruct every pseudo-panel from shared date indices/common signs and
    # independently refit every local regression. No within-draw recentering.
    raw_reference <- matrix(NA_real_, N, fit$B)
    sample_means <- matrix(NA_real_, N, fit$B)
    for (b in seq_len(fit$B)) {
      noise <- sweep(expected_pool[, fit$bootstrap_indices[, b], drop = FALSE],
                     2L, fit$bootstrap_signs[, b], "*")
      pseudo <- expected_null + noise
      raw_reference[, b] <- reference_fit(pseudo, factors, h)$delta
      sample_means[, b] <- rowMeans(noise)
    }
    stopifnot(max(abs(sample_means)) > 0.05)
    assert_near(fit$bootstrap_statistics, raw_reference, "Bootstrap draws against full QR refits")
    centered_reference <- sweep(raw_reference, 1L, expected_bias, "-")
    assert_near(fit$p_values, upper_tail(raw_reference, reference$delta), "Plus-one raw-null tail")
    assert_near(fit$centered_p_values, upper_tail(centered_reference, reference$delta),
                "Plus-one centered diagnostic tail")
    stopifnot(identical(fit$selected, bh_reference(fit$p_values, gamma)),
              identical(fit$centered_selected, bh_reference(fit$centered_p_values, gamma)),
              all(fit$p_values >= 1 / (fit$B + 1)), all(fit$p_values <= 1))
    assert_near(fit$p_values * (fit$B + 1) - 1, fit$exceedances, "Raw tail counts")
    expected_se2 <- rowMeans(reference$residuals^2 %*% t(reference$kernel_rows)) / TT
    assert_near(fit$standard_errors^2, expected_se2, "Original analytic comparator variance")
    assert_near(fit$normal_p_values, stats::pnorm(reference$delta / sqrt(expected_se2),
                                                lower.tail = FALSE), "Normal comparator p-values")
    stopifnot(identical(fit$normal_selected, bh_reference(fit$normal_p_values, gamma)))
    section(paste(resampling, "independent QR refits, null curves, df, tails and BH"))
  }

  # A nonzero centered residual pool creates exactly representable tied draws.
  # Verify the actual computational kernel includes equality, not only >.
  stopifnot(identical(upper_tail(matrix(c(0, 0, 1), 1L), 0), 1),
            identical(upper_tail(matrix(c(0, 0, 1), 1L), 1), 0.5))
  tie_case <- .ccy_boot_batch(draw = matrix(c(1L, 2L, 1L, 1L, 2L, 2L), 2L),
    residual_pool = matrix(c(-1, 1), 1L), a = c(.5, .5),
    bootstrap_bias = 0, delta_hat = 0, resampling = "residual",
    keep_bootstrap = TRUE, keep_resamples = TRUE)
  stopifnot(identical(tie_case$statistics, matrix(c(0, -1, 1), 1L)),
            identical(tie_case$raw_hits, 2), identical(tie_case$centered_hits, 2))
  for (resampling in names(checked_fits)) {
    base <- checked_fits[[resampling]]
    no_df <- fit_args(Y, factors, resampling, df_adjust = FALSE)
    assert_near(no_df$residual_pool * base$residual_scale, base$residual_pool,
                "Degrees-of-freedom scaling")
    assert_near(sweep(no_df$bootstrap_statistics, 1L, no_df$bootstrap_bias, "-") *
                  base$residual_scale,
                sweep(base$bootstrap_statistics, 1L, base$bootstrap_bias, "-"),
                "Bootstrap draws respect degrees-of-freedom scaling")
    stopifnot(identical(no_df$residual_scale, 1))
  }
  section("Degrees-of-freedom correction changes only centered residual scale")

  # This exact equivariance permits one generated panel and one bootstrap
  # distribution to evaluate all nine shifts plus the global-null diagnostic.
  signals <- c(.02, .04, .06, .08, .10, .13, .175, .23, .30)
  truth <- c(-.25, 0, 0, .10, .20, .30)
  for (resampling in names(checked_fits)) {
    base <- checked_fits[[resampling]]
    shifts <- lapply(signals, function(s) (truth > 0) * (s - signals[1L]) * sqrt(200 / TT))
    shifts[[length(shifts) + 1L]] <- -truth
    for (shift in shifts) {
      shifted <- fit_args(Y + shift, factors, resampling)
      assert_near(shifted$delta_hat, base$delta_hat + shift, "Constant-shift estimator equivariance")
      assert_near(shifted$fitted_null, base$fitted_null, "Shift-invariant mean-alpha null")
      assert_near(shifted$residual_pool, base$residual_pool, "Shift-invariant residual pool")
      assert_near(shifted$bootstrap_statistics, base$bootstrap_statistics, "Shift-invariant bootstrap draws")
      assert_near(shifted$standard_errors, base$standard_errors, "Shift-invariant normal SE")
      assert_near(shifted$p_values, upper_tail(base$bootstrap_statistics, base$delta_hat + shift),
                  "Reused-null-distribution p-values")
    }
  }
  section("All nine signal shifts and the global-null shift agree with refitting")

  for (resampling in names(checked_fits)) {
    base <- checked_fits[[resampling]]
    rechunked <- fit_args(Y, factors, resampling, batch_size = 1L)
    longer <- fit_args(Y, factors, resampling, B = 79L, batch_size = 17L)
    stopifnot(identical(base$bootstrap_indices, rechunked$bootstrap_indices),
              identical(base$bootstrap_signs, rechunked$bootstrap_signs),
              identical(base$p_values, rechunked$p_values),
              identical(base$centered_p_values, rechunked$centered_p_values),
              identical(base$bootstrap_indices, longer$bootstrap_indices[, seq_len(base$B)]),
              identical(base$bootstrap_signs, longer$bootstrap_signs[, seq_len(base$B)]))
    assert_near(base$bootstrap_statistics, rechunked$bootstrap_statistics, "Batch-size invariant draws")
    assert_near(base$bootstrap_statistics, longer$bootstrap_statistics[, seq_len(base$B)],
                "B-prefix invariant draws")
  }
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(1234)
  caller_seed <- .Random.seed
  caller_kind <- RNGkind()
  invisible(fit_args(Y, factors, B = 19L))
  stopifnot(identical(caller_seed, .Random.seed), identical(caller_kind, RNGkind()))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(fit_args(Y, factors, B = 19L))
  stopifnot(!exists(".Random.seed", .GlobalEnv, inherits = FALSE),
            identical(caller_kind, RNGkind()))
  section("Chunk/B-prefix invariance and caller RNG state restoration")

  test_cores <- if (identical(Sys.getenv("TVFUNDSELECT_PARALLEL_TESTS"), "true")) 25L else 2L
  for (resampling in names(checked_fits)) {
    serial <- fit_args(Y, factors, resampling, B = 127L, batch_size = 5L)
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(7719)
    before_parallel_seed <- .Random.seed
    before_parallel_kind <- RNGkind()
    cap_names <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
    before_parallel_caps <- Sys.getenv(cap_names, unset = NA_character_)
    parallel_fit <- fit_args(Y, factors, resampling, B = 127L, batch_size = 5L, cores = test_cores)
    stopifnot(identical(serial$bootstrap_indices, parallel_fit$bootstrap_indices),
              identical(serial$bootstrap_signs, parallel_fit$bootstrap_signs),
              identical(serial$p_values, parallel_fit$p_values),
              identical(serial$centered_p_values, parallel_fit$centered_p_values),
              identical(serial$selected, parallel_fit$selected),
              identical(serial$centered_selected, parallel_fit$centered_selected),
              parallel_fit$cores_requested == test_cores, parallel_fit$cores_used == test_cores,
              length(unique(parallel_fit$worker_pids)) == test_cores,
              length(unique(parallel_fit$computing_worker_pids)) == test_cores,
              identical(before_parallel_seed, .Random.seed),
              identical(before_parallel_kind, RNGkind()),
              identical(before_parallel_caps, Sys.getenv(cap_names, unset = NA_character_)))
    assert_near(serial$bootstrap_statistics, parallel_fit$bootstrap_statistics,
                "Serial/25-worker bootstrap draws")
    section(paste(resampling, "serial vs actual 25-worker equality"))
  }

  # Check the exact shift shortcut on the paper's N=200 panels at both T values.
  for (paper_T in c(200L, 400L)) {
    panel <- TVFundSelect::simulate_fund_panel(N = 200L, TT = paper_T,
      example = 3L, distribution = "t3", signal_shift = .02, seed = 90100L + paper_T)
    common <- list(factors = panel$factors, B = 99L, seed = 777L, cores = 1L,
                   resampling = "signed_residual", keep_bootstrap = TRUE,
                   keep_resamples = TRUE)
    base <- do.call(ccy_bootstrap_bh, c(list(Y = panel$Y), common))
    for (shift in list((panel$delta > 0) * (.30 - .02) * sqrt(200 / paper_T), -panel$delta)) {
      shifted <- do.call(ccy_bootstrap_bh, c(list(Y = panel$Y + shift), common))
      assert_near(shifted$bootstrap_statistics, base$bootstrap_statistics,
                  "Paper-DGP shift-invariant bootstrap distribution", tolerance = 1e-8)
      assert_near(shifted$p_values, upper_tail(base$bootstrap_statistics, base$delta_hat + shift),
                  "Paper-DGP reused-draw p-values")
    }
  }
  section("Paper N=200, T=200/400 t3 panels: signal/global-null reuse matches new fits")

  expect_error(fit_args(Y, factors, B = 0L))
  expect_error(fit_args(Y, factors, B = 1.5))
  expect_error(fit_args(Y, factors, cores = 0L))
  expect_error(fit_args(Y, factors, cores = 1.5))
  expect_error(fit_args(Y, factors, batch_size = 0L))
  expect_error(ccy_bootstrap_bh(Y, factors, seed = NA_integer_, cores = 1L))
  expect_error(fit_args(Y, factors, df_adjust = NA))
  expect_error(ccy_bootstrap_bh(Y, factors, keep_resamples = NA, cores = 1L))
  expect_error(ccy_bootstrap_bh(Y, factors, gamma = .7, cores = 1L))
  expect_error(ccy_bootstrap_bh(Y, factors, bandwidth = 0, cores = 1L))
  expect_error(ccy_bootstrap_bh(Y, factors, bandwidth = 1 / TT, cores = 1L))
  expect_error(ccy_bootstrap_bh(Y, cbind(factors, factors[, 1L]), bandwidth = h, cores = 1L))
  bad_Y <- Y
  bad_Y[1L, 1L] <- NA_real_
  expect_error(ccy_bootstrap_bh(bad_Y, factors, cores = 1L))
  expect_error(ccy_bootstrap_bh(Y, factors[-1L, , drop = FALSE], cores = 1L))
  section("Invalid inputs and unidentified local fits fail explicitly")
  cat("ALL CCY BOOTSTRAP TESTS PASSED.\n")
  invisible(TRUE)
}

run_ccy_bootstrap_tests()
