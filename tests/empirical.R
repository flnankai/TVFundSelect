library(TVFundSelect)

near <- function(x, y, tolerance = 1e-9) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance))
must_error <- function(x) stopifnot(inherits(tryCatch(force(x), error = identity), "error"))
set.seed(7613)
TT <- 120L; N <- 12L
F <- matrix(rnorm(TT * 4L, sd = .02), TT, 4L)
rf <- rep(.002, TT)
raw <- matrix(.001, N, TT) + matrix(rep(drop(F %*% c(.8, .2, -.1, .15)), each = N), N, TT) +
  matrix(rnorm(N * TT, sd = .015), N, TT)
raw[1L, 3L] <- 0
fixture <- tempfile("fund_inputs_"); dir.create(fixture)
utils::write.table(raw, file.path(fixture, "fund.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
utils::write.table(F, file.path(fixture, "factor.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
utils::write.table(100 * rf, file.path(fixture, "rf.csv"), row.names = FALSE, col.names = FALSE, sep = ",")
data <- read_fund_panel(file.path(fixture, "fund.csv"), file.path(fixture, "factor.csv"), file.path(fixture, "rf.csv"))
stopifnot(!data$observed[1L, 3L], is.na(data$returns[1L, 3L]), near(data$rf, rf),
          identical(dim(data$returns), c(N, TT)), length(data$input_md5) == 3L)

# Independent intercept/HAC calculation, without calling the implementation's
# score formula through a duplicate function.
y <- raw[2L, ] - rf
X <- cbind(1, F)
lm_ref <- stats::lm(y ~ F)
coef_ref <- stats::coef(lm_ref)[1L]
bread <- solve(crossprod(X)); meat <- matrix(0, 5, 5)
obs_scores <- X * stats::residuals(lm_ref)
for (t in seq_len(TT)) meat <- meat + tcrossprod(obs_scores[t, ])
for (lag in 1:6) for (t in (lag + 1L):TT) {
  term <- tcrossprod(obs_scores[t, ], obs_scores[t - lag, ])
  meat <- meat + (1 - lag / 7) * (term + t(term))
}
cov_ref <- bread %*% meat %*% bread * TT / (TT - 5)
ev <- portfolio_alpha(raw[2L, ], F, rf, hac_lag = 6L)
stopifnot(near(ev$annual_alpha, 12 * coef_ref), near(ev$annual_se, 12 * sqrt(cov_ref[1L, 1L])))
cash <- portfolio_alpha(rf, F, rf)
stopifnot(cash$annual_alpha == 0, is.na(cash$p_value), cash$all_cash_or_zero_spread)

seed <- .Random.seed
result <- reproduce_empirical(data, lookbacks = 60L, holding = 20L, knot_candidates = 0L,
  cores = 1L, benchmarks = "FF4", sensitivity_knots = c(0L, 1L))
stopifnot(identical(.Random.seed, seed), result$fresh_computation,
  nrow(result$window_summary) == 6L, nrow(result$monthly_returns) == 120L,
  nrow(result$paired_alpha) == 1L, nrow(result$sensitivity_summary) == 4L,
  all(result$window_summary$hold_start > result$window_summary$train_end),
  !any(c("FDP", "TDP") %in% names(result$performance_summary)))

# Future data cannot change formation-date decisions or training eligibility.
changed <- data
changed$returns[, 61:TT] <- changed$returns[, 61:TT] + .2
changed$observed[2L, 61L] <- FALSE; changed$returns[2L, 61L] <- NA_real_
a <- TVFundSelect:::.emp_window_worker(1L, result$task_schedule, data, result$config)
b <- TVFundSelect:::.emp_window_worker(1L, result$task_schedule, changed, result$config)
stopifnot(identical(a$eligible, b$eligible), near(a$scores$p_value, b$scores$p_value))
held <- TVFundSelect:::.emp_hold(2L, 61:63, changed)
stopifnot(near(held$return, rf[61:63]), all(held$cash_weight == 1), sum(held$missing_events) == 1)

# Fresh fitted generator + fresh draws: no archive of residuals or results.
cal <- reproduce_calibration(data, replications = 2L, cores = 1L, lookbacks = 60L,
  train_end = 100L, shifts = c(0, .004), methods = c("TV-LS-BH", "CC-LS-BH"),
  signal_fraction = .2, knot_candidates = 0L)
stopifnot(identical(.Random.seed, seed), cal$fresh_computation,
  nrow(cal$replication_metrics) == 8L,
  all(cal$replication_metrics$FDP[cal$replication_metrics$m1 == 0] ==
      cal$replication_metrics$any_rejection[cal$replication_metrics$m1 == 0]),
  all(is.na(cal$calibration_summary$TDP[cal$calibration_summary$m1 == 0])),
  identical(cal$audit_draws[[1L]]$index, cal$audit_draws[[3L]]$index),
  identical(cal$audit_draws[[1L]]$signs, cal$audit_draws[[3L]]$signs))
design <- cal$calibrated_designs[[1L]]
stopifnot(max(abs(rowMeans(design$errors))) < 1e-12,
          max(abs(rowMeans(design$alpha_centered))) < 1e-12)
ref <- stats::lm.fit(cbind(1, make_tv_design(design$factors, 0L, 2L)$Z),
  t(sweep(data$returns[design$rows, design$train, drop = FALSE], 2L, data$rf[design$train], "-")))
stopifnot(near(t(ref$fitted.values) - design$original_alpha, design$null_mean))

must_error(reproduce_empirical(data, lookbacks = 121L, cores = 1L))
must_error(reproduce_calibration(data, train_end = 121L, cores = 1L))
must_error(portfolio_alpha(y, F, hac_lag = TT))
if (identical(Sys.getenv("TVFUNDSELECT_TEST_PARALLEL"), "true")) {
  reference <- reproduce_empirical(data, lookbacks = 60L, holding = 2L, knot_candidates = 0L,
    cores = 1L, benchmarks = "FF4", sensitivity = FALSE)
  parallel_fit <- reproduce_empirical(data, lookbacks = 60L, holding = 2L, knot_candidates = 0L,
    cores = 25L, benchmarks = "FF4", sensitivity = FALSE)
  stopifnot(identical(reference$fund_scores, parallel_fit$fund_scores),
    identical(reference$monthly_returns, parallel_fit$monthly_returns),
    identical(reference$strategy_alpha, parallel_fit$strategy_alpha), identical(.Random.seed, seed))
  cal_args <- list(data = data, replications = 13L, lookbacks = 60L, train_end = 100L,
    shifts = c(0, .004), methods = c("TV-LS-BH", "CC-LS-BH"),
    signal_fraction = .2, knot_candidates = 0L)
  ca <- do.call(reproduce_calibration, c(cal_args, list(cores = 1L)))
  cb <- do.call(reproduce_calibration, c(cal_args, list(cores = 25L)))
  stopifnot(identical(ca$replication_metrics, cb$replication_metrics),
    identical(ca$calibration_summary, cb$calibration_summary),
    identical(ca$audit_draws, cb$audit_draws), identical(.Random.seed, seed))
  cat("Empirical rolling and calibration tests also passed with 25 PSOCK workers.\n")
}
cat("Fresh empirical tests passed: raw-input units/missingness, independent HAC, rolling fits, no lookahead, holding rule, fresh generator/draws and RNG restoration.\n")
