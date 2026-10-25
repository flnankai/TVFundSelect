# Fresh empirical computations for Sections 5.3--5.5 of the accompanying paper.
# Historical result CSVs and cached method decisions are never inputs.

.emp_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) ||
      x < minimum || x > .Machine$integer.max) stop(name, " must be an integer >= ", minimum, ".")
  as.integer(x)
}

.emp_seed_state <- function() {
  list(had = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
       seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
         get(".Random.seed", envir = .GlobalEnv) else NULL, kind = RNGkind())
}

.emp_restore_seed <- function(state) {
  do.call(RNGkind, as.list(state$kind))
  if (state$had) assign(".Random.seed", state$seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    rm(".Random.seed", envir = .GlobalEnv)
}

.emp_map <- function(X, FUN, cores, ...) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  cores <- .emp_integer(cores, "cores")
  vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit({
    if (anyNA(old)) Sys.unsetenv(vars[is.na(old)])
    if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(vars)), vars)))
  if (cores == 1L || length(X) == 1L) return(lapply(X, FUN, ...))
  cl <- parallel::makePSOCKcluster(min(cores, length(X)), setup_timeout = 30)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  initialize <- function(paths, ctype) {
    if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ctype))
    .libPaths(paths)
    loadNamespace("TVFundSelect")
    NULL
  }
  environment(initialize) <- baseenv()
  parallel::clusterCall(cl, initialize, .libPaths(), Sys.getlocale("LC_CTYPE"))
  parallel::parLapply(cl, X, FUN, ...)
}

.emp_output_dir <- function(path) {
  if (is.null(path)) return(NULL)
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop("out_dir must be NULL or one nonempty path.")
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)))
    stop("Refusing to overwrite a nonempty output directory: ", path)
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop("Cannot create out_dir.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.emp_validate_data <- function(data) {
  if (!is.list(data) || !all(c("returns", "factors", "rf", "observed", "dates", "fund_id") %in% names(data)))
    stop("data must be a fund panel made by read_fund_panel().")
  y <- data$returns; f <- data$factors
  if (!is.matrix(y) || !is.numeric(y) || nrow(y) < 2L || ncol(y) < 3L ||
      !is.matrix(f) || !is.numeric(f) || nrow(f) != ncol(y) ||
      ncol(f) < 1L || any(!is.finite(f)) || length(data$rf) != ncol(y) ||
      any(!is.finite(data$rf)) || !identical(dim(data$observed), dim(y)) ||
      !is.logical(data$observed) || anyNA(data$observed) ||
      any(!is.finite(y[data$observed])) || any(y[data$observed] <= -1) ||
      length(data$dates) != ncol(y) || anyNA(data$dates) ||
      any(diff(as.numeric(data$dates)) <= 0) || length(data$fund_id) != nrow(y) ||
      anyNA(data$fund_id) || anyDuplicated(data$fund_id)) stop("Invalid fund panel dimensions or values.")
  invisible(data)
}

#' Read raw fund returns, factors and risk-free returns
#' @export
read_fund_panel <- function(fund_file, factor_file, rf_file, start_date = "1987-02-01",
    dates = NULL, zero_is_missing = TRUE, rf_in_percent = TRUE, header = FALSE) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  paths <- c(fund_file, factor_file, rf_file)
  if (length(paths) != 3L || anyNA(paths) || !all(file.exists(paths))) stop("All three input CSV files must exist.")
  if (!is.logical(zero_is_missing) || length(zero_is_missing) != 1L || is.na(zero_is_missing) ||
      !is.logical(rf_in_percent) || length(rf_in_percent) != 1L || is.na(rf_in_percent)) stop("Invalid input-unit flags.")
  read_mat <- function(path) {
    x <- as.matrix(utils::read.csv(path, header = header, check.names = FALSE))
    if (!is.numeric(x)) stop("Input CSV must contain only numeric values: ", basename(path))
    x
  }
  raw <- read_mat(fund_file); factors <- read_mat(factor_file); rf_matrix <- read_mat(rf_file)
  if (nrow(rf_matrix) != 1L && ncol(rf_matrix) != 1L) stop("rf_file must contain one row or one column.")
  rf <- as.numeric(rf_matrix) / if (rf_in_percent) 100 else 1
  observed <- is.finite(raw) & if (zero_is_missing) raw != 0 else TRUE
  returns <- raw; returns[!observed] <- NA_real_
  if (is.null(dates)) dates <- seq(as.Date(start_date), by = "month", length.out = ncol(raw))
  dates <- as.Date(dates)
  if (ncol(factors) == 4L) colnames(factors) <- c("Mkt-RF", "SMB", "HML", "MOM")
  out <- list(returns = returns, observed = observed, factors = factors, rf = rf,
    dates = dates, fund_id = sprintf("row_%04d", seq_len(nrow(raw))),
    input_md5 = stats::setNames(unname(tools::md5sum(paths)), c("fund.csv", "factor.csv", "rf.csv")),
    input_convention = list(zero_is_missing = zero_is_missing, rf_in_percent = rf_in_percent,
      fund_returns = "decimal net returns", factors = "decimal; no intercept", response = "net return minus RF"))
  .emp_validate_data(out)
  class(out) <- c("fund_panel", "list")
  out
}

.emp_hold <- function(selected, hold, data) {
  if (!length(selected)) return(data.frame(month_index = hold, return = data$rf[hold],
    cash_weight = 1, missing_events = 0L))
  sleeves <- rep(1 / length(selected), length(selected)); inactive <- rep(FALSE, length(selected))
  rows <- vector("list", length(hold))
  for (j in seq_along(hold)) {
    tt <- hold[j]; missing <- !data$observed[selected, tt]
    first_missing <- missing & !inactive; inactive <- inactive | missing
    r <- data$returns[selected, tt]; r[inactive] <- data$rf[tt]
    previous <- sum(sleeves); cash <- sum(sleeves[inactive]) / previous
    sleeves <- sleeves * (1 + r)
    rows[[j]] <- data.frame(month_index = tt, return = sum(sleeves) / previous - 1,
      cash_weight = cash, missing_events = sum(first_missing))
  }
  do.call(rbind, rows)
}

.emp_window_worker <- function(i, tasks, data, config) {
  task <- tasks[i, , drop = FALSE]
  train <- seq.int(task$train_first, task$train_last); hold <- seq.int(task$hold_first, task$hold_last)
  eligible <- which(rowSums(!data$observed[, train, drop = FALSE]) == 0L)
  if (length(eligible) < 2L) stop("Fewer than two training-complete funds in task ", i)
  Y <- sweep(data$returns[eligible, train, drop = FALSE], 2L, data$rf[train], "-")
  F <- data$factors[train, , drop = FALSE]
  fits <- fit_fund_methods(Y, F, gamma = config$gamma, knot_candidates = config$knot_candidates,
    degree = config$degree, tolerance = config$tolerance, max_iter = config$max_iter,
    methods = config$methods)
  monthly <- windows <- scores <- sensitivity <- list()
  for (method in names(fits)) {
    fit <- fits[[method]]; selected <- eligible[fit$selected]
    check <- which(stats::p.adjust(fit$p_values, "BH") <= config$gamma)
    if (!identical(as.integer(check), as.integer(fit$selected))) stop("BH identity failed.")
    held <- .emp_hold(selected, hold, data)
    key <- data.frame(L = task$L, block = task$block, task_id = i, method = method)
    monthly[[method]] <- cbind(key[rep(1L, nrow(held)), , drop = FALSE], held,
      date = as.character(data$dates[hold]), rf = data$rf[hold], n_selected = length(selected), n_eligible = length(eligible))
    windows[[method]] <- cbind(key, train_start = as.character(data$dates[min(train)]),
      train_end = as.character(data$dates[max(train)]), hold_start = as.character(data$dates[min(hold)]),
      hold_end = as.character(data$dates[max(hold)]), n_train = length(train), n_hold = length(hold),
      n_eligible = length(eligible), n_selected = length(selected),
      n_internal = if (is.null(fit$n_internal)) NA_integer_ else fit$n_internal)
    scores[[method]] <- cbind(key[rep(1L, length(eligible)), , drop = FALSE], fund_row = eligible,
      fund_id = data$fund_id[eligible], statistic = fit$statistics, p_value = fit$p_values,
      selected = eligible %in% selected)
    if (config$sensitivity && method %in% c("TV-SS-BH", "TV-LS-BH")) {
      for (nk in config$sensitivity_knots) {
        alternate <- fit_fund_methods(Y, F, gamma = config$gamma, knot_candidates = nk,
          degree = config$degree, tolerance = config$tolerance, max_iter = config$max_iter,
          methods = method)[[1L]]
        alt <- eligible[alternate$selected]; u <- union(selected, alt)
        sensitivity[[length(sensitivity) + 1L]] <- cbind(key, knots = nk,
          n_selected = length(alt), baseline_selected = length(selected),
          union_size = length(u), intersection_size = length(intersect(selected, alt)),
          overlap_active = if (length(u)) length(intersect(selected, alt)) / length(u) else NA_real_)
      }
    }
  }
  list(monthly = do.call(rbind, monthly), windows = do.call(rbind, windows),
    scores = do.call(rbind, scores), sensitivity = do.call(rbind, sensitivity), eligible = eligible)
}

#' Estimate portfolio alpha and Bartlett HAC uncertainty
#' @export
portfolio_alpha <- function(returns, factors, rf = rep(0, length(returns)),
    benchmark = c("FF4", "TV-FF4"), hac_lag = 6L, n_internal = 3L, degree = 2L,
    spread = FALSE) {
  benchmark <- match.arg(benchmark)
  if (!is.numeric(returns) || any(!is.finite(returns)) || !is.matrix(factors) ||
      !is.numeric(factors) || nrow(factors) != length(returns) || any(!is.finite(factors)) ||
      length(rf) != length(returns) || any(!is.finite(rf))) stop("Invalid evaluation returns or factors.")
  lag <- .emp_integer(hac_lag, "hac_lag", 0L)
  y <- returns - if (spread) 0 else rf
  Z <- if (benchmark == "TV-FF4") build_tv_design(factors, n_internal, degree)$Z else factors
  X <- cbind(intercept = 1, Z); n <- nrow(X); k <- ncol(X)
  if (n <= k || lag >= n) stop("Evaluation regression needs positive residual df and hac_lag < n.")
  fit <- stats::lm.fit(X, y)
  if (fit$rank != k) stop("Evaluation regression is rank deficient.")
  inv <- solve(crossprod(X)); score <- drop(X %*% inv[, 1L]) * fit$residuals
  variance <- sum(score^2)
  if (lag > 0L) for (j in seq_len(lag))
    variance <- variance + 2 * (1 - j / (lag + 1)) * sum(score[(j + 1L):n] * score[seq_len(n - j)])
  variance <- variance * n / (n - k)
  all_cash <- max(abs(y)) < 1e-14
  alpha <- if (all_cash) 0 else unname(fit$coefficients[1L])
  se <- if (all_cash) NA_real_ else sqrt(max(variance, 0))
  vol <- stats::sd(fit$residuals) * sqrt(12)
  data.frame(benchmark = benchmark, hac_lag = lag, n_months = n, n_coefficients = k,
    annual_alpha = 12 * alpha, annual_se = 12 * se,
    ci_low = 12 * (alpha - stats::qnorm(.975) * se), ci_high = 12 * (alpha + stats::qnorm(.975) * se),
    p_value = if (!is.finite(se) || se == 0) NA_real_ else 2 * stats::pnorm(-abs(alpha / se)),
    annual_residual_volatility = vol,
    information_ratio = if (all_cash || vol < 1e-14) NA_real_ else 12 * alpha / vol,
    all_cash_or_zero_spread = all_cash)
}

#' Recompute rolling fund selections and economic evaluation from raw inputs
#' @export
reproduce_empirical <- function(data, lookbacks = c(120L, 180L), holding = 6L,
    methods = c("TV-SS-BH", "TV-LS-BH"), gamma = .1, knot_candidates = 0:3,
    degree = 2L, tolerance = 1e-6, max_iter = 10000L, cores = 25L,
    benchmarks = c("FF4", "TV-FF4"), hac_lags = 6L, sensitivity = TRUE,
    sensitivity_knots = 0:3, out_dir = NULL) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  .emp_validate_data(data)
  holding <- .emp_integer(holding, "holding")
  cores <- .emp_integer(cores, "cores")
  if (!length(lookbacks) || any(!is.finite(lookbacks)) || any(lookbacks != floor(lookbacks)) ||
      any(lookbacks < 5) || any(lookbacks >= ncol(data$returns)) || anyDuplicated(lookbacks)) stop("Invalid lookbacks.")
  if (!length(methods) || anyDuplicated(methods) || any(!methods %in% five_method_names)) stop("Invalid methods.")
  if (!length(benchmarks) || any(!benchmarks %in% c("FF4", "TV-FF4")) ||
      any(!is.finite(hac_lags)) || any(hac_lags < 0) || any(hac_lags != floor(hac_lags))) stop("Invalid evaluation settings.")
  if (!is.logical(sensitivity) || length(sensitivity) != 1L || is.na(sensitivity)) stop("Invalid sensitivity flag.")
  out_dir <- .emp_output_dir(out_dir)
  state <- .emp_seed_state(); on.exit(.emp_restore_seed(state), add = TRUE)
  tasks <- list()
  for (L in as.integer(lookbacks)) {
    starts <- seq.int(1L, ncol(data$returns) - L, by = holding)
    for (b in seq_along(starts)) {
      first <- starts[b]
      tasks[[length(tasks) + 1L]] <- data.frame(L = L, block = b, train_first = first,
        train_last = first + L - 1L, hold_first = first + L,
        hold_last = min(first + L + holding - 1L, ncol(data$returns)))
    }
  }
  tasks <- do.call(rbind, tasks); tasks$task_id <- seq_len(nrow(tasks))
  config <- list(gamma = gamma, knot_candidates = knot_candidates, degree = degree,
    tolerance = tolerance, max_iter = max_iter, methods = methods, sensitivity = sensitivity,
    sensitivity_knots = sensitivity_knots, lookbacks = lookbacks, holding = holding,
    convention = "excess", benchmarks = benchmarks, hac_lags = hac_lags, cores = cores)
  runs <- .emp_map(seq_len(nrow(tasks)), .emp_window_worker, cores, tasks = tasks, data = data, config = config)
  monthly <- do.call(rbind, lapply(runs, `[[`, "monthly"))
  windows <- do.call(rbind, lapply(runs, `[[`, "windows"))
  scores <- do.call(rbind, lapply(runs, `[[`, "scores"))
  sensitivity_rows <- do.call(rbind, lapply(runs, `[[`, "sensitivity"))
  performance <- evaluation <- pairs <- list()
  for (L in lookbacks) for (method in methods) {
    ix <- monthly$L == L & monthly$method == method
    x <- monthly[ix, ]; x <- x[order(x$month_index), ]
    w <- windows[windows$L == L & windows$method == method, ]
    wealth <- cumprod(1 + x$return); drawdown <- wealth / cummax(c(1, wealth))[-1L] - 1
    monthly$wealth[ix] <- wealth; monthly$drawdown[ix] <- drawdown
    excess <- x$return - x$rf
    performance[[length(performance) + 1L]] <- data.frame(L = L, method = method,
      start = x$date[1L], end = tail(x$date, 1L), n_months = nrow(x),
      terminal_wealth = tail(wealth, 1L), annualized_return = tail(wealth, 1L)^(12 / nrow(x)) - 1,
      annualized_volatility = stats::sd(x$return) * sqrt(12),
      sharpe = if (stats::sd(excess) > 0) mean(excess) / stats::sd(excess) * sqrt(12) else NA_real_,
      max_drawdown = min(drawdown), mean_selected = mean(w$n_selected),
      active_windows = sum(w$n_selected > 0), total_windows = nrow(w),
      average_cash_weight = mean(x$cash_weight), holding_missing_events = sum(x$missing_events))
    for (bm in benchmarks) for (lag in hac_lags) {
      ev <- portfolio_alpha(x$return, data$factors[x$month_index, , drop = FALSE],
        data$rf[x$month_index], bm, lag)
      evaluation[[length(evaluation) + 1L]] <- cbind(L = L, method = method, ev)
    }
  }
  if (all(c("TV-SS-BH", "TV-LS-BH") %in% methods)) for (L in lookbacks) {
    a <- monthly[monthly$L == L & monthly$method == "TV-SS-BH", ]; a <- a[order(a$month_index), ]
    b <- monthly[monthly$L == L & monthly$method == "TV-LS-BH", ]; b <- b[order(b$month_index), ]
    if (!identical(a$month_index, b$month_index)) stop("Paired OOS dates do not agree.")
    for (bm in benchmarks) for (lag in hac_lags)
      pairs[[length(pairs) + 1L]] <- cbind(L = L, contrast = "TV-SS-BH minus TV-LS-BH",
        portfolio_alpha(a$return - b$return, data$factors[a$month_index, , drop = FALSE],
          benchmark = bm, hac_lag = lag, spread = TRUE))
  }
  sensitivity_summary <- NULL
  if (!is.null(sensitivity_rows)) {
    groups <- split(seq_len(nrow(sensitivity_rows)), interaction(sensitivity_rows$L,
      sensitivity_rows$method, sensitivity_rows$knots, drop = TRUE))
    sensitivity_summary <- do.call(rbind, lapply(groups, function(ix) {
      x <- sensitivity_rows[ix, ]; cbind(x[1L, c("L", "method", "knots")],
        mean_selected = mean(x$n_selected), active_windows = sum(x$n_selected > 0), total_windows = nrow(x),
        overlap_active = if (any(x$union_size > 0)) mean(x$overlap_active, na.rm = TRUE) else NA_real_)
    }))
  }
  out <- list(monthly_returns = monthly, window_summary = windows, fund_scores = scores,
    performance_summary = do.call(rbind, performance), strategy_alpha = do.call(rbind, evaluation),
    paired_alpha = do.call(rbind, pairs), sensitivity_by_window = sensitivity_rows,
    sensitivity_summary = sensitivity_summary, task_schedule = tasks, config = config,
    input_md5 = data$input_md5, R_version = R.version.string, fresh_computation = TRUE,
    interpretation = "Real fund truth is unknown; no empirical FDP/TDP is reported.")
  class(out) <- c("fund_empirical_reproduction", "list")
  if (!is.null(out_dir)) {
    for (nm in names(out)) if (is.data.frame(out[[nm]]))
      utils::write.csv(out[[nm]], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE, na = "NA")
    saveRDS(out, file.path(out_dir, "empirical_recomputed.rds"))
    plot_empirical(out, file.path(out_dir, "rolling_selection.pdf"))
    writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  }
  out
}

#' Plot freshly recomputed rolling selections and wealth
#' @export
plot_empirical <- function(x, file) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  if (!inherits(x, "fund_empirical_reproduction")) stop("x must come from reproduce_empirical().")
  grDevices::pdf(file, width = 10, height = 7)
  old <- graphics::par(mfrow = c(2, 2))
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  for (L in unique(x$monthly_returns$L)) {
    z <- x$monthly_returns[x$monthly_returns$L == L, ]; w <- x$window_summary[x$window_summary$L == L, ]
    ms <- unique(z$method); times <- sort(unique(z$month_index))
    wealth <- vapply(ms, function(m) z$wealth[z$method == m][order(z$month_index[z$method == m])],
      numeric(length(times)))
    wealth <- matrix(wealth, nrow = length(times), ncol = length(ms))
    graphics::matplot(times, wealth, type = "l", lty = seq_along(ms), col = seq_along(ms),
      xlab = "Month index", ylab = "Wealth", main = paste(L, "month training"))
    graphics::legend("topleft", ms, lty = seq_along(ms), col = seq_along(ms), bty = "n", cex = .7)
    blocks <- sort(unique(w$block))
    counts <- vapply(ms, function(m) w$n_selected[w$method == m][order(w$block[w$method == m])],
      numeric(length(blocks)))
    counts <- matrix(counts, nrow = length(blocks), ncol = length(ms))
    graphics::matplot(blocks, counts, type = "l", lty = seq_along(ms), col = seq_along(ms),
      xlab = "Formation block", ylab = "Selected funds", main = paste(L, "month training"))
  }
  invisible(file)
}
