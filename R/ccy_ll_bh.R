# Cai, Cheng and Yan (2018), Sections 2.2 and 2.4: local-linear coefficient
# estimation, Epanechnikov kernel, normal-reference bandwidth, and average alpha.
# The default long-run variance follows the cited Ang--Kristensen (2012)
# Theorem 2 / Eq. (17). The one-sided normal p-values and BH step are added here.
# Uses namespace validation and the common bh_select() helper.

.ccy_kernel_cache <- new.env(parent = emptyenv())

ccy_kernel_design <- function(TT, bandwidth) {
  key <- paste(TT, sprintf("%.17g", bandwidth), sep = ":")
  if (exists(key, envir = .ccy_kernel_cache, inherits = FALSE)) {
    return(get(key, envir = .ccy_kernel_cache, inherits = FALSE))
  }
  u <- seq_len(TT) / TT
  distance <- outer(u, u, function(target, observation) (observation - target) / bandwidth)
  kernel <- 0.75 * pmax(1 - distance^2, 0)
  weight_sum <- rowSums(kernel)
  normalized_kernel <- kernel / weight_sum
  out <- list(u = u, distance = distance, kernel = kernel,
              normalized_kernel = normalized_kernel,
              variance_weights = colMeans(normalized_kernel),
              support = lapply(seq_len(TT), function(t) which(kernel[t, ] > 0)))
  # Simulation workers use only T=200 and T=400. Limit storage if a caller
  # supplies many different bandwidths in a long interactive session.
  if (length(ls(.ccy_kernel_cache, all.names = TRUE)) >= 8L) {
    rm(list = ls(.ccy_kernel_cache, all.names = TRUE), envir = .ccy_kernel_cache)
  }
  assign(key, out, envir = .ccy_kernel_cache)
  out
}

ccy_ll_design <- function(factors, bandwidth, keep_alpha_smoother = FALSE) {
  TT <- nrow(factors)
  X <- cbind(intercept = 1, factors)
  P <- ncol(X)
  KD <- ccy_kernel_design(TT, bandwidth)
  fitted_smoother <- matrix(0, TT, TT)
  alpha_smoother <- if (keep_alpha_smoother) matrix(0, TT, TT) else NULL
  average_alpha_weights <- numeric(TT)
  min_local_observations <- TT
  max_reproduction_error <- 0
  for (t in seq_len(TT)) {
    index <- KD$support[[t]]
    min_local_observations <- min(min_local_observations, length(index))
    if (length(index) < 2L * P) {
      stop("Too few observations in a local-linear kernel window at t = ", t,
           ". Increase the bandwidth.", call. = FALSE)
    }
    x <- X[index, , drop = FALSE]
    scaled_time <- KD$distance[t, index]
    # Scaling the slope regressors by 1/h changes slope coefficients only;
    # it leaves all fitted values and alpha estimates exactly unchanged.
    D <- cbind(x, x * scaled_time)
    weights <- KD$kernel[t, index]
    gram <- crossprod(D, D * weights)
    rhs <- matrix(0, 2L * P, 2L)
    rhs[1L, 1L] <- 1
    rhs[seq_len(P), 2L] <- X[t, ]
    coefficients <- tryCatch(solve(gram, rhs), error = function(e) NULL)
    if (is.null(coefficients) || any(!is.finite(coefficients))) {
      stop("Singular local-linear regression at t = ", t, ".", call. = FALSE)
    }
    smoother <- (D %*% coefficients) * weights
    average_alpha_weights[index] <- average_alpha_weights[index] + smoother[, 1L] / TT
    fitted_smoother[t, index] <- smoother[, 2L]
    if (keep_alpha_smoother) alpha_smoother[t, index] <- smoother[, 1L]
    # A local-linear fit must reproduce all the columns in [X, (u-u_t)X].
    reproduction <- crossprod(smoother, D) - t(rhs)
    max_reproduction_error <- max(max_reproduction_error, abs(reproduction))
  }
  if (!is.finite(max_reproduction_error) || max_reproduction_error > 1e-7) {
    stop("Local-linear reproduction check failed.", call. = FALSE)
  }
  list(fitted_smoother = fitted_smoother,
       alpha_smoother = alpha_smoother,
       average_alpha_weights = average_alpha_weights,
       variance_weights = KD$variance_weights,
       min_local_observations = min_local_observations,
       max_reproduction_error = max_reproduction_error)
}

ccy_ll_bh <- function(Y, factors, gamma = 0.1, bandwidth = NULL,
                       keep_curves = FALSE) {
  .fund_panel(Y, factors, gamma)
  if (!is.logical(keep_curves) || length(keep_curves) != 1L || is.na(keep_curves)) {
    stop("keep_curves must be TRUE or FALSE.", call. = FALSE)
  }
  assert_finite_matrix(Y, "Y")
  assert_finite_matrix(factors, "factors")
  TT <- ncol(Y)
  if (TT < 5L || nrow(factors) != TT || nrow(Y) < 1L) {
    stop("Y and factors must have compatible nonempty dimensions, with T >= 5.",
         call. = FALSE)
  }
  bandwidth_rule <- if (is.null(bandwidth)) "CCY-normal-reference" else "user-specified"
  if (is.null(bandwidth)) {
    bandwidth <- 1.06 * stats::sd(seq_len(TT) / TT) * TT^(-1 / 5)
  }
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("bandwidth must be one finite positive number on the normalized time scale.",
         call. = FALSE)
  }
  design <- ccy_ll_design(factors, bandwidth, keep_alpha_smoother = keep_curves)
  fitted_values <- Y %*% t(design$fitted_smoother)
  residuals <- Y - fitted_values
  delta_hat <- drop(Y %*% design$average_alpha_weights)
  variance_weights <- design$variance_weights
  # AK (2012), Eq. (17), Delta=1, then Theorem 2: integrate the local
  # residual variance over target times and divide by T for Var(delta_hat).
  # Average normalized kernel rows before multiplication; no T x T x N array.
  long_run_variance <- drop((residuals^2) %*% variance_weights)
  standard_errors <- sqrt(long_run_variance / TT)
  if (any(!is.finite(standard_errors)) || any(standard_errors <= 0)) {
    stop("A local-linear long-run variance is zero or non-finite.", call. = FALSE)
  }
  statistics <- delta_hat / standard_errors
  p_values <- stats::pnorm(statistics, lower.tail = FALSE)
  bh <- bh_select(p_values, gamma)
  list(selected = bh$selected, k = bh$k, bh_threshold = bh$threshold,
       statistics = statistics, p_values = p_values, delta_hat = delta_hat,
       standard_errors = standard_errors, long_run_variance = long_run_variance,
       bandwidth = bandwidth, bandwidth_rule = bandwidth_rule,
       kernel = "Epanechnikov", variance_estimator = "ak2012",
       min_local_observations = design$min_local_observations,
       max_reproduction_error = design$max_reproduction_error,
       average_alpha_weights = design$average_alpha_weights,
       residual_variance_weights = variance_weights,
       residuals = if (keep_curves) residuals else NULL,
       alpha_hat = if (keep_curves) Y %*% t(design$alpha_smoother) else NULL,
       fitted_values = if (keep_curves) fitted_values else NULL,
       gamma = gamma,
       method = "CCY-BH")
}
