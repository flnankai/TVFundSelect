# Constant-coefficient least-squares alpha inference followed by one-sided BH.
# Y is N x T; factors is T x p and does not contain an intercept.
# Uses the namespace's assert_finite_matrix() and bh_select() helpers.
# No caller-side source() operations are required.
#
# delta_hat = Y h/(h'h), h = M_F 1. The default variance uses the full
# intercept-plus-factor regression residuals with denominator T-p-1,
# matching the residual-df convention of the TV-LS-BH comparison.
# The separate "paper" option uses RSS/T, the original D-BH convention.

cc_ls_bh <- function(
    Y,
    factors,
    gamma = 0.1,
    variance_estimator = c("residual_df", "paper")) {
  .fund_panel(Y, factors, gamma, min_time = 2L)
  variance_estimator <- match.arg(variance_estimator)
  assert_finite_matrix(Y, "Y")
  assert_finite_matrix(factors, "factors")
  if (length(gamma) != 1L || !is.numeric(gamma) || !is.finite(gamma) ||
      gamma <= 0 || gamma >= 0.5) {
    stop("gamma must be one finite number in (0, 0.5).", call. = FALSE)
  }
  TT <- ncol(Y)
  pp <- ncol(factors)
  if (nrow(Y) < 1L || TT < 2L || nrow(factors) != TT) {
    stop("Y and factors must have compatible nonempty dimensions, with T >= 2.",
         call. = FALSE)
  }
  full_design_rank <- pp + 1L
  ols_residual_df <- TT - full_design_rank
  if (ols_residual_df <= 0L) {
    stop("The full constant-coefficient regression has no residual degrees of freedom.",
         call. = FALSE)
  }

  if (pp > 0L) {
    qf <- qr(factors, tol = 1e-10, LAPACK = FALSE)
    if (qf$rank != pp) {
      stop("The observed constant factor design is rank deficient.", call. = FALSE)
    }
    if (qr(cbind(1, factors), tol = 1e-10, LAPACK = FALSE)$rank != full_design_rank) {
      stop("The observed factors span the intercept, so alpha is not identified.",
           call. = FALSE)
    }
    Q <- qr.Q(qf, complete = FALSE)
    h <- drop(rep(1, TT) - Q %*% crossprod(Q, rep(1, TT)))
    restricted_residuals <- Y - (Y %*% Q) %*% t(Q)
  } else {
    h <- rep(1, TT)
    restricted_residuals <- Y
  }
  omega <- sum(h^2)
  if (!is.finite(omega) || omega <= TT * sqrt(.Machine$double.eps)) {
    stop("The residualized constant-regression intercept has zero or non-finite norm.",
         call. = FALSE)
  }
  delta_hat <- drop(Y %*% h) / omega
  # FWL: these unrestricted residuals remove the intercept AND all factors.
  residuals <- restricted_residuals - delta_hat %o% h
  rss <- rowSums(residuals^2)
  numerical_zero <- (100 * .Machine$double.eps)^2 * pmax(rowSums(Y^2), 1)
  if (any(!is.finite(rss)) || any(rss <= numerical_zero)) {
    stop("At least one constant-regression residual variance is zero or numerically degenerate.",
         call. = FALSE)
  }

  sigma2_df <- rss / ols_residual_df
  sigma2_paper <- rss / TT
  standard_errors_df <- sqrt(sigma2_df / omega)
  standard_errors_paper <- sqrt(sigma2_paper / omega)
  statistics_df <- delta_hat / standard_errors_df
  statistics_paper <- delta_hat / standard_errors_paper
  p_values_df <- stats::pnorm(statistics_df, lower.tail = FALSE)
  p_values_paper <- stats::pnorm(statistics_paper, lower.tail = FALSE)
  bh_df <- bh_select(p_values_df, gamma = gamma)
  bh_paper <- bh_select(p_values_paper, gamma = gamma)
  use_df <- variance_estimator == "residual_df"
  bh <- if (use_df) bh_df else bh_paper

  full_design <- cbind(intercept = 1, factors)
  normal_equations <- residuals %*% full_design
  normal_equation_scale <- 1 + sqrt(rowSums(Y^2)) %o% sqrt(colSums(full_design^2))
  max_normal_equation_residual <- max(abs(normal_equations) / normal_equation_scale)
  if (!is.finite(max_normal_equation_residual) || max_normal_equation_residual > 1e-8) {
    stop("Constant-regression OLS normal equations failed their scaled numerical check.",
         call. = FALSE)
  }

  list(
    selected = bh$selected, k = bh$k, bh_threshold = bh$threshold,
    delta_hat = delta_hat,
    sigma2 = if (use_df) sigma2_df else sigma2_paper,
    standard_errors = if (use_df) standard_errors_df else standard_errors_paper,
    statistics = if (use_df) statistics_df else statistics_paper,
    p_values = if (use_df) p_values_df else p_values_paper,
    sigma2_df = sigma2_df, sigma2_paper = sigma2_paper,
    standard_errors_df = standard_errors_df,
    standard_errors_paper = standard_errors_paper,
    statistics_df = statistics_df, statistics_paper = statistics_paper,
    p_values_df = p_values_df, p_values_paper = p_values_paper,
    selected_df = bh_df$selected, selected_paper = bh_paper$selected,
    rss = rss, residuals = residuals, omega = omega, omega_bar = omega / TT,
    design_rank = pp, full_design_rank = full_design_rank,
    residual_df = ols_residual_df, ols_residual_df = ols_residual_df,
    max_normal_equation_residual = max_normal_equation_residual,
    variance_estimator = variance_estimator, gamma = gamma,
    method = if (use_df) "CC-LS-BH" else "CC-LS-BH-paper"
  )
}
