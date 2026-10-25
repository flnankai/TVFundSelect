# Least-squares average-alpha inference in the same centered spline model as
# TV-SS-BH. Shared design-selection and BH helpers resolve in the namespace.
#
# TV-D-BH extends Section 2 of arXiv:2411.14016 to the time-varying design Z:
#   delta_hat = Y M_Z 1 / (1' M_Z 1), sigma_paper^2 = RSS/T,
#   z = delta_hat * sqrt(1' M_Z 1) / sigma_paper, p = 1-Phi(z).
# The conventional residual-df version is explicitly a separate variant.

tv_ols_bh <- function(
    Y,
    factors,
    gamma = 0.1,
    knot_candidates = 0:3,
    degree = 2L,
    bic_penalty = c("paper", "standard"),
    variance_estimator = c("residual_df", "paper")) {
  .fund_panel(Y, factors, gamma)
  bic_penalty <- match.arg(bic_penalty)
  variance_estimator <- match.arg(variance_estimator)
  chosen <- select_tv_design_bic(
    Y, factors, knot_candidates = knot_candidates, degree = degree,
    penalty = bic_penalty
  )
  TT <- ncol(Y)
  q <- chosen$design$rank
  full_design_rank <- q + 1L
  ols_residual_df <- TT - full_design_rank
  if (ols_residual_df <= 0L) {
    stop("The full OLS regression has no residual degrees of freedom.", call. = FALSE)
  }

  # Frisch-Waugh projection: do not confuse mean(Y M_Z) with delta_hat.
  Q <- chosen$residual$Q
  a <- drop(rep(1, TT) - Q %*% crossprod(Q, rep(1, TT)))
  omega <- sum(a^2)
  if (!is.finite(omega) || omega <= TT * sqrt(.Machine$double.eps)) {
    stop("The residualized OLS intercept has zero or non-finite norm.", call. = FALSE)
  }
  delta_hat <- drop(Y %*% a) / omega
  # These are unrestricted residuals, after fitting the intercept AND Z.
  residuals <- chosen$residual$X - delta_hat %o% a
  rss <- rowSums(residuals^2)
  numerical_zero <- (100 * .Machine$double.eps)^2 * pmax(rowSums(Y^2), 1)
  if (any(!is.finite(rss)) || any(rss <= numerical_zero)) {
    stop("At least one OLS residual variance is zero or numerically degenerate.", call. = FALSE)
  }
  sigma2_df <- rss / ols_residual_df
  sigma2_paper <- rss / TT
  standard_errors_df <- sqrt(sigma2_df / omega)
  standard_errors_paper <- sqrt(sigma2_paper / omega)
  statistics_df <- delta_hat / standard_errors_df
  statistics_paper <- delta_hat / standard_errors_paper
  p_values_df <- stats::pnorm(statistics_df, lower.tail = FALSE)
  p_values_paper <- stats::pnorm(statistics_paper, lower.tail = FALSE)
  bh_df <- bh_select(p_values_df, gamma)
  bh_paper <- bh_select(p_values_paper, gamma)
  use_df <- variance_estimator == "residual_df"
  bh <- if (use_df) bh_df else bh_paper

  full_design <- cbind(intercept = 1, chosen$design$Z)
  normal_equations <- residuals %*% full_design
  normal_equation_scale <- 1 + sqrt(rowSums(Y^2)) %o% sqrt(colSums(full_design^2))
  max_normal_equation_residual <- max(abs(normal_equations) / normal_equation_scale)
  if (!is.finite(max_normal_equation_residual) || max_normal_equation_residual > 1e-8) {
    stop("OLS normal equations failed their scaled numerical check.", call. = FALSE)
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
    rss = rss, residuals = residuals,
    omega = omega, omega_bar = omega / TT,
    design_rank = q, full_design_rank = full_design_rank,
    ols_residual_df = ols_residual_df,
    n_internal = chosen$design$n_internal,
    spline_degree = chosen$design$degree,
    spline_columns = chosen$design$L,
    bic_table = chosen$candidates,
    max_normal_equation_residual = max_normal_equation_residual,
    variance_estimator = variance_estimator, gamma = gamma,
    method = if (use_df) "TV-LS-BH" else "TV-D-BH",
    implementation_label = if (use_df) "TV-OLS-BH-df" else "TV-D-BH"
  )
}
