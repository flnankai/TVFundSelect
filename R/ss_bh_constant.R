# Original constant-coefficient SS-BH of Wang, Zhao, Feng and Wang (2025).
# Uses the shared namespace's joint spatial-median solver and BH helper.
# Y is N x T; factors is T x p and contains the observed factors only.
# Reference: arXiv:2411.14016, Section 2, equations (4)-(7), Algorithm 1.

ss_bh_constant <- function(
    Y,
    factors,
    gamma = 0.1,
    tolerance = 1e-6,
    max_iter = 10000L) {
  .fund_panel(Y, factors, gamma, min_funds = 2L)
  assert_finite_matrix(Y, "Y")
  assert_finite_matrix(factors, "factors")
  TT <- ncol(Y)
  pp <- ncol(factors)
  if (nrow(factors) != TT) {
    stop("nrow(factors) must equal ncol(Y).", call. = FALSE)
  }
  if (TT <= pp + 1L) {
    stop("The constant factor regression has no residual degrees of freedom.", call. = FALSE)
  }

  # The restricted projection includes F only. Adding an intercept here would
  # remove the alpha signal before the spatial location was estimated.
  if (pp > 0L) {
    qf <- qr(factors, tol = 1e-10, LAPACK = FALSE)
    if (qf$rank != pp) {
      stop("The observed constant factor design is rank deficient.", call. = FALSE)
    }
    augmented_rank <- qr(cbind(1, factors), tol = 1e-10, LAPACK = FALSE)$rank
    if (augmented_rank != pp + 1L) {
      stop("The observed factors span the intercept, so alpha is not identified.", call. = FALSE)
    }
    restricted <- residualize_no_intercept(
      Y, list(Z = factors, qr = qf, rank = pp)
    )
    X <- restricted$X
    omega_bar <- restricted$omega_bar
  } else {
    X <- Y
    omega_bar <- 1
  }

  spatial <- fit_joint_spatial_median(
    X, tolerance = tolerance, max_iter = max_iter
  )
  radii <- spatial$radii
  zeta_minus1 <- mean(1 / radii)
  zeta1 <- mean(radii)
  zeta2 <- mean(radii^2)
  eta_direct <- 1 - 2 * (1 - omega_bar) * zeta_minus1 * zeta1 +
    (1 - omega_bar) * zeta_minus1^2 * zeta2
  # Algebraically identical sum of nonnegative terms avoids cancellation.
  eta <- omega_bar + (1 - omega_bar) * mean((1 - zeta_minus1 * radii)^2)
  if (!is.finite(eta) || eta <= 0) {
    stop("The original SS-BH radial normalizer is non-positive or non-finite.", call. = FALSE)
  }
  if (!isTRUE(all.equal(eta, eta_direct, tolerance = 1e-9))) {
    stop("Direct and stable SS-BH radial normalizers disagree.", call. = FALSE)
  }
  varsigma <- nrow(Y) * zeta_minus1^2 / eta
  statistics <- sqrt(TT * varsigma) * spatial$theta / spatial$d
  p_values <- stats::pnorm(statistics, lower.tail = FALSE)
  bh <- bh_select(p_values, gamma = gamma)

  list(
    selected = bh$selected,
    k = bh$k,
    bh_threshold = bh$threshold,
    statistics = statistics,
    p_values = p_values,
    theta = spatial$theta,
    d = spatial$d,
    varsigma = varsigma,
    variance_theta = spatial$d^2 / (TT * varsigma),
    eta = eta,
    eta_direct = eta_direct,
    zeta = c(minus1 = zeta_minus1, plus1 = zeta1, plus2 = zeta2),
    omega_bar = omega_bar,
    spatial = spatial,
    design_rank = pp,
    calibration = "none",
    score_multiplier = 1,
    gamma = gamma,
    method = "CC-SS-BH", implementation_label = "SS-BH-constant"
  )
}
