# Core implementation of TV-SS-BH for the model
#   y_it = delta_i + lambda_i' z_t + error_it.
#
# Matrices follow the paper's orientation: Y is N x T and factors are T x p.

assert_finite_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix.", call. = FALSE)
  }
}

helmert_basis <- function(L) {
  if (L < 2L) stop("The spline basis must have at least two columns.", call. = FALSE)
  qr.Q(qr(stats::contr.helmert(L)), complete = FALSE)
}

make_bspline_basis <- function(u, n_internal = 2L, degree = 2L) {
  if (!is.numeric(u) || !length(u) || any(!is.finite(u)) || any(u < 0) || any(u > 1)) {
    stop("u must contain finite values in [0, 1].", call. = FALSE)
  }
  .fund_knots(n_internal, degree)
  n_internal <- as.integer(n_internal)
  degree <- as.integer(degree)
  if (n_internal < 0L || degree < 1L) {
    stop("n_internal must be nonnegative and degree must be positive.", call. = FALSE)
  }
  knots <- if (n_internal == 0L) NULL else {
    seq(0, 1, length.out = n_internal + 2L)[-c(1L, n_internal + 2L)]
  }
  B <- splines::bs(
    u,
    knots = knots,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  unname(as.matrix(B))
}

build_tv_design <- function(factors, n_internal = 2L, degree = 2L) {
  assert_finite_matrix(factors, "factors")
  TT <- nrow(factors)
  if (TT < 5L) stop("At least five time observations are required.", call. = FALSE)
  u <- seq_len(TT) / TT
  B <- make_bspline_basis(u, n_internal = n_internal, degree = degree)
  L <- ncol(B)

  B_tilde <- sweep(B, 2L, colMeans(B), FUN = "-")
  A_L <- helmert_basis(L)
  C <- B_tilde %*% A_L

  factor_blocks <- lapply(seq_len(ncol(factors)), function(j) {
    sweep(B, 1L, factors[, j], FUN = "*")
  })
  Z <- do.call(cbind, c(list(C), factor_blocks))
  colnames(Z) <- c(
    paste0("alpha_c", seq_len(ncol(C))),
    unlist(lapply(seq_len(ncol(factors)), function(j) {
      paste0("f", j, "_b", seq_len(L))
    }))
  )

  qz <- qr(Z, tol = 1e-10, LAPACK = FALSE)
  if (qz$rank != ncol(Z)) {
    stop(
      sprintf("Expanded spline design is rank deficient: rank %d < %d.", qz$rank, ncol(Z)),
      call. = FALSE
    )
  }
  q_aug <- qr(cbind(1, Z), tol = 1e-10, LAPACK = FALSE)
  if (q_aug$rank != ncol(Z) + 1L) {
    stop("Expanded design spans the intercept, so average alpha is not identified.", call. = FALSE)
  }

  list(
    Z = Z,
    B = B,
    B_tilde = B_tilde,
    C = C,
    A_L = A_L,
    L = L,
    rank = qz$rank,
    qr = qz,
    n_internal = as.integer(n_internal),
    degree = as.integer(degree)
  )
}

residualize_no_intercept <- function(Y, design) {
  assert_finite_matrix(Y, "Y")
  if (ncol(Y) != nrow(design$Z)) {
    stop("ncol(Y) must equal nrow(Z).", call. = FALSE)
  }
  Q <- qr.Q(design$qr, complete = FALSE)
  if (ncol(Q) > design$rank) Q <- Q[, seq_len(design$rank), drop = FALSE]
  X <- Y - (Y %*% Q) %*% t(Q)
  one_proj <- drop(crossprod(Q, rep(1, nrow(Q))))
  omega_bar <- 1 - sum(one_proj^2) / nrow(Q)
  if (!is.finite(omega_bar) || omega_bar <= sqrt(.Machine$double.eps)) {
    stop("The residualized intercept has zero or non-finite norm.", call. = FALSE)
  }
  list(X = X, Q = Q, omega_bar = omega_bar)
}

select_tv_design_bic <- function(
    Y,
    factors,
    knot_candidates = 0:3,
    degree = 2L,
    penalty = c("paper", "standard")) {
  assert_finite_matrix(Y, "Y")
  assert_finite_matrix(factors, "factors")
  .fund_knots(knot_candidates, degree)
  penalty <- match.arg(penalty)
  N <- nrow(Y)
  TT <- ncol(Y)
  if (N < 1L || TT < 5L) stop("Y must have at least one fund and five dates.", call. = FALSE)
  if (nrow(factors) != TT) stop("factors and Y have incompatible time dimensions.", call. = FALSE)

  fits <- lapply(as.integer(knot_candidates), function(nk) {
    tryCatch({
      design <- build_tv_design(factors, n_internal = nk, degree = degree)
      residual <- residualize_no_intercept(Y, design)
      mse <- sum(residual$X^2) / (N * TT)
      # Zhao-Wang's BIC counts (p + 1)(n + q) = (p + 1)L terms,
      # even though the draft removes one redundant centered-alpha column.
      paper_df <- (ncol(factors) + 1L) * design$L
      penalty_value <- if (penalty == "paper") {
        log(N * TT) * paper_df / (N * TT)
      } else {
        log(N * TT) * design$rank / TT
      }
      list(
        ok = TRUE,
        bic = log(mse) + penalty_value,
        mse = mse,
        design = design,
        residual = residual
      )
    }, error = function(e) {
      list(ok = FALSE, bic = Inf, message = conditionMessage(e), n_internal = nk)
    })
  })
  bics <- vapply(fits, `[[`, numeric(1), "bic")
  if (all(!is.finite(bics))) {
    msgs <- vapply(fits, function(x) x$message %||% "unknown error", character(1))
    stop("No admissible spline design: ", paste(unique(msgs), collapse = "; "), call. = FALSE)
  }
  best <- which.min(bics)
  out <- fits[[best]]
  out$candidates <- data.frame(
    n_internal = as.integer(knot_candidates),
    bic = bics,
    admissible = vapply(fits, `[[`, logical(1), "ok")
  )
  out$penalty <- penalty
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

spatial_equation_residuals <- function(X, theta, d) {
  standardized <- sweep(sweep(X, 1L, theta, FUN = "-"), 1L, d, FUN = "/")
  radii <- sqrt(colSums(standardized^2))
  if (any(!is.finite(radii)) || any(radii <= .Machine$double.eps)) {
    stop("A standardized residual radius is zero or non-finite.", call. = FALSE)
  }
  U <- sweep(standardized, 2L, radii, FUN = "/")
  N <- nrow(X)
  list(
    standardized = standardized,
    radii = radii,
    U = U,
    R_loc = sqrt(N) * max(abs(rowMeans(U))),
    R_scale = max(abs(N * rowMeans(U^2) - 1))
  )
}

fit_joint_spatial_median <- function(
    X,
    tolerance = 1e-6,
    max_iter = 10000L,
    min_scale = 1e-12) {
  assert_finite_matrix(X, "X")
  .fund_scalar(tolerance, "tolerance", 0, 0.5, closed_lower = FALSE, closed_upper = FALSE)
  .fund_scalar(max_iter, "max_iter", 1, .Machine$integer.max, integer = TRUE)
  .fund_scalar(min_scale, "min_scale", 0, closed_lower = FALSE)
  N <- nrow(X)
  TT <- ncol(X)
  if (N < 2L || TT < 3L) stop("X must have at least two rows and three columns.", call. = FALSE)
  if (!is.finite(tolerance) || tolerance <= 0 || tolerance >= 0.5) {
    stop("tolerance must be in (0, 0.5).", call. = FALSE)
  }

  theta <- rowMeans(X)
  d2 <- rowMeans((X - rowMeans(X))^2)
  if (any(!is.finite(d2)) || any(d2 <= min_scale)) {
    stop("At least one component has zero or non-finite initial variance.", call. = FALSE)
  }
  d2 <- exp(log(d2) - mean(log(d2)))

  last <- NULL
  for (iteration in seq_len(as.integer(max_iter))) {
    d <- sqrt(d2)
    eq <- spatial_equation_residuals(X, theta, d)
    last <- eq
    if (max(eq$R_loc, eq$R_scale) <= tolerance) {
      return(list(
        theta = theta,
        d = d,
        d2 = d2,
        radii = eq$radii,
        U = eq$U,
        standardized = eq$standardized,
        converged = TRUE,
        iterations = iteration - 1L,
        R_loc = eq$R_loc,
        R_scale = eq$R_scale,
        tolerance = tolerance
      ))
    }

    # Both updates deliberately use the same current standardized directions.
    theta_new <- theta + d * rowSums(eq$U) / sum(1 / eq$radii)
    d2_new <- N * d2 * rowMeans(eq$U^2)
    if (any(!is.finite(d2_new)) || any(d2_new <= min_scale)) {
      stop("The diagonal-scale update became zero or non-finite.", call. = FALSE)
    }
    d2_new <- exp(log(d2_new) - mean(log(d2_new)))
    theta <- theta_new
    d2 <- d2_new
  }

  stop(
    sprintf(
      "Joint spatial median did not converge in %d iterations (R_loc=%.3g, R_scale=%.3g).",
      max_iter, last$R_loc, last$R_scale
    ),
    call. = FALSE
  )
}

bh_select <- function(p_values, gamma = 0.1) {
  if (!is.numeric(p_values) || !is.null(dim(p_values)) ||
      any(!is.finite(p_values)) || any(p_values < 0) || any(p_values > 1)) {
    stop("p_values must lie in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 || gamma >= 0.5) {
    stop("gamma must lie in (0, 0.5).", call. = FALSE)
  }
  N <- length(p_values)
  if (!N) return(list(selected = integer(0), k = 0L, threshold = 0))
  ord <- order(p_values)
  eligible <- which(p_values[ord] <= gamma * seq_len(N) / N)
  k <- if (length(eligible)) max(eligible) else 0L
  selected <- if (k > 0L) which(p_values <= gamma * k / N) else integer(0)
  list(selected = selected, k = k, threshold = if (k > 0L) gamma * k / N else 0)
}

tv_ss_bh <- function(
    Y,
    factors,
    gamma = 0.1,
    knot_candidates = 0:3,
    degree = 2L,
    bic_penalty = c("paper", "standard"),
    tolerance = 1e-6,
    max_iter = 10000L,
    calibration = c("residual_df", "none")) {
  .fund_panel(Y, factors, gamma, min_funds = 2L)
  bic_penalty <- match.arg(bic_penalty)
  calibration <- match.arg(calibration)
  selected_design <- select_tv_design_bic(
    Y = Y,
    factors = factors,
    knot_candidates = knot_candidates,
    degree = degree,
    penalty = bic_penalty
  )
  spatial <- fit_joint_spatial_median(
    selected_design$residual$X,
    tolerance = tolerance,
    max_iter = max_iter
  )

  radii <- spatial$radii
  zeta_minus1 <- mean(1 / radii)
  zeta1 <- mean(radii)
  zeta2 <- mean(radii^2)
  omega_bar <- selected_design$residual$omega_bar

  eta_direct <- 1 - 2 * (1 - omega_bar) * zeta_minus1 * zeta1 +
    (1 - omega_bar) * zeta_minus1^2 * zeta2
  eta_stable <- omega_bar + (1 - omega_bar) * mean((1 - zeta_minus1 * radii)^2)
  if (!is.finite(eta_stable) || eta_stable <= 0) {
    stop("The radial-moment normalizer is non-positive or non-finite.", call. = FALSE)
  }
  if (!isTRUE(all.equal(eta_direct, eta_stable, tolerance = 1e-9))) {
    stop("Direct and stable evaluations of eta disagree.", call. = FALSE)
  }

  TT <- ncol(Y)
  residual_df <- TT - selected_design$design$rank
  if (residual_df <= 0L) stop("The nuisance design exhausts residual degrees of freedom.", call. = FALSE)
  df_ratio <- residual_df / TT
  varsigma_uncorrected <- nrow(Y) * zeta_minus1^2 / eta_stable
  statistics_uncorrected <- sqrt(TT * varsigma_uncorrected) * spatial$theta / spatial$d
  # q is rank(Z), excluding the intercept. The author-approved operational
  # calibration multiplies the variance of theta by T/(T-q), equivalently
  # multiplying the inverse-variance factor by (T-q)/T. Apply it exactly once.
  inverse_variance_multiplier <- if (calibration == "residual_df") df_ratio else 1
  score_multiplier <- sqrt(inverse_variance_multiplier)
  varsigma <- inverse_variance_multiplier * varsigma_uncorrected
  statistics <- score_multiplier * statistics_uncorrected
  p_values <- stats::pnorm(statistics, lower.tail = FALSE)
  bh <- bh_select(p_values, gamma = gamma)

  list(
    selected = bh$selected,
    k = bh$k,
    bh_threshold = bh$threshold,
    p_values = p_values,
    statistics = statistics,
    statistics_uncorrected = statistics_uncorrected,
    theta = spatial$theta,
    d = spatial$d,
    varsigma = varsigma,
    varsigma_uncorrected = varsigma_uncorrected,
    variance_theta = spatial$d^2 / (TT * varsigma),
    residual_df = residual_df,
    df_ratio = df_ratio,
    score_multiplier = score_multiplier,
    calibration = calibration,
    eta = eta_stable,
    eta_direct = eta_direct,
    zeta = c(minus1 = zeta_minus1, plus1 = zeta1, plus2 = zeta2),
    omega_bar = omega_bar,
    spatial = spatial,
    n_internal = selected_design$design$n_internal,
    spline_degree = selected_design$design$degree,
    spline_columns = selected_design$design$L,
    design_rank = selected_design$design$rank,
    bic_table = selected_design$candidates,
    gamma = gamma,
    method = if (calibration == "residual_df") "TV-SS-BH" else "TV-SS-BH-uncalibrated",
    implementation_label = if (calibration == "residual_df") "TV-SS-BH-df" else "TV-SS-BH"
  )
}

discovery_metrics <- function(selected, delta) {
  if (!is.numeric(delta) || !is.null(dim(delta)) || !length(delta) ||
      any(!is.finite(delta))) stop("delta must be a nonempty finite numeric vector.", call. = FALSE)
  if (!is.numeric(selected) || any(!is.finite(selected)) ||
      any(selected != floor(selected))) stop("selected must contain finite integer indices.", call. = FALSE)
  selected <- unique(as.integer(selected))
  if (length(selected) && (min(selected) < 1L || max(selected) > length(delta))) {
    stop("selected contains an invalid index.", call. = FALSE)
  }
  alternatives <- which(delta > 0)
  nulls <- which(delta <= 0)
  true_discoveries <- sum(selected %in% alternatives)
  false_discoveries <- sum(selected %in% nulls)
  list(
    FDP = false_discoveries / max(length(selected), 1L),
    TDP = true_discoveries / max(length(alternatives), 1L),
    selected = length(selected),
    true_discoveries = true_discoveries,
    false_discoveries = false_discoveries,
    alternatives = length(alternatives),
    nulls = length(nulls)
  )
}
