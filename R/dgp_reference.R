# Data-generating processes adapted from Zhao and Wang (2026), with the
# positive/negative average-alpha configuration of Wang et al. (2025).

logistic_G <- function(z, kappa1 = 2, kappa2 = 2) {
  stats::plogis(kappa1 * (z - kappa2))
}

simulate_ar_garch_factor <- function(
    TT,
    mu,
    ar,
    omega,
    garch,
    arch,
    burn_in = 50L) {
  # Index 1 represents t = -burn_in.  Simulate through t = T and retain
  # t = 1, ..., T, matching the paper's stated burn-in convention.
  total <- as.integer(TT + burn_in + 1L)
  innovation <- stats::rnorm(total)
  factor <- numeric(total)
  h <- numeric(total)
  factor[1L] <- 0
  h[1L] <- 1
  for (tt in 2:total) {
    h[tt] <- omega + garch * h[tt - 1L] + arch * h[tt - 1L] * innovation[tt - 1L]^2
    factor[tt] <- mu + ar * (factor[tt - 1L] - mu) + sqrt(h[tt]) * innovation[tt]
  }
  factor[(burn_in + 2L):total]
}

simulate_reference_factors <- function(TT, example = 3L, burn_in = 50L) {
  example <- as.integer(example)
  if (!example %in% 1:3) stop("example must be 1, 2, or 3.", call. = FALSE)
  specs <- data.frame(
    mu = c(0.34, 0.04, 0.06),
    ar = c(0.05, 0.07, 0.04),
    omega = c(0.32, 0.33, 0.26),
    garch = c(0.67, 0.51, 0.72),
    arch = c(0.13, 0.03, 0.05)
  )
  p <- if (example == 1L) 1L else 3L
  factors <- vapply(seq_len(p), function(j) {
    simulate_ar_garch_factor(
      TT = TT,
      mu = specs$mu[j],
      ar = specs$ar[j],
      omega = specs$omega[j],
      garch = specs$garch[j],
      arch = specs$arch[j],
      burn_in = burn_in
    )
  }, numeric(TT))
  if (is.null(dim(factors))) factors <- matrix(factors, ncol = 1L)
  colnames(factors) <- c("MKT", "SMB", "HML")[seq_len(p)]
  factors
}

simulate_reference_beta <- function(TT, example = 3L, burn_in = 50L) {
  example <- as.integer(example)
  u <- seq_len(TT) / TT
  if (example == 1L) {
    return(matrix(logistic_G(10 * u, 2, 2), ncol = 1L, dimnames = list(NULL, "MKT")))
  }
  if (example == 3L) {
    a <- c(0.5, 0.1, 0.2)
    b <- c(0.5, 0.5, 0.5)
    beta <- vapply(seq_len(3L), function(j) a[j] * logistic_G(10 * u, 2, 2) + b[j], numeric(TT))
    colnames(beta) <- c("MKT", "SMB", "HML")
    return(beta)
  }

  total <- TT + burn_in + 1L
  eps <- stats::rnorm(total)
  z <- numeric(total)
  sigma2 <- numeric(total)
  sigma2[1L] <- 1
  z[1L] <- 0
  for (tt in 2:total) {
    sigma2[tt] <- 0.1 + 0.3 * sigma2[tt - 1L]
    z[tt] <- 0.5 * z[tt - 1L] + sqrt(sigma2[tt]) * eps[tt]
  }
  z <- z[(burn_in + 2L):total]
  a <- c(0.8, 0.5, 0.6)
  b <- c(0.3, 0.1, 0.2)
  beta <- vapply(seq_len(3L), function(j) a[j] + b[j] * z, numeric(TT))
  colnames(beta) <- c("MKT", "SMB", "HML")
  beta
}

toeplitz_correlation <- function(N, rho = 0.5) {
  if (!is.finite(rho) || abs(rho) >= 1) stop("abs(rho) must be less than one.", call. = FALSE)
  stats::toeplitz(rho^(0:(N - 1L)))
}

simulate_reference_errors <- function(
    N,
    TT,
    rho = 0.5,
    distribution = c("normal", "t3", "mixture", "icm"),
    mixture_probability = 0.9) {
  distribution <- match.arg(distribution)
  Sigma <- toeplitz_correlation(N, rho)
  R <- chol(Sigma)

  if (distribution == "normal") {
    E <- matrix(stats::rnorm(TT * N), TT, N) %*% R
  } else if (distribution == "t3") {
    df <- 3
    gaussian <- matrix(stats::rnorm(TT * N), TT, N) %*% R
    # Unit-covariance revision: dividing the ordinary multivariate t by
    # sqrt(df/(df-2)) makes Cov(epsilon_t) exactly Sigma, not 3 * Sigma.
    E <- gaussian / sqrt(stats::rchisq(TT, df = df) / (df - 2))
  } else if (distribution == "mixture") {
    gaussian <- matrix(stats::rnorm(TT * N), TT, N) %*% R
    scale <- ifelse(stats::runif(TT) <= mixture_probability, 1, 3)
    E <- gaussian * scale / sqrt(mixture_probability + 9 * (1 - mixture_probability))
  } else {
    df <- 3
    components <- matrix(stats::rt(TT * N, df = df), TT, N) / sqrt(df / (df - 2))
    E <- components %*% R
  }
  t(E)
}

generate_average_alphas <- function(
    N,
    nonnull_fraction = 0.1,
    signal_shift = 0,
    positive_range = c(0.02, 0.34),
    negative_range = c(-0.30, -0.20),
    positive_signal_scale = 1) {
  if (length(positive_range) != 2L || any(!is.finite(positive_range)) ||
      positive_range[1L] >= positive_range[2L] ||
      !is.finite(positive_signal_scale) || positive_signal_scale <= 0 ||
      positive_range[1L] + signal_shift <= 0) {
    stop("Positive average-alpha range and scale must define strictly positive alternatives.", call. = FALSE)
  }
  n_each <- as.integer(floor(nonnull_fraction * N))
  if (n_each < 1L || 2L * n_each > N) {
    stop("nonnull_fraction must produce disjoint, nonempty positive and negative sets.", call. = FALSE)
  }
  positive <- sample.int(N, n_each)
  negative <- sample(setdiff(seq_len(N), positive), n_each)
  delta <- numeric(N)
  delta[positive] <- positive_signal_scale *
    (stats::runif(n_each, positive_range[1L], positive_range[2L]) + signal_shift)
  delta[negative] <- stats::runif(n_each, negative_range[1L], negative_range[2L])
  list(delta = delta, positive = positive, negative = negative)
}

make_alpha_paths <- function(
    delta,
    TT,
    deviation = c("centered_logistic", "none", "sine"),
    deviation_amplitude = 0.3) {
  deviation <- match.arg(deviation)
  alpha <- matrix(delta, nrow = length(delta), ncol = TT)
  if (deviation == "none" || deviation_amplitude == 0) return(alpha)
  u <- seq_len(TT) / TT
  g <- if (deviation == "centered_logistic") logistic_G(10 * u, 2, 2) else sin(2 * pi * u)
  g <- g - mean(g)
  amplitudes <- stats::runif(length(delta), -deviation_amplitude, deviation_amplitude)
  alpha + amplitudes %o% g
}

simulate_tv_fund_panel <- function(
    N = 200L,
    TT = 200L,
    example = 3L,
    rho = 0.5,
    distribution = c("normal", "t3", "mixture", "icm"),
    nonnull_fraction = 0.1,
    signal_shift = 0,
    alpha_deviation = c("centered_logistic", "none", "sine"),
    alpha_deviation_amplitude = 0.3,
    burn_in = 50L,
    positive_range = c(0.02, 0.34),
    positive_signal_scale = NULL) {
  distribution <- match.arg(distribution)
  alpha_deviation <- match.arg(alpha_deviation)
  if (is.null(positive_signal_scale)) positive_signal_scale <- sqrt(200 / TT)
  factors <- simulate_reference_factors(TT, example = example, burn_in = burn_in)
  beta <- simulate_reference_beta(TT, example = example, burn_in = burn_in)
  signal <- generate_average_alphas(
    N = N,
    nonnull_fraction = nonnull_fraction,
    signal_shift = signal_shift,
    positive_range = positive_range,
    positive_signal_scale = positive_signal_scale
  )
  alpha <- make_alpha_paths(
    signal$delta,
    TT = TT,
    deviation = alpha_deviation,
    deviation_amplitude = alpha_deviation_amplitude
  )
  errors <- simulate_reference_errors(
    N = N,
    TT = TT,
    rho = rho,
    distribution = distribution
  )
  systematic <- rowSums(beta * factors)
  Y <- alpha + matrix(systematic, nrow = N, ncol = TT, byrow = TRUE) + errors

  list(
    Y = Y,
    factors = factors,
    beta = beta,
    alpha = alpha,
    delta = signal$delta,
    positive = signal$positive,
    negative = signal$negative,
    errors = errors,
    settings = list(
      N = N,
      TT = TT,
      example = example,
      rho = rho,
      distribution = distribution,
      nonnull_fraction = nonnull_fraction,
      signal_shift = signal_shift,
      positive_range = positive_range,
      positive_signal_scale = positive_signal_scale,
      alpha_deviation = alpha_deviation,
      alpha_deviation_amplitude = alpha_deviation_amplitude
    )
  )
}
