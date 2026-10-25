library(TVFundSelect)
near <- function(a, b, tolerance = 2e-8) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance))
}
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

set.seed(4311)
N <- 12L
TT <- 120L
u <- seq_len(TT) / TT
F <- matrix(0.25 + rnorm(TT), ncol = 1L)
Y <- matrix(rnorm(N * TT), N, TT) +
  outer(seq(-0.2, 0.3, length.out = N), rep(1, TT)) +
  outer(rep(1, N), F[, 1L] * (0.5 + 0.3 * sin(2 * pi * u)))
fits <- fit_fund_methods(Y, F, knot_candidates = 1L)
stopifnot(identical(names(fits), c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH")))
standalone <- list(
  `TV-SS-BH` = tv_ss_bh(Y, F, knot_candidates = 1L),
  `TV-LS-BH` = tv_ls_bh(Y, F, knot_candidates = 1L),
  `CC-SS-BH` = cc_ss_bh(Y, F),
  `CC-LS-BH` = cc_ls_bh(Y, F),
  `CCY-BH` = ccy_bh(Y, F))
for (name in names(standalone)) {
  stopifnot(identical(standalone[[name]]$method, name),
            identical(standalone[[name]]$statistics, fits[[name]]$statistics),
            identical(standalone[[name]]$p_values, fits[[name]]$p_values),
            identical(standalone[[name]]$selected, fits[[name]]$selected))
}
stopifnot(identical(standalone[["TV-LS-BH"]]$implementation_label, "TV-OLS-BH-df"),
          identical(standalone[["CC-SS-BH"]]$implementation_label, "SS-BH-constant"),
          identical(standalone[["CCY-BH"]]$implementation_label, "CCY-BH"))
for (name in names(fits)) {
  fit <- fits[[name]]
  stopifnot(identical(fit$method, name), length(fit$p_values) == N,
            all(fit$p_values >= 0 & fit$p_values <= 1),
            near(fit$p_values, pnorm(fit$statistics, lower.tail = FALSE)),
            identical(fit$selected, bh_select(fit$p_values)$selected))
}

# Independent full OLS regression and intercept covariance, including the
# reduced centered spline block: the number of nuisance parameters excludes 1.
d <- make_tv_design(F, n_internal = 1L)
D <- cbind(1, d$Z)
reference <- lm.fit(D, t(Y))
rss <- colSums(reference$residuals^2)
sigma2 <- rss / (TT - ncol(D))
var_intercept <- solve(crossprod(D))[1L, 1L]
stopifnot(near(fits[["TV-LS-BH"]]$delta_hat, reference$coefficients[1L, ]),
          near(fits[["TV-LS-BH"]]$sigma2, sigma2),
          near(fits[["TV-LS-BH"]]$standard_errors, sqrt(sigma2 * var_intercept)),
          fits[["TV-LS-BH"]]$ols_residual_df == TT - d$rank - 1L)
constant_D <- cbind(1, F)
constant_reference <- lm.fit(constant_D, t(Y))
constant_sigma2 <- colSums(constant_reference$residuals^2) / (TT - ncol(constant_D))
stopifnot(near(fits[["CC-LS-BH"]]$delta_hat, constant_reference$coefficients[1L, ]),
          near(fits[["CC-LS-BH"]]$sigma2, constant_sigma2))

# Joint estimating equations, and each df calibration applied once only.
tvss <- fits[["TV-SS-BH"]]
stopifnot(tvss$spatial$R_loc <= 1e-6, tvss$spatial$R_scale <= 1e-6,
          near(tvss$statistics, sqrt((TT - d$rank) * tvss$varsigma_uncorrected) *
                 tvss$theta / tvss$d),
          near(tvss$variance_theta, tvss$d^2 / (TT * tvss$varsigma)))
ccss <- fits[["CC-SS-BH"]]
stopifnot(identical(ccss$calibration, "none"),
          near(ccss$statistics, sqrt(TT * ccss$varsigma) * ccss$theta / ccss$d))

# Local-linear fitted residuals and the literal T^-2 kernel variance objective.
h <- 0.4
local <- ccy_bh(Y, F, bandwidth = h, keep_curves = TRUE)
kernel <- 0.75 * pmax(1 - (outer(u, u, "-") / h)^2, 0)
kernel <- kernel / rowSums(kernel)
literal_var <- vapply(seq_len(N), function(i) {
  sum(sweep(kernel, 2L, local$residuals[i, ]^2, "*")) / TT^2
}, numeric(1))
stopifnot(near(local$standard_errors^2, literal_var),
          near(local$delta_hat, rowMeans(local$alpha_hat)),
          near(local$bandwidth, h))

# BIC is based on restricted residuals, not full/intercept OLS residuals.
chosen <- select_tv_design_bic(Y, F, knot_candidates = 0:2)
scores <- vapply(0:2, function(k) {
  z <- make_tv_design(F, k)$Z
  res <- lm.fit(z, t(Y))$residuals
  L <- k + 3L
  log(sum(res^2) / (N * TT)) + 2 * L * log(N * TT) / (N * TT)
}, numeric(1))
stopifnot(near(chosen$candidates$bic, scores),
          chosen$design$n_internal == which.min(scores) - 1L)

stopifnot(identical(bh_select(numeric())$selected, integer()),
          near(discovery_metrics(c(1L, 2L), c(0.2, 0, -0.1))$FDP, 0.5),
          near(discovery_metrics(c(1L, 2L), c(0.2, 0, -0.1))$TDP, 1))
fails(tv_ss_bh(Y[1L, , drop = FALSE], F))
fails(tv_ls_bh(Y, F, gamma = c(0.1, 0.2)))
fails(tv_ls_bh(Y, F, knot_candidates = 0.4))
fails(make_tv_design(F, n_internal = -1L))
fails(make_tv_design(F, degree = 2.1))
fails(fit_fund_methods(Y, F, methods = c("CC-LS-BH", "CC-LS-BH")))
fails(ccy_bh(Y, F, bandwidth = 0))
fails(bh_select(c(0.1, NA_real_)))
fails(discovery_metrics(1.2, rep(0, N)))
cat("Core methods: independent regression, spatial equations, kernel variance, BIC and validation passed.\n")
