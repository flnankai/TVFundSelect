# Canonical paper-facing entry points. The historical implementation names
# remain internal so that archived replication modules retain their semantics.

.fund_method_label <- function(fit, label) {
  fit$implementation_label <- fit$implementation_label %||% fit$method
  fit$method <- label
  fit
}

tv_ls_bh <- function(Y, factors, gamma = 0.1, knot_candidates = 0:3,
                     degree = 2L, bic_penalty = c("paper", "standard")) {
  fit <- tv_ols_bh(Y, factors, gamma = gamma, knot_candidates = knot_candidates,
            degree = degree, bic_penalty = match.arg(bic_penalty),
            variance_estimator = "residual_df")
  .fund_method_label(fit, "TV-LS-BH")
}

cc_ss_bh <- function(Y, factors, gamma = 0.1, tolerance = 1e-6,
                     max_iter = 10000L) {
  fit <- ss_bh_constant(Y, factors, gamma = gamma, tolerance = tolerance,
                 max_iter = max_iter)
  .fund_method_label(fit, "CC-SS-BH")
}

ccy_bh <- function(Y, factors, gamma = 0.1, bandwidth = NULL,
    keep_curves = FALSE, B = 9999L, seed = 20261027L, cores = 25L,
    batch_size = 500L, keep_bootstrap = FALSE, keep_resamples = FALSE) {
  keep_curves <- .ccy_boot_flag(keep_curves, "keep_curves")
  keep_resamples <- .ccy_boot_flag(keep_resamples, "keep_resamples")
  fit <- ccy_bootstrap_bh(Y, factors, gamma = gamma, bandwidth = bandwidth,
    B = B, seed = seed, resampling = "residual", df_adjust = TRUE,
    batch_size = batch_size, cores = cores, keep_bootstrap = keep_bootstrap,
    keep_resamples = keep_curves || keep_resamples)
  if (keep_curves) {
    fit$alpha_hat <- sweep(fit$alpha_null, 1L, fit$delta_hat, "+")
    fit$residuals <- fit$original_residuals
  }
  fit$implementation_label <- "CCY-RB-BH"
  .fund_method_label(fit, "CCY-BH")
}

fit_fund_methods <- function(Y, factors, gamma = 0.1, knot_candidates = 0:3,
                             degree = 2L, tolerance = 1e-6,
                             max_iter = 10000L, ccy_bandwidth = NULL,
                             methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"),
                             ccy_B = 9999L, ccy_seed = 20261027L, ccy_cores = 25L) {
  fit_five_methods(Y, factors, gamma = gamma, knot_candidates = knot_candidates,
                   degree = degree, tolerance = tolerance, max_iter = max_iter,
                   ccy_bandwidth = ccy_bandwidth, methods = methods,
                   ccy_B = ccy_B, ccy_seed = ccy_seed, ccy_cores = ccy_cores)
}

make_tv_design <- function(factors, n_internal = 2L, degree = 2L) {
  build_tv_design(factors, n_internal = n_internal, degree = degree)
}
