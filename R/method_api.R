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
                   keep_curves = FALSE) {
  fit <- ccy_ll_bh(Y, factors, gamma = gamma, bandwidth = bandwidth,
            keep_curves = keep_curves)
  .fund_method_label(fit, "CCY-BH")
}

fit_fund_methods <- function(Y, factors, gamma = 0.1, knot_candidates = 0:3,
                             degree = 2L, tolerance = 1e-6,
                             max_iter = 10000L, ccy_bandwidth = NULL,
                             methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH")) {
  fit_five_methods(Y, factors, gamma = gamma, knot_candidates = knot_candidates,
                   degree = degree, tolerance = tolerance, max_iter = max_iter,
                   ccy_bandwidth = ccy_bandwidth, methods = methods)
}

make_tv_design <- function(factors, n_internal = 2L, degree = 2L) {
  build_tv_design(factors, n_internal = n_internal, degree = degree)
}
