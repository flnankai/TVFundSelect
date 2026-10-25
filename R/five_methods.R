# Unified public labels for the five comparators, resolved in the namespace.
# Historical source files and results keep their original labels for provenance.
five_method_names <- c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH")

fit_five_methods <- function(Y, factors, gamma = .1, knot_candidates = 0:3,
    degree = 2L, tolerance = 1e-6, max_iter = 10000L,
    ccy_bandwidth = NULL, methods = five_method_names) {
  if (!is.character(methods) || !length(methods) || anyNA(methods) || anyDuplicated(methods) ||
      any(!methods %in% five_method_names)) stop("Invalid or duplicated method names.")
  fits <- list(
    `TV-SS-BH` = function() tv_ss_bh(Y, factors, gamma = gamma,
      knot_candidates = knot_candidates, degree = degree,
      tolerance = tolerance, max_iter = max_iter),
    `TV-LS-BH` = function() tv_ols_bh(Y, factors, gamma = gamma,
      knot_candidates = knot_candidates, degree = degree,
      variance_estimator = "residual_df"),
    `CC-SS-BH` = function() ss_bh_constant(Y, factors, gamma = gamma,
      tolerance = tolerance, max_iter = max_iter),
    `CC-LS-BH` = function() cc_ls_bh(Y, factors, gamma = gamma,
      variance_estimator = "residual_df"),
    `CCY-BH` = function() ccy_ll_bh(Y, factors, gamma = gamma,
      bandwidth = ccy_bandwidth))
  out <- lapply(methods, function(name) {
    fit <- fits[[name]]()
    fit$implementation_label <- fit$implementation_label %||% fit$method
    fit$method <- name
    fit
  })
  setNames(out, methods)
}
