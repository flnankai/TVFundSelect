# Shared validation. Core numerical operations retain the archived arithmetic.
.fund_scalar <- function(x, name, lower, upper = Inf, integer = FALSE,
                         closed_lower = TRUE, closed_upper = TRUE) {
  invalid <- !is.numeric(x) || length(x) != 1L || !is.finite(x)
  if (!invalid) {
    invalid <- (if (closed_lower) x < lower else x <= lower) ||
      (if (closed_upper) x > upper else x >= upper) ||
      (integer && x != floor(x))
  }
  if (invalid) stop("Invalid ", name, ".", call. = FALSE)
  invisible(TRUE)
}

.fund_panel <- function(Y, factors, gamma, min_funds = 1L, min_time = 5L) {
  assert_finite_matrix(Y, "Y")
  assert_finite_matrix(factors, "factors")
  if (nrow(Y) < min_funds || ncol(Y) < min_time || nrow(factors) != ncol(Y)) {
    stop("Y must be N by T and factors T by p with compatible nonempty dimensions.",
         call. = FALSE)
  }
  .fund_scalar(gamma, "gamma (must lie in (0, 0.5))", 0, 0.5,
               closed_lower = FALSE, closed_upper = FALSE)
  invisible(TRUE)
}

.fund_knots <- function(knot_candidates, degree) {
  .fund_scalar(degree, "degree", 1, .Machine$integer.max, integer = TRUE)
  if (!is.numeric(knot_candidates) || !length(knot_candidates) ||
      any(!is.finite(knot_candidates)) || any(knot_candidates < 0) ||
      any(knot_candidates != floor(knot_candidates)) ||
      any(knot_candidates > .Machine$integer.max) || anyDuplicated(knot_candidates)) {
    stop("knot_candidates must contain distinct nonnegative integers.", call. = FALSE)
  }
  invisible(TRUE)
}
