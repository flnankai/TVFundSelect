library(TVFundSelect)
named_mask <- stats::setNames(c(TRUE, TRUE, FALSE, TRUE), paste0("V", 1:4))
stopifnot(identical(TVFundSelect:::.cc_longest_span(named_mask),
                    c(first = 1L, last = 2L, T = 2L)))
near_cc <- function(a, b, tolerance = 2e-9) {
  isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tolerance))
}
fail_cc <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

# Independent literal regression reference: no package QA, geometry, state,
# correlation or bootstrap helper is used in this calculation.
cc_reference <- function(Y, factors, L) {
  TT <- ncol(Y); N <- nrow(Y); k <- ncol(factors) + 1L
  knots <- if (L > 3L) seq(0, 1, length.out = L - 1L)[-c(1L, L - 1L)] else numeric()
  basis <- splines::bs(seq_len(TT)/TT, knots = knots, degree = 2L,
                       Boundary.knots = c(0, 1), intercept = TRUE)
  X <- cbind(1, factors)
  D <- do.call(cbind, lapply(seq_len(k), function(j) basis * X[, j]))
  restricted <- lm.fit(X, t(Y)); unrestricted <- lm.fit(D, t(Y))
  rss0 <- colSums(restricted$residuals^2); rss1 <- colSums(unrestricted$residuals^2)
  C <- (rss0 - rss1) / rss0
  q <- k * (L - 1L); m <- TT - k; nu <- TT - k * L
  Z <- (m * C - q) / sqrt(2 * q * nu / (m + 2))
  residual <- sweep(unrestricted$residuals, 2L, colMeans(unrestricted$residuals), "-")
  corr <- stats::cor(residual)
  B2 <- max(N, nu/(nu - 1) * (sum(corr^2) - N^2/nu))
  B4 <- N
  for (i in seq_len(N)) for (j in seq_len(N)) {
    if (i != j && abs(corr[i, j]) > 2 * sqrt(log(2*N)/nu)) B4 <- B4 + corr[i, j]^4
  }
  V4 <- (8/q)*B2 + (2 + 4/q)*B4
  list(C = C, Z = Z, B2 = B2, B4 = B4, V4 = V4,
       statistic = c(S1 = sum(Z)/sqrt(B2), S2 = (sum(Z^2) - N)/sqrt(V4), MAX = max(Z)),
       fitted_null = restricted$fitted.values, residual = residual)
}

set.seed(6441)
N <- 8L; TT <- 64L
F <- matrix(rnorm(TT * 2L), TT, 2L)
Y <- matrix(rnorm(N * TT), N, TT) +
  outer(seq(-0.2, 0.3, length.out = N), rep(1, TT)) +
  outer(rep(1, N), 0.4 * F[, 1L] + 0.1 * sin(seq_len(TT)/TT * 2*pi))
L <- 4L; B <- 19L
saved_seed_cc <- .Random.seed; saved_kind_cc <- RNGkind()
thread_vars_cc <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
saved_threads_cc <- Sys.getenv(thread_vars_cc, unset = NA_character_)
fit <- constancy_test(Y, F, L = L, B = B, seed = 9331L, cores = 1L)
ref <- cc_reference(Y, F, L)
stopifnot(near_cc(fit$statistics[c("S1", "S2", "MAX")], ref$statistic),
          near_cc(fit$B2, ref$B2), near_cc(fit$B4, ref$B4), near_cc(fit$V4, ref$V4),
          identical(.Random.seed, saved_seed_cc), identical(RNGkind(), saved_kind_cc),
          identical(Sys.getenv(thread_vars_cc, unset = NA_character_), saved_threads_cc))

# Generate the actual null response and fully refit both regressions per draw.
# In particular residual correlations and normalization are not held fixed.
boot_reference <- t(vapply(seq_len(B), function(b) {
  pseudo <- t(ref$fitted_null + ref$residual[fit$bootstrap_indices[, b], , drop = FALSE])
  cc_reference(pseudo, F, L)$statistic[c("S2", "MAX")]
}, numeric(2L)))
stopifnot(near_cc(fit$bootstrap_statistics, boot_reference))
hits_ref <- colSums(sweep(boot_reference, 2L, ref$statistic[c("S2", "MAX")], ">="))
p_ref <- c(S1 = pnorm(ref$statistic["S1"], lower.tail = FALSE), (1 + hits_ref)/(B + 1))
terms <- ifelse(p_ref <= 0.5, 1/tan(pi * p_ref), -1/tan(pi * (1 - p_ref)))
terms[p_ref == 0.5] <- 0
cauchy_stat <- if (any(p_ref == 1)) -Inf else if (any(p_ref == 0)) Inf else mean(terms)
stopifnot(near_cc(fit$p_values[1:3], p_ref),
          near_cc(fit$p_values["Cauchy"], pcauchy(cauchy_stat, lower.tail = FALSE)),
          identical(as.numeric(fit$exceedances), as.numeric(hits_ref)))

# The individual Ma reference differs from the finite-restriction panel score.
single <- individual_constancy_test(Y[1L, ], F, L = L, B = B, seed = 9332L)
individual_Z <- ((TT - ncol(F) - 1L) * ref$C[1L] - (ncol(F) + 1L)*(L - 1L)) /
  sqrt(2 * (ncol(F) + 1L)*(L - 1L))
stopifnot(near_cc(single$C, ref$C[1L]), near_cc(single$statistic, individual_Z),
          near_cc(single$p_normal, pnorm(individual_Z, lower.tail = FALSE)))
old_kind_cc <- RNGkind()
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(9332L)
indices <- matrix(sample.int(TT, TT * B, replace = TRUE), TT, B)
do.call(RNGkind, as.list(old_kind_cc))
assign(".Random.seed", saved_seed_cc, envir = .GlobalEnv)
individual_boot <- vapply(seq_len(B), function(b) {
  pseudo <- ref$fitted_null[, 1L] + ref$residual[indices[, b], 1L]
  cc_reference(rbind(pseudo, pseudo), F, L)$C[1L]
}, numeric(1L))
stopifnot(near_cc(single$bootstrap_statistics, individual_boot),
          near_cc(single$p_bootstrap, (1 + sum(individual_boot >= single$C))/(B + 1)),
          identical(.Random.seed, saved_seed_cc))

# Underflow is explicitly distinguishable from a mathematical p-value of zero.
set.seed(6442)
strong_T <- 200L; strong_F <- matrix(rnorm(strong_T), strong_T, 1L)
strong_y <- 50 * (seq_len(strong_T)/strong_T)^2 + rnorm(strong_T, sd = 0.01)
strong <- individual_constancy_test(strong_y, strong_F, L = 3L, B = 3L, seed = 9333L)
stopifnot(strong$p_normal == 0, is.finite(strong$log_p_normal), strong$log_p_normal < -700)
longest_cc <- getFromNamespace(".cc_longest_span", "TVFundSelect")
stopifnot(identical(unname(longest_cc(c(TRUE, TRUE, FALSE, TRUE, TRUE))), c(1L, 2L, 2L)),
          longest_cc(rep(FALSE, 8L))["T"] == 0L)
fails <- list(
  quote(constancy_test(Y, F, L = 30L, B = 1L, cores = 1L)),
  quote(constancy_test(Y, cbind(F[, 1L], F[, 1L]), L = 3L, B = 1L, cores = 1L)),
  quote(constancy_test(Y, F, L = 2L, B = 1L, cores = 1L)),
  quote(constancy_test(Y, F, B = 0L, cores = 1L)),
  quote(constancy_test(Y, F, B = 1L, seed = Inf, cores = 1L)),
  quote(constancy_test(Y[, -1L], F, B = 1L, cores = 1L)),
  quote(individual_constancy_test(c(NA_real_, Y[1L, -1L]), F, B = 1L)),
  quote(individual_constancy_test(matrix(Y[1L, ], 1L), F, B = 1L))
)
for (expr in fails) fail_cc(eval(expr))

# Opt-in only: ordinary R CMD check does not spawn 25 local processes.
if (identical(Sys.getenv("TVFUNDSELECT_PARALLEL_TESTS"), "true")) {
  seed_before_cc_parallel <- .Random.seed; kind_before_cc_parallel <- RNGkind()
  serial <- constancy_test(Y, F, L = L, B = 25L, seed = 9334L, cores = 1L)
  parallel_fit <- constancy_test(Y, F, L = L, B = 25L, seed = 9334L, cores = 25L)
  stopifnot(identical(serial$bootstrap_indices, parallel_fit$bootstrap_indices),
            identical(serial$bootstrap_statistics, parallel_fit$bootstrap_statistics),
            identical(serial$p_values, parallel_fit$p_values),
            identical(serial$statistics, parallel_fit$statistics),
            identical(.Random.seed, seed_before_cc_parallel),
            identical(RNGkind(), kind_before_cc_parallel),
            identical(Sys.getenv(thread_vars_cc, unset = NA_character_), saved_threads_cc))
}
cat("Constancy: independent full-QR null refits, recomputed bootstrap normalization, tails and validation passed.\n")
