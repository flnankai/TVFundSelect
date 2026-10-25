# Appendix B: base-R, fixed-design QR implementation of the archived diagnostics.
# No native compiler, local cache or original result file is needed.
.cc_geometry <- function(factors, L) {
  TT <- nrow(factors); k <- ncol(factors) + 1L
  L <- .tf_integer(L, "L", 3L)
  knots <- if (L > 3L) seq(0, 1, length.out = L - 1L)[-c(1L, L - 1L)] else numeric()
  B <- splines::bs(seq_len(TT)/TT, knots = knots, degree = 2L,
                   Boundary.knots = c(0, 1), intercept = TRUE)
  X <- cbind(1, factors)
  D <- do.call(cbind, lapply(seq_len(k), function(j) B * X[, j]))
  scaled_qr <- function(x) {
    sc <- sqrt(colSums(x^2))
    if (any(!is.finite(sc)) || any(sc <= 0)) stop("Degenerate constancy design.")
    qr(sweep(x, 2L, sc, "/"), tol = 1e-10, LAPACK = FALSE)
  }
  qrX <- scaled_qr(X); qrD <- scaled_qr(D)
  if (qrX$rank != k || qrD$rank != k * L || TT - k * L <= 1L) {
    stop("The unrestricted constancy design has insufficient rank or residual degrees of freedom.")
  }
  QD <- qr.Q(qrD); QX <- qr.Q(qrX)
  completion <- qr.Q(qr(crossprod(QD, QX)), complete = TRUE)
  QA <- QD %*% completion[, seq.int(k + 1L, k * L), drop = FALSE]
  list(qrX = qrX, qrD = qrD, QA = QA, q = k*(L-1L), m = TT-k,
       nu = TT-k*L, T = TT, L = L, k = k)
}

.cc_cauchy <- function(p) {
  # Preserve the archived endpoint rule; retain log S1 separately for underflow.
  if (any(p == 1)) stat <- -Inf else if (any(p == 0)) stat <- Inf else {
    terms <- ifelse(p <= 0.5, 1/tan(pi*p), -1/tan(pi*(1-p)))
    terms[p == 0.5] <- 0
    stat <- mean(terms)
  }
  c(statistic = stat, p_value = stats::pcauchy(stat, lower.tail = FALSE))
}

.cc_state <- function(Yt, geometry) {
  U <- colSums(crossprod(geometry$QA, Yt)^2)
  residual <- qr.resid(geometry$qrD, Yt)
  W <- colSums(residual^2)
  if (any(!is.finite(U)) || any(!is.finite(W)) || any(W <= 0)) {
    stop("Nonfinite or degenerate constancy residuals.")
  }
  residual <- sweep(residual, 2L, colMeans(residual), "-")
  covariance <- crossprod(residual)/geometry$nu
  scales <- sqrt(diag(covariance))
  if (any(!is.finite(scales)) || any(scales <= 0)) stop("Zero constancy residual scale.")
  R <- covariance/outer(scales, scales)
  R[] <- pmax(-1, pmin(1, R)); diag(R) <- 1
  N <- ncol(Yt); q <- geometry$q; m <- geometry$m; nu <- geometry$nu
  B2 <- max(N, nu/(nu-1)*(sum(R^2)-N^2/nu))
  off <- abs(R[upper.tri(R)])
  B4 <- N + 2*sum(off[off > 2*sqrt(log(2*N)/nu)]^4)
  V4 <- (8/q)*B2 + (2+4/q)*B4
  C <- U/(U+W)
  Z <- (m*C-q)/sqrt(2*q*nu/(m+2))
  statistic <- c(S1 = sum(Z)/sqrt(B2), S2 = (sum(Z^2)-N)/sqrt(V4), MAX = max(Z))
  list(statistic = statistic, C = C, Z = Z, B2 = B2, B4 = B4,
       V4 = V4, residual = residual, geometry = geometry)
}

.cc_check_panel <- function(Y, factors) {
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 1L || ncol(Y) < 6L ||
      any(!is.finite(Y)) || !is.matrix(factors) || !is.numeric(factors) ||
      nrow(factors) != ncol(Y) || ncol(factors) < 1L || any(!is.finite(factors))) {
    stop("Y must be a finite N-by-T matrix and factors a finite T-by-p matrix.")
  }
}

.cc_boot_draw <- function(b, residual, geometry, indices) {
  # Adding the fitted constant null is annihilated by both fixed-design
  # residual projections. This is algebraically identical to full refitting;
  # correlation normalizations ARE recomputed in every generated panel.
  state <- .cc_state(residual[indices[, b], , drop = FALSE], geometry)
  unname(state$statistic[c("S2", "MAX")])
}

constancy_test <- function(Y, factors, L = NULL, B = 1999L, seed = 20281020L,
                           cores = 25L) {
  .cc_check_panel(Y, factors)
  B <- .tf_integer(B, "B"); cores <- .tf_integer(cores, "cores")
  if (is.null(L)) L <- 3L + floor(ncol(Y)^0.2)
  geometry <- .cc_geometry(factors, L)
  state <- .cc_state(t(Y), geometry)
  indices <- .tf_with_seed(seed, function() matrix(
    sample.int(ncol(Y), ncol(Y)*as.double(B), replace = TRUE), ncol(Y), B))
  parts <- .tf_map(seq_len(B), .cc_boot_draw, cores = cores,
    residual = state$residual, geometry = geometry, indices = indices)
  bootstrap <- matrix(unlist(parts), nrow = B, ncol = 2L, byrow = TRUE,
                      dimnames = list(NULL, c("S2", "MAX")))
  hits <- colSums(sweep(bootstrap, 2L, state$statistic[c("S2", "MAX")], ">="))
  pv <- c(S1 = stats::pnorm(state$statistic["S1"], lower.tail = FALSE),
          (1+hits)/(B+1))
  names(pv)[1L] <- "S1"
  combined <- .cc_cauchy(pv)
  out <- list(statistics = c(state$statistic, Cauchy = unname(combined["statistic"])),
    p_values = c(pv, Cauchy = unname(combined["p_value"])),
    log_p_S1 = unname(stats::pnorm(state$statistic["S1"], lower.tail = FALSE, log.p = TRUE)),
    bootstrap_statistics = bootstrap, bootstrap_indices = indices,
    exceedances = hits, B = B, seed = seed, cores = cores, L = L,
    N = nrow(Y), T = ncol(Y), q = geometry$q, nu = geometry$nu,
    B2 = state$B2, B4 = state$B4, V4 = state$V4,
    null = "All alpha and factor-loading functions are constant",
    calibration = "S1 normal; S2/MAX iid whole-date residual bootstrap; equal-weight Cauchy",
    computation = "new draws; exact fixed-design QR coefficient refits; correlations recomputed in every draw")
  out
}

individual_constancy_test <- function(y, factors, L = NULL, B = 1999L,
                                      seed = 20261020L) {
  if (!is.numeric(y) || !is.null(dim(y))) stop("y must be a numeric time-series vector.")
  .cc_check_panel(matrix(y, nrow = 1L), factors)
  if (is.null(L)) L <- 3L + floor(length(y)^0.2)
  B <- .tf_integer(B, "B")
  des <- .cc_geometry(factors, L)
  U <- sum(crossprod(des$QA, y)^2)
  residual <- drop(qr.resid(des$qrD, y))
  W <- sum(residual^2)
  if (!is.finite(W) || W <= 0) stop("Degenerate individual constancy residuals.")
  C <- U/(U+W); Z <- (des$m*C-des$q)/sqrt(2*des$q)
  residual <- residual-mean(residual)
  indices <- .tf_with_seed(seed, function() matrix(
    sample.int(length(y), length(y)*as.double(B), replace = TRUE), length(y), B))
  pseudo <- matrix(residual[indices], length(y), B)
  Us <- colSums(crossprod(des$QA, pseudo)^2)
  Ws <- colSums(qr.resid(des$qrD, pseudo)^2)
  Cs <- Us/(Us+Ws)
  if (any(!is.finite(Cs)) || any(Ws <= 0)) stop("Degenerate individual bootstrap draw.")
  hits <- sum(Cs >= C)
  list(C = C, statistic = Z, p_normal = stats::pnorm(Z, lower.tail = FALSE),
       log_p_normal = stats::pnorm(Z, lower.tail = FALSE, log.p = TRUE),
       p_bootstrap = (1+hits)/(B+1), exceedances = hits, bootstrap_statistics = Cs,
       B = B, seed = seed, T = length(y), L = L, q = des$q, nu = des$nu,
       null = "The alpha and all factor loadings are jointly constant")
}

.cc_longest_span <- function(observed) {
  observed <- unname(observed)
  z <- rle(observed); ends <- cumsum(z$lengths); starts <- ends-z$lengths+1L
  keep <- which(z$values)
  if (!length(keep)) return(c(first = NA_integer_, last = NA_integer_, T = 0L))
  at <- keep[which.max(z$lengths[keep])]
  c(first = starts[at], last = ends[at], T = z$lengths[at])
}

.cc_individual_job <- function(job, data, B) {
  ii <- seq.int(job$first, job$last)
  y <- data$returns[job$fund, ii]-data$rf[ii]
  fit <- individual_constancy_test(y, data$factors[ii, , drop = FALSE], B = B, seed = job$seed)
  data.frame(sample = job$sample, fund_id = data$fund_id[job$fund],
    fund_row = job$fund, start = as.character(data$dates[job$first]),
    end = as.character(data$dates[job$last]), T = fit$T, L = fit$L, nu = fit$nu,
    C = fit$C, Z_Ma = fit$statistic, p_normal = fit$p_normal,
    log_p_normal = fit$log_p_normal, p_bootstrap = fit$p_bootstrap,
    B = B, seed = job$seed, status = "computed")
}

reproduce_constancy <- function(data, B = 1999L, cores = 25L, seed = 20261019L,
    panel_lengths = c(ncol(data$returns), 180L, 120L), common_length = 180L,
    min_residual_df = 30L, fund_rows = seq_len(nrow(data$returns)), out_dir = NULL) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  .emp_validate_data(data)
  B <- .tf_integer(B, "B"); cores <- .tf_integer(cores, "cores")
  seed <- .tf_integer(seed, "seed", 0L)
  min_residual_df <- .tf_integer(min_residual_df, "min_residual_df", 2L)
  common_length <- .tf_integer(common_length, "common_length", 6L)
  TT <- ncol(data$returns)
  if (common_length > TT || any(panel_lengths > TT) || any(panel_lengths < 6L) ||
      any(!is.finite(panel_lengths)) || any(panel_lengths != floor(panel_lengths)) ||
      any(!is.finite(fund_rows)) || any(fund_rows != floor(fund_rows)) ||
      any(fund_rows < 1L | fund_rows > nrow(data$returns)) || anyDuplicated(fund_rows)) {
    stop("Invalid panel lengths or fund_rows.")
  }
  jobs <- list(); unavailable <- list(); audit <- list()
  k <- ncol(data$factors)+1L
  for (i in fund_rows) {
    span <- .cc_longest_span(data$observed[i, ])
    Li <- 3L+floor(span["T"]^0.2)
    nu <- unname(span["T"]-k*Li)
    eligible <- nu >= min_residual_df
    audit[[length(audit)+1L]] <- data.frame(fund_row = i, fund_id = data$fund_id[i],
      first = span["first"], last = span["last"], T = span["T"], L = Li, nu = nu,
      eligible = eligible)
    if (eligible) jobs[[length(jobs)+1L]] <- list(sample = "longest_history", fund = i,
      first = span["first"], last = span["last"], seed = seed+i)
    else unavailable[[length(unavailable)+1L]] <- data.frame(fund_row = i,
      fund_id = data$fund_id[i], T = span["T"], nu = nu, status = "insufficient residual df")
  }
  rows <- seq.int(TT-common_length+1L, TT)
  complete <- fund_rows[rowSums(!data$observed[fund_rows, rows, drop = FALSE]) == 0L]
  for (i in complete) jobs[[length(jobs)+1L]] <- list(sample = paste0("last_", common_length),
    fund = i, first = min(rows), last = max(rows), seed = seed+10000L+i)
  individual <- if (length(jobs)) do.call(rbind, .tf_map(jobs, .cc_individual_job,
    cores = cores, data = data, B = B)) else data.frame()
  summary <- list()
  if (nrow(individual)) {
    individual$p_normal_BH <- individual$p_bootstrap_BH <- NA_real_
    for (s in unique(individual$sample)) {
      ix <- which(individual$sample == s)
      individual$p_normal_BH[ix] <- stats::p.adjust(individual$p_normal[ix], "BH")
      individual$p_bootstrap_BH[ix] <- stats::p.adjust(individual$p_bootstrap[ix], "BH")
      for (cal in c("normal", "bootstrap")) {
        p <- individual[[paste0("p_", cal)]][ix]
        summary[[length(summary)+1L]] <- data.frame(sample = s, calibration = cal,
          N = length(ix), reject_05 = sum(p <= .05), reject_10 = sum(p <= .1),
          reject_BH10 = sum(stats::p.adjust(p, "BH") <= .1))
      }
    }
  }
  panels <- list(); panel_rows <- list()
  for (j in seq_along(unique(panel_lengths))) {
    len <- unique(panel_lengths)[j]; rows <- seq.int(TT-len+1L, TT)
    ids <- fund_rows[rowSums(!data$observed[fund_rows, rows, drop = FALSE]) == 0L]
    if (!length(ids)) stop("No complete funds in the panel of length ", len, ".")
    Y <- sweep(data$returns[ids, rows, drop = FALSE], 2L, data$rf[rows], "-")
    fit <- constancy_test(Y, data$factors[rows, , drop = FALSE], B = B,
                          cores = cores, seed = seed+20000L+j)
    label <- paste0("last_", len); panels[[label]] <- fit
    panel_rows[[j]] <- data.frame(sample = label, N = length(ids), T = len,
      method = names(fit$p_values), statistic = unname(fit$statistics),
      p_value = unname(fit$p_values), log_p_S1 = fit$log_p_S1, B = B)
  }
  out <- list(individual = individual,
    individual_summary = if (length(summary)) do.call(rbind, summary) else data.frame(),
    panel = do.call(rbind, panel_rows), panels = panels,
    sample_audit = do.call(rbind, audit),
    unavailable = if (length(unavailable)) do.call(rbind, unavailable) else data.frame(),
    config = list(B = B, cores = cores, seed = seed, computation = "fresh from raw returns",
      min_residual_df = min_residual_df, common_length = common_length,
      panel_lengths = panel_lengths))
  if (!is.null(out_dir)) {
    .tf_save_tables(out, out_dir, "constancy")
    grDevices::pdf(file.path(out_dir, "Figure_7_constancy.pdf"), width = 9, height = 4)
    oldpar <- graphics::par(mfrow = c(1, 2))
    for (s in unique(individual$sample)) {
      graphics::hist(individual$p_normal[individual$sample == s], breaks = seq(0, 1, .05),
        main = s, xlab = "Coefficient-constancy p-value", col = "lightblue")
    }
    graphics::par(oldpar); grDevices::dev.off()
  }
  out
}
