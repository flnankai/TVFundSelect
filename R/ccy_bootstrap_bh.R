# CCY average-alpha bootstrap inference used by the public CCY-BH method.
# Intermediate bootstrap matrices are kept in memory only.
#
# With A the conditional-alpha smoother, a = colMeans(A), S the fitted-return
# smoother and delta_hat = Y a, the fitted null mean is S y_i - delta_hat_i 1.
# Its estimated conditional alpha has average zero, while time variation in
# alpha and all estimated factor loadings are retained. Reapplying the local
# estimator gives bootstrap mean a' S y_i - delta_hat_i, which estimates the
# smoothing bias. The primary p-value retains this bootstrap null mean. A
# separately labelled centered-pivot diagnostic removes it.
#
# Whole residual vectors are sampled across dates, using one common index for
# every fund. Optional common signs match the paper's empirical calibration
# design and impose symmetry. These iid schemes do not preserve serial error
# dependence or the date-specific heteroskedasticity of a fixed-date wild scheme.

.ccy_boot_integer <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > .Machine$integer.max) {
    stop(name, " must be one finite integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(x)
}

.ccy_boot_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

# The master draws all random indices/signs. This worker has no random-number
# generator calls, so scheduling and worker count cannot alter the experiment.
.ccy_boot_batch <- function(draw, residual_pool, a, bootstrap_bias, delta_hat,
                            resampling, keep_bootstrap, keep_resamples) {
  TT <- nrow(draw)
  count <- ncol(draw)
  if (resampling == "signed_residual") {
    indices <- (draw - 1L) %% TT + 1L
    signs <- matrix(ifelse(draw <= TT, 1L, -1L), TT, count)
  } else {
    indices <- draw
    signs <- matrix(1L, TT, count)
  }
  W <- matrix(0, TT, count)
  offsets <- TT * (seq_len(count) - 1L)
  for (t in seq_len(TT)) {
    pos <- indices[t, ] + offsets
    W[pos] <- W[pos] + a[t] * signs[t, ]
  }
  centered_draws <- residual_pool %*% W
  null_draws <- sweep(centered_draws, 1L, bootstrap_bias, "+")
  if (any(!is.finite(null_draws))) {
    stop("A bootstrap average-alpha statistic is nonfinite.", call. = FALSE)
  }
  list(raw_hits = rowSums(sweep(null_draws, 1L, delta_hat, ">=")),
       centered_hits = rowSums(sweep(centered_draws, 1L, delta_hat, ">=")),
       statistics = if (keep_bootstrap) null_draws else NULL,
       indices = if (keep_resamples) indices else NULL,
       signs = if (keep_resamples) signs else NULL,
       worker_pid = Sys.getpid())
}

# cores controls parallel bootstrap batches for ONE panel. If an outer Monte
# Carlo loop already uses 25 workers, call this function with cores=1 to avoid
# nested parallelism. Windows uses 25 independent PSOCK R worker processes.
ccy_bootstrap_bh <- function(Y, factors, gamma = 0.1, bandwidth = NULL,
    B = 9999L, seed = 20261027L,
    resampling = c("residual", "signed_residual"), df_adjust = TRUE,
    batch_size = 500L, keep_bootstrap = FALSE, keep_resamples = FALSE,
    cores = 25L) {
  if (!requireNamespace("TVFundSelect", quietly = TRUE)) {
    stop("Install TVFundSelect before sourcing and using this experimental function.",
         call. = FALSE)
  }
  ns <- asNamespace("TVFundSelect")
  get(".fund_panel", ns)(Y, factors, gamma)
  B <- .ccy_boot_integer(B, "B")
  seed <- .ccy_boot_integer(seed, "seed", 0L)
  batch_size <- .ccy_boot_integer(batch_size, "batch_size")
  cores <- .ccy_boot_integer(cores, "cores")
  cores_used <- min(cores, B)
  df_adjust <- .ccy_boot_flag(df_adjust, "df_adjust")
  keep_bootstrap <- .ccy_boot_flag(keep_bootstrap, "keep_bootstrap")
  keep_resamples <- .ccy_boot_flag(keep_resamples, "keep_resamples")
  resampling <- match.arg(resampling)
  N <- nrow(Y)
  TT <- ncol(Y)
  bandwidth_rule <- if (is.null(bandwidth)) "CCY-normal-reference" else "user-specified"
  if (is.null(bandwidth)) {
    bandwidth <- 1.06 * stats::sd(seq_len(TT) / TT) * TT^(-1 / 5)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("bandwidth must be one finite positive number.", call. = FALSE)
  }
  if (keep_bootstrap && as.double(N) * B > 25000000) {
    stop("Keeping bootstrap statistics would exceed 25 million cells; use keep_bootstrap=FALSE.",
         call. = FALSE)
  }
  if (keep_resamples && as.double(TT) * B > 25000000) {
    stop("Keeping bootstrap indices would exceed 25 million cells; use keep_resamples=FALSE.",
         call. = FALSE)
  }
  # Bound temporary matrices independently of the total bootstrap count.
  batch_used <- max(1L, min(batch_size, B,
                           floor(2000000 / (3 * TT + 2 * N))))
  if (cores_used > 1L) {
    # Bound the aggregate working matrices in one dispatch wave and ensure
    # enough batches to engage every requested worker when B >= cores.
    batch_used <- max(1L, min(batch_used, floor(B / cores_used),
                             floor(4000000 / (cores_used * (3 * TT + 2 * N)))))
  }
  design <- get("ccy_ll_design", ns)(factors, bandwidth,
                                     keep_alpha_smoother = keep_resamples)
  S <- design$fitted_smoother
  M <- diag(TT) - S
  a <- design$average_alpha_weights
  if (abs(sum(a) - 1) > 1e-7 || max(abs(S %*% rep(1, TT) - 1)) > 1e-7) {
    stop("Local-linear constant-reproduction check failed.", call. = FALSE)
  }
  fitted <- Y %*% t(S)
  residuals <- Y - fitted
  delta_hat <- drop(Y %*% a)
  null_mean <- sweep(fitted, 1L, delta_hat, "-")
  bootstrap_bias <- drop(null_mean %*% a)

  # S is generally neither symmetric nor idempotent. After centering residuals,
  # their iid-homoskedastic noise degrees of freedom are ||C(I-S)||_F^2,
  # C=I-11'/T. T-tr(S) is not the residual degrees of freedom for this smoother.
  residual_df_uncentered <- sum(M * M)
  residual_df <- residual_df_uncentered - sum(colSums(M)^2) / TT
  if (!is.finite(residual_df) || residual_df <= 1e-10) {
    stop("The centered local-linear residual degrees of freedom are nonpositive.",
         call. = FALSE)
  }
  residual_scale <- if (df_adjust) sqrt(TT / residual_df) else 1
  residual_pool <- sweep(residuals, 1L, rowMeans(residuals), "-") * residual_scale
  if (any(!is.finite(residual_pool)) || any(rowSums(residual_pool^2) <= 0)) {
    stop("A centered residual pool is zero or nonfinite.", call. = FALSE)
  }
  # The original normal-reference result is returned solely as a comparator.
  # Its denominator is not used to obtain either bootstrap p-value.
  normal_se <- sqrt(drop((residuals^2) %*% design$variance_weights) / TT)
  if (any(!is.finite(normal_se)) || any(normal_se <= 0)) {
    stop("An original CCY standard error is zero or nonfinite.", call. = FALSE)
  }
  normal_statistics <- delta_hat / normal_se
  normal_p_values <- stats::pnorm(normal_statistics, lower.tail = FALSE)
  normal_bh <- get("bh_select", ns)(normal_p_values, gamma)
  null_alpha <- if (keep_resamples) {
    sweep(Y %*% t(design$alpha_smoother), 1L, delta_hat, "-")
  } else NULL
  if (keep_resamples && max(abs(rowMeans(null_alpha))) >
      1e-8 * max(1, max(abs(Y)))) {
    stop("The fitted null alpha does not have average zero.", call. = FALSE)
  }

  # Save the caller's complete RNG state. No bootstrap randomness is consumed
  # outside this scope, and an error restores state as well as normal return.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()
  caps <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  old_caps <- Sys.getenv(caps, unset = NA_character_)
  cl <- NULL
  on.exit({
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    absent <- is.na(old_caps)
    if (any(absent)) Sys.unsetenv(caps[absent])
    if (any(!absent)) do.call(Sys.setenv, as.list(old_caps[!absent]))
  }, add = TRUE)
  model <- list(residual_pool = residual_pool, a = a,
    bootstrap_bias = bootstrap_bias, delta_hat = delta_hat,
    resampling = resampling, keep_bootstrap = keep_bootstrap,
    keep_resamples = keep_resamples)
  worker_pids <- Sys.getpid()
  if (cores_used > 1L) {
    do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(caps)), caps)))
    cl <- parallel::makePSOCKcluster(cores_used, master = "127.0.0.1",
                                    setup_timeout = 30)
    initialize <- function(model, batch_function, caps, ctype) {
      if (.Platform$OS.type == "windows") {
        changed <- suppressWarnings(Sys.setlocale("LC_CTYPE", ctype))
        if (!nzchar(changed)) stop("Could not initialize worker locale.")
      }
      do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(caps)), caps)))
      assign(".ccy_boot_model", model, envir = .GlobalEnv)
      environment(batch_function) <- baseenv()
      assign(".ccy_boot_batch", batch_function, envir = .GlobalEnv)
      Sys.getpid()
    }
    environment(initialize) <- baseenv()
    worker_pids <- unlist(parallel::clusterCall(cl, initialize, model,
      .ccy_boot_batch, caps, Sys.getlocale("LC_CTYPE")), use.names = FALSE)
    worker <- function(draw) {
      do.call(get(".ccy_boot_batch", envir = .GlobalEnv),
              c(list(draw = draw), get(".ccy_boot_model", envir = .GlobalEnv)))
    }
    environment(worker) <- baseenv()
  }
  # Seed only after worker startup; startup is allowed to use randomness without
  # affecting the bootstrap stream, and on.exit restores the caller's RNG.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  raw_hits <- centered_hits <- numeric(N)
  bootstrap_statistics <- if (keep_bootstrap) matrix(NA_real_, N, B) else NULL
  bootstrap_indices <- if (keep_resamples) matrix(NA_integer_, TT, B) else NULL
  bootstrap_signs <- if (keep_resamples) matrix(NA_integer_, TT, B) else NULL
  starts <- seq.int(1L, B, by = batch_used)
  used_worker_pids <- integer()
  for (wave_first in seq.int(1L, length(starts), by = cores_used)) {
    wave <- seq.int(wave_first, min(length(starts), wave_first + cores_used - 1L))
    # One integer represents both index and sign. All random draws originate on
    # the master in bootstrap order, independently of chunking and scheduling.
    jobs <- lapply(starts[wave], function(first) {
      count <- min(batch_used, B - first + 1L)
      matrix(sample.int(if (resampling == "signed_residual") 2L * TT else TT,
                        TT * count, replace = TRUE), TT, count)
    })
    results <- if (cores_used == 1L) {
      lapply(jobs, function(draw) do.call(.ccy_boot_batch, c(list(draw = draw), model)))
    } else parallel::parLapply(cl, jobs, worker)
    for (j in seq_along(results)) {
      answer <- results[[j]]
      first <- starts[wave[j]]
      at <- seq.int(first, length.out = ncol(jobs[[j]]))
      raw_hits <- raw_hits + answer$raw_hits
      centered_hits <- centered_hits + answer$centered_hits
      used_worker_pids <- union(used_worker_pids, answer$worker_pid)
      if (keep_bootstrap) bootstrap_statistics[, at] <- answer$statistics
      if (keep_resamples) {
        bootstrap_indices[, at] <- answer$indices
        bootstrap_signs[, at] <- answer$signs
      }
    }
  }
  p_values <- (1 + raw_hits) / (B + 1)
  centered_p_values <- (1 + centered_hits) / (B + 1)
  bh <- get("bh_select", ns)(p_values, gamma)
  centered_bh <- get("bh_select", ns)(centered_p_values, gamma)
  list(method = "CCY-Boot-BH", selected = bh$selected, k = bh$k,
       bh_threshold = bh$threshold, p_values = p_values,
       statistics = delta_hat, delta_hat = delta_hat,
       centered_method = "CCY-Boot-Centered-BH",
       centered_p_values = centered_p_values,
       centered_selected = centered_bh$selected, centered_k = centered_bh$k,
       centered_bh_threshold = centered_bh$threshold,
       exceedances = raw_hits, centered_exceedances = centered_hits,
       bootstrap_bias = bootstrap_bias, bootstrap_null_mean = bootstrap_bias,
       standard_errors = normal_se, normal_statistics = normal_statistics,
       normal_p_values = normal_p_values, normal_selected = normal_bh$selected,
       normal_bh_threshold = normal_bh$threshold,
       residual_df = residual_df, residual_df_uncentered = residual_df_uncentered,
       residual_scale = residual_scale, smoother_trace = sum(diag(S)),
       bandwidth = bandwidth, bandwidth_rule = bandwidth_rule,
       average_alpha_weights = a, residual_variance_weights = design$variance_weights,
       B = B, seed = seed, gamma = gamma, resampling = resampling,
       df_adjust = df_adjust, batch_size = batch_used,
       cores_requested = cores, cores_used = cores_used,
       worker_pids = worker_pids, computing_worker_pids = used_worker_pids,
       p_value_resolution = 1 / (B + 1),
       first_BH_threshold = gamma / N,
       minimum_BH_discovery_count = ceiling(N / (gamma * (B + 1))),
       bootstrap_statistics = bootstrap_statistics,
       bootstrap_indices = bootstrap_indices, bootstrap_signs = bootstrap_signs,
       fitted_null = if (keep_resamples) null_mean else NULL,
       alpha_null = null_alpha,
       residual_pool = if (keep_resamples) residual_pool else NULL,
       original_residuals = if (keep_resamples) residuals else NULL,
       package_version = as.character(utils::packageVersion("TVFundSelect")),
       interpretation = paste("Fixed-factor, fixed-bandwidth average-alpha null bootstrap.",
         "Primary raw-null draws retain estimated smoothing bias; centered draws are a diagnostic.",
         "Plus-one tails are finite Monte Carlo estimates, not a proof of finite-sample BH FDR control."))
}
