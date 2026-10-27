library(TVFundSelect)
near <- function(x, y, tolerance = 1e-12) isTRUE(all.equal(x, y, tolerance = tolerance))
must_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
output_files <- function(path, ...) {
  restore <- getFromNamespace(".tf_utf8_scope", "TVFundSelect")()
  on.exit(restore(), add = TRUE)
  list.files(path, ...)
}
remove_temporary_output <- function(path) {
  restore <- getFromNamespace(".tf_utf8_scope", "TVFundSelect")()
  on.exit(restore(), add = TRUE)
  resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
  temporary_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  stopifnot(startsWith(resolved, paste0(temporary_root, "/")))
  unlink(path, recursive = TRUE)
}
set.seed(20261025)
saved_seed <- .Random.seed
saved_kind <- RNGkind()
panel <- simulate_fund_panel(N = 40L, TT = 80L, example = 1L,
                            distribution = "t3", signal_shift = 0.1, seed = 99L)
again <- simulate_fund_panel(N = 40L, TT = 80L, example = 1L,
                            distribution = "t3", signal_shift = 0.1, seed = 99L)
stopifnot(identical(panel, again), identical(.Random.seed, saved_seed),
          identical(RNGkind(), saved_kind), identical(dim(panel$Y), c(40L, 80L)),
          max(abs(rowMeans(panel$alpha) - panel$delta)) < 1e-12,
          all(apply(panel$alpha, 1L, function(a) diff(range(a))) > 1e-8),
          length(panel$positive) == 4L, length(panel$negative) == 4L,
          length(intersect(panel$positive, panel$negative)) == 0L,
          all(panel$delta[panel$positive] > 0), all(panel$delta[panel$negative] < 0),
          panel$settings$positive_signal_scale == sqrt(200 / 80))

# An algebraic, not sample-moment, check of the heavy-tail covariance scaling.
error_generator <- getFromNamespace("simulate_reference_errors", "TVFundSelect")
set.seed(77)
t3 <- error_generator(4L, 15L, distribution = "t3")
set.seed(77)
gaussian <- matrix(rnorm(60L), 15L, 4L) %*% chol(toeplitz(0.5^(0:3)))
scale <- sqrt(rchisq(15L, df = 3L))
stopifnot(near(t3, t(gaussian / scale)))

set.seed(905)
before <- .Random.seed
settings <- list(reps = 2L, examples = 1L, times = 80L,
  distributions = c("normal", "t3"), signals = c(0.04, 0.10), N = 40L,
  seed = 40001L, plots = FALSE, chunk_size = 1L, ccy_B = 99L)
serial <- do.call(reproduce_simulations, c(settings, list(cores = 1L)))
stopifnot(identical(.Random.seed, before), nrow(serial$tasks) == 8L,
  nrow(serial$raw) == 40L, nrow(serial$summary) == 20L,
  all(serial$raw$success), all(serial$raw$FDP >= 0 & serial$raw$FDP <= 1),
  all(serial$raw$TDP >= 0 & serial$raw$TDP <= 1),
  !anyDuplicated(serial$tasks$task_seed), isTRUE(serial$config$fresh_computation))
for (i in seq_len(nrow(serial$summary))) {
  s <- serial$summary[i, ]
  r <- serial$raw[serial$raw$method == s$method &
    serial$raw$distribution == s$distribution & serial$raw$signal_shift == s$signal_shift, ]
  stopifnot(near(s$mean_FDP, mean(r$FDP)), near(s$mean_TDP, mean(r$TDP)),
    near(s$mcse_FDP, sd(r$FDP) / sqrt(nrow(r))))
}

# Independently regenerate one returned seed and recompute all five selections.
first <- serial$tasks[1L, ]
fresh <- simulate_fund_panel(N = 40L, TT = first$T, example = first$example,
  distribution = first$distribution, signal_shift = first$signal_shift,
  seed = first$task_seed)
fits <- fit_fund_methods(fresh$Y, fresh$factors, ccy_B = 99L,
  ccy_seed = serial$config$ccy_seed + 2029L * first$task_id, ccy_cores = 1L)
for (name in names(fits)) {
  selected <- fits[[name]]$selected
  correct <- sum(fresh$delta[selected] > 0)
  row <- serial$raw[serial$raw$task_id == 1L & serial$raw$method == name, ]
  stopifnot(near(row$FDP, (length(selected) - correct) / max(1L, length(selected))),
    near(row$TDP, correct / length(fresh$positive)))
}

bandwidth <- reproduce_bandwidth(reps = 2L, cores = 1L, examples = 1L,
  times = 80L, distributions = "normal", signals = 0.10, N = 40L,
  seed = 40002L, plots = FALSE, multipliers = c(1, 2, 4, 8, 16), chunk_size = 1L, ccy_B = 99L)
stopifnot(nrow(bandwidth$raw) == 10L, nrow(bandwidth$summary) == 5L,
  all(bandwidth$raw$bandwidth == bandwidth$raw$bandwidth_multiplier *
        1.06 * sd((1:80) / 80) * 80^(-1/5)))

check_output <- function() {
  out <- file.path(tempfile("fresh-simulation-"), "\u4e2d\u6587")
  on.exit(remove_temporary_output(out), add = TRUE)
  arguments <- list(reps = 1L, cores = 1L, examples = 1L, times = 80L,
    distributions = "normal", signals = 0.1, N = 20L,
    methods = c("TV-LS-BH", "CC-LS-BH"), output_dir = out, plots = TRUE)
  calculated <- do.call(reproduce_simulations, arguments)
  stopifnot(identical(sort(output_files(out)), c("fdp_tdp.pdf", "simulation_summary.txt")))
  must_error(do.call(reproduce_simulations, arguments))
  must_error(do.call(reproduce_simulations, c(arguments, list(resume = TRUE))))
  changed <- arguments
  changed$seed <- 99L
  must_error(do.call(reproduce_simulations, c(changed, list(resume = TRUE))))
  invisible(TRUE)
}
check_output()

if (identical(Sys.getenv("TVFUNDSELECT_PARALLEL_TESTS"), "true")) {
  parallel <- do.call(reproduce_simulations, c(settings, list(cores = 25L)))
  stopifnot(identical(serial$raw, parallel$raw), identical(serial$summary, parallel$summary),
    identical(serial$tasks, parallel$tasks), identical(.Random.seed, before))
  bandwidth25 <- reproduce_bandwidth(reps = 2L, cores = 25L, examples = 1L,
    times = 80L, distributions = "normal", signals = 0.10, N = 40L,
    seed = 40002L, plots = FALSE, multipliers = c(1, 2, 4, 8, 16), chunk_size = 1L, ccy_B = 99L)
  stopifnot(identical(bandwidth$raw, bandwidth25$raw),
    identical(bandwidth$summary, bandwidth25$summary))
  check_parallel_unicode <- function() {
    out <- file.path(tempfile("parallel-output-"), "\u91cd\u73b0\u7ed3\u679c")
    old_locale <- Sys.getlocale("LC_CTYPE")
    on.exit(remove_temporary_output(out), add = TRUE)
    p <- do.call(reproduce_simulations, c(settings, list(cores = 2L, output_dir = out)))
    stopifnot(identical(serial$raw, p$raw), identical(serial$summary, p$summary),
      identical(Sys.getlocale("LC_CTYPE"), old_locale),
      dir.exists(p$output_dir), file.exists(file.path(p$output_dir, "simulation_summary.txt")))
    must_error(do.call(reproduce_simulations,
      c(settings, list(cores = 2L, output_dir = out, resume = TRUE))))
    stopifnot(!length(output_files(out, pattern = "[.]rds$", recursive = TRUE)))
  }
  check_parallel_unicode()
}

must_error(simulate_fund_panel(N = 2))
must_error(simulate_fund_panel(example = 4))
must_error(simulate_fund_panel(rho = 1))
must_error(simulate_fund_panel(signal_shift = -0.02))
must_error(reproduce_simulations(reps = 0))
must_error(reproduce_simulations(cores = 0))
must_error(reproduce_simulations(examples = c(1, 1)))
must_error(reproduce_simulations(signals = c(0.1, 0.1)))
must_error(reproduce_simulations(seed = .Machine$integer.max))
must_error(reproduce_simulations(resume = TRUE))
cat("Fresh simulation, exact DGP scaling, summaries, restart guards and optional 25-worker checks passed.\n")
