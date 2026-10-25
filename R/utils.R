# Transactional, deterministic utilities shared by fresh-computation workflows.
.tf_utf8_scope <- function() {
  old <- Sys.getlocale("LC_CTYPE")
  if (.Platform$OS.type == "windows" && !grepl("utf-?8", old, ignore.case = TRUE)) {
    chosen <- suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
    if (!nzchar(chosen)) chosen <- suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
    if (!nzchar(chosen)) stop("A UTF-8 character locale is needed for Windows file workflows.")
  }
  function() {
    if (!identical(Sys.getlocale("LC_CTYPE"), old)) suppressWarnings(Sys.setlocale("LC_CTYPE", old))
    invisible(NULL)
  }
}

.tf_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(name, " must be a finite integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.tf_with_seed <- function(seed, FUN) {
  seed <- .tf_integer(seed, "seed", 0L)
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  FUN()
}

.tf_map <- function(X, FUN, cores = 25L, ...) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  cores <- .tf_integer(cores, "cores")
  vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  old <- Sys.getenv(vars, unset = NA_character_)
  on.exit({
    absent <- is.na(old)
    if (any(absent)) Sys.unsetenv(vars[absent])
    if (any(!absent)) do.call(Sys.setenv, as.list(old[!absent]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(vars)), vars)))
  if (!length(X)) return(list())
  if (cores == 1L || length(X) == 1L) return(lapply(X, FUN, ...))
  cl <- parallel::makePSOCKcluster(min(cores, length(X)), setup_timeout = 30)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  initialize <- function(paths, ctype) {
    if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ctype))
    .libPaths(paths)
    loadNamespace("TVFundSelect")
    NULL
  }
  environment(initialize) <- baseenv()
  parallel::clusterCall(cl, initialize, .libPaths(), Sys.getlocale("LC_CTYPE"))
  parallel::parLapply(cl, X, FUN, ...)
}

.tf_output_dir <- function(out_dir) {
  if (is.null(out_dir)) return(invisible(NULL))
  if (!is.character(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)) {
    stop("out_dir must be NULL or one nonempty path.", call. = FALSE)
  }
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Refusing to overwrite a nonempty output directory: ", out_dir, call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir)) stop("Cannot create output directory.", call. = FALSE)
  invisible(out_dir)
}

.tf_save_tables <- function(result, out_dir, stem) {
  if (is.null(out_dir)) return(invisible(NULL))
  .tf_output_dir(out_dir)
  saveRDS(result, file.path(out_dir, paste0(stem, ".rds")))
  for (nm in names(result)) if (is.data.frame(result[[nm]])) {
    utils::write.csv(result[[nm]], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)
  }
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  invisible(out_dir)
}
