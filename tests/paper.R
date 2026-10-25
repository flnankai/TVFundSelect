library(TVFundSelect)
set.seed(123)
before <- .Random.seed
before_locale <- Sys.getlocale()
d <- example_fund_data(N = 12L, TT = 120L)
stopifnot(identical(.Random.seed,before),isTRUE(d$synthetic),
          identical(dim(d$returns),c(12L,120L)),nrow(paper_targets()) == 13L,
          paper_settings()$replications == 1000L)
fails_paper <- function(x) stopifnot(inherits(tryCatch(force(x),error=identity),"error"))
fails_paper(reproduce_paper(output_dir=tempfile(),components="empirical"))
fails_paper(reproduce_paper(components="not_a_component"))

# A small end-to-end newly computed workflow, not archived plotting.
out <- tempfile("fresh_paper_")
r <- reproduce_paper(output_dir=out,quick=TRUE,cores=1L)
stopifnot(isTRUE(r$manifest$fresh_computation),isTRUE(r$manifest$synthetic_data),
          identical(.Random.seed,before),length(r$simulations$tasks$task_seed)==8L,
          identical(Sys.getlocale(),before_locale),
          nrow(r$calibration$replication_metrics)>0L,
          nrow(r$empirical$fund_scores)>0L,
          all(file.exists(file.path(out,paste0("Table_",1:9,".csv")))))
fails_paper(reproduce_paper(output_dir=out,quick=TRUE,cores=1L))
stopifnot(file.exists(file.path(out,"run_manifest.rds")))
cat("Fresh article orchestration tests passed.\n")
