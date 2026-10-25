library(TVFundSelect)
args <- commandArgs(TRUE)
quick <- "--quick" %in% args
args <- setdiff(args, "--quick")
if (!length(args)) stop("Usage: Rscript 00_all.R OUT_DIR [DATA_DIR] [CORES] [--quick]")
data <- if (length(args) >= 2L) {
  d <- args[2L]
  read_fund_panel(file.path(d,"fund.csv"),file.path(d,"factor.csv"),file.path(d,"rf.csv"))
} else NULL
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 25L
reproduce_paper(data, output_dir = args[1L], cores = cores, quick = quick)
