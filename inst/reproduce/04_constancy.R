library(TVFundSelect)
args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript 04_constancy.R DATA_DIR [OUT_DIR] [CORES] [B]")
d <- args[1L]
out <- if (length(args) >= 2L) args[2L] else "recomputed/constancy"
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 25L
B <- if (length(args) >= 4L) as.integer(args[4L]) else 1999L
data <- read_fund_panel(file.path(d, "fund.csv"), file.path(d, "factor.csv"), file.path(d, "rf.csv"))
reproduce_constancy(data, B = B, cores = cores, out_dir = out)
