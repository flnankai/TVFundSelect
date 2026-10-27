# Fresh computation of the six Table 4 cases: fitted TV mean + independent
# sign-symmetrized whole residual-month resampling. No stored results are inputs.
# Usage: Rscript 05_calibration.R INPUT_DIRECTORY OUTPUT_DIRECTORY [REPS] [CORES]
library(TVFundSelect)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply raw CSV input directory and NEW output directory.")
replications <- if (length(args) > 2L) as.integer(args[3L]) else 1000L
cores <- if (length(args) > 3L) as.integer(args[4L]) else 25L
data <- read_fund_panel(file.path(args[1L], "fund.csv"), file.path(args[1L], "factor.csv"),
                        file.path(args[1L], "rf.csv"))
result <- reproduce_calibration(data, replications = replications, cores = cores, out_dir = args[2L])
print(result$calibration_summary)
