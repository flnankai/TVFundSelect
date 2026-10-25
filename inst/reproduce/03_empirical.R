# Fresh computation of manuscript Tables 6--9 and Figure 8 from raw CSV inputs.
# Usage: Rscript 03_empirical.R INPUT_DIRECTORY OUTPUT_DIRECTORY [CORES]
library(TVFundSelect)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply a directory containing fund.csv/factor.csv/rf.csv, and a NEW output directory.")
cores <- if (length(args) > 2L) as.integer(args[3L]) else 25L
data <- read_fund_panel(file.path(args[1L], "fund.csv"), file.path(args[1L], "factor.csv"),
                        file.path(args[1L], "rf.csv"))
result <- reproduce_empirical(data, cores = cores, out_dir = args[2L])
print(result$performance_summary)
print(result$strategy_alpha)
