#!/usr/bin/env Rscript
# Refit CCY-BH on NEW independent panels for every bandwidth multiplier.
# Rscript 02_bandwidth.R --quick
# Rscript 02_bandwidth.R --reps=1000 --cores=25 --out=paper_bandwidth
library(TVFundSelect)
args <- commandArgs(trailingOnly = TRUE)
opt <- list(reps = 1000L, cores = 25L, out = "paper_bandwidth")
quick <- "--quick" %in% args
for (arg in setdiff(args, "--quick")) {
  pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
  if (length(pair) != 2L || !pair[1L] %in% names(opt)) stop("Invalid option: ", arg)
  opt[[pair[1L]]] <- if (pair[1L] == "out") pair[2L] else as.integer(pair[2L])
}
if (quick) {
  result <- reproduce_bandwidth(reps = 2L, cores = min(opt$cores, 2L),
    examples = 1L, times = 80L, distributions = "normal",
    signals = c(0.04, 0.10), N = 40L, output_dir = opt$out, chunk_size = 1L)
} else {
  result <- reproduce_bandwidth(reps = opt$reps, cores = opt$cores, output_dir = opt$out)
}
print(result$summary)
