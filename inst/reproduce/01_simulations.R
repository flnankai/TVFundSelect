#!/usr/bin/env Rscript
# Full default = freshly generated 216000 panels, 25 PSOCK workers.
# Rscript 01_simulations.R --quick
# Rscript 01_simulations.R --reps=1000 --cores=25 --out=paper_simulations
library(TVFundSelect)
args <- commandArgs(trailingOnly = TRUE)
opt <- list(reps = 1000L, cores = 25L, out = "paper_simulations")
quick <- "--quick" %in% args
for (arg in setdiff(args, "--quick")) {
  pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
  if (length(pair) != 2L || !pair[1L] %in% names(opt)) stop("Invalid option: ", arg)
  opt[[pair[1L]]] <- if (pair[1L] == "out") pair[2L] else as.integer(pair[2L])
}
if (quick) {
  result <- reproduce_simulations(reps = 2L, cores = min(opt$cores, 2L),
    examples = 1L, times = 80L, distributions = c("normal", "t3"),
    signals = c(0.04, 0.10), N = 40L, output_dir = opt$out, chunk_size = 1L)
} else {
  result <- reproduce_simulations(reps = opt$reps, cores = opt$cores, output_dir = opt$out)
}
print(result$summary)
