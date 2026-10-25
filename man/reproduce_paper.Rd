\name{reproduce_paper}
\alias{reproduce_paper}
\alias{paper_settings}
\alias{paper_targets}
\alias{example_fund_data}
\title{Recalculate Article Results from Generators and Raw Data}
\description{An end-to-end workflow that generates new panels, refits every
model and produces new tables and plots. It does not load archived result tables,
models, selections or original-author checkpoints.}
\usage{
paper_settings()
paper_targets()
example_fund_data(N = 30L, TT = 240L, seed = 20261025L)
reproduce_paper(data = NULL, output_dir = "paper_recomputed",
  components = c("simulations", "bandwidth", "constancy", "calibration", "empirical"),
  cores = 25L, simulation_reps = 1000L, calibration_reps = 1000L,
  bootstrap_reps = 1999L, seed = 420261010L, quick = FALSE)
}
\arguments{
\item{data}{Raw empirical input from \code{read_fund_panel}. Required for a full
empirical run. Quick mode may create clearly labeled synthetic demonstration data.}
\item{output_dir}{New or empty output directory. No existing results are overwritten.}
\item{components}{Subsets of the five named fresh-computation modules.}
\item{cores}{Maximum local PSOCK workers, default 25.}
\item{simulation_reps}{New replications per simulation setting, default 1000.}
\item{calibration_reps}{New replications per empirical calibration case, default 1000.}
\item{bootstrap_reps}{New constancy bootstrap draws, default 1999.}
\item{seed}{Main simulation seed or synthetic demonstration seed. Other modules
retain their documented independent article-aligned defaults.}
\item{quick}{Use small fresh simulation/bootstrap counts and synthetic empirical
data if none is supplied. This demonstrates the software, not paper evidence.}
\item{N}{Number of synthetic demonstration funds.}
\item{TT}{Number of demonstration months.}
}
\details{\code{paper_targets} maps each article table and figure to its computing
function. The full Monte Carlo defaults are deliberately expensive: 216000 main
panels with five method fits each, a separate bandwidth experiment, and empirical
calibration and bootstrap work. Numerical failures are reported, not dropped.
Simulation seeds produce the same statistical design but need not reproduce the
legacy Monte Carlo table digit for digit. Use the supplied R/version/configuration
records to compare uncertainty. Raw CRSP-derived inputs are not bundled.
\code{example_fund_data} creates raw-format data for software examples only;
they cannot reproduce the empirical tables.}
\value{Settings and targets return a list and data frame. The demonstration
generator returns a synthetic raw-data list. \code{reproduce_paper} invisibly
returns all freshly computed module outputs and a manifest, and writes component
results and Table 1--9 CSVs corresponding to the requested targets.}
\examples{
paper_targets()
x <- example_fund_data(N = 12, TT = 120)
dim(x$returns)
\dontrun{
# Recalculates a small complete workflow. No archived results are read.
reproduce_paper(output_dir = "fresh_demo", quick = TRUE, cores = 2)
# Full article calculations after providing authorized raw empirical data:
reproduce_paper(data = raw_panel, output_dir = "fresh_full", cores = 25)
}
}
