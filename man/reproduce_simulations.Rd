\name{reproduce_simulations}
\alias{reproduce_simulations}
\alias{reproduce_bandwidth}
\alias{plot_simulation_results}
\alias{print.fund_simulation}
\title{Recompute the Simulation Studies from Newly Generated Panels}
\description{Fresh Monte Carlo calculations for the five-method comparison in
Section 4. A separate optional bandwidth diagnostic is retained for compatibility
but is not included in the revised article or default paper workflow. Every
task generates observations and refits the requested methods. Archived article
summaries are not inputs to these functions.}
\usage{
reproduce_simulations(reps = 1000L, cores = 25L, examples = 1:3,
  times = c(200L, 400L),
  distributions = c("normal", "t3", "mixture", "icm"),
  signals = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.175, 0.23, 0.30),
  N = 200L, gamma = 0.1, seed = 420261010L, rho = 0.5,
  nonnull_fraction = 0.1, knot_candidates = 0:3, degree = 2L,
  tolerance = 1e-6, max_iter = 10000L,
  methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"),
  output_dir = NULL, plots = TRUE, resume = FALSE, chunk_size = 20L,
  ccy_B = 9999L, ccy_seed = 20261027L)

reproduce_bandwidth(reps = 1000L, cores = 25L, examples = 1:3,
  times = c(200L, 400L),
  distributions = c("normal", "t3", "mixture", "icm"),
  signals = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.175, 0.23, 0.30),
  multipliers = c(1, 2, 4, 8, 16), N = 200L, gamma = 0.1,
  seed = 1120261010L, rho = 0.5, nonnull_fraction = 0.1,
  output_dir = NULL, plots = TRUE, resume = FALSE, chunk_size = 20L,
  ccy_B = 9999L, ccy_seed = 20261027L)

plot_simulation_results(x, file = NULL, overwrite = FALSE)

\method{print}{fund_simulation}(x, ...)
}
\arguments{
  \item{reps}{Independent generated panels per setting. Defaults to 1000.}
  \item{cores}{Number of local PSOCK workers; default 25. Use one for serial
  computation. BLAS/OpenMP environment limits are one per worker.}
  \item{examples}{Distinct loading examples drawn from 1, 2 and 3.}
  \item{times}{Distinct sample lengths. The article uses 200 and 400.}
  \item{distributions}{Distinct error-model labels from
  \code{c("normal", "t3", "mixture", "icm")}.}
  \item{signals}{Distinct positive-average-alpha shifts.}
  \item{N}{Number of funds; the article uses 200.}
  \item{gamma}{Nominal BH level, strictly between zero and one half.}
  \item{seed}{Master seed used to assign a reproducible integer seed to each
  panel before dispatch; task seeds are returned. The two defaults separate
  the main and bandwidth experiments.}
  \item{rho}{Toeplitz error-correlation parameter.}
  \item{nonnull_fraction}{Fraction in each positive and negative signal set.}
  \item{knot_candidates}{Distinct nonnegative integer internal-knot counts.}
  \item{degree}{Spline polynomial degree, two in the article.}
  \item{tolerance, max_iter}{Spatial-median/scale convergence controls.}
  \item{methods}{Distinct paper method labels to refit on each generated panel.}
  \item{multipliers}{Positive distinct multiples of the CCY reference bandwidth.}
  \item{output_dir}{Optional new output directory. With \code{NULL}, all
  calculations still take place but no files are written.}
  \item{plots}{Whether to write freshly calculated curves when
  \code{output_dir} is supplied.}
  \item{resume}{Must be false. Intermediate checkpoints are not written or reused.}
  \item{chunk_size}{Panels per in-memory dispatch chunk.}
  \item{ccy_B}{CCY bootstrap draws per panel, default 9999.}
  \item{ccy_seed}{Base CCY seed; each task uses this plus 2029 times its task ID.}
  \item{x}{A result returned by either simulation runner.}
  \item{file}{Optional PDF path; \code{NULL} draws on the current device.}
  \item{overwrite}{Whether to replace an existing figure, default false.}
  \item{...}{Unused printing arguments.}
}
\details{
With the defaults there are 216 settings and 216000 independently generated
panels. The main experiment performs five fits per panel. The bandwidth
experiment instead fits CCY-BH at each bandwidth multiplier on the same panel;
its panels are independent of the main experiment because of its separate seed.
For each task, the seed is \code{seed + 1009 * task_id}, where task IDs enumerate
\code{expand.grid(distribution, signal_shift, replication, T, example)} with
the first coordinate varying fastest. This documented schedule produces the
same statistical experiment as the manuscript but is not the legacy batch seed
schedule. Newly computed estimates should be compared within Monte Carlo
uncertainty, rather than expected to match every printed digit.

The actual task seeds and full configuration are returned; scheduling, worker
count and chunk size do not change a task's random observations. Changing the
specified grid changes task numbering and hence the associated seeds. The
caller RNG state and thread environment settings are restored. A fitting
failure stops the computation with the task seed and setting; failures are
never silently deleted when averaging results.
On Windows, file-output operations temporarily use a UTF-8 character locale
and restore the caller's exact \code{LC_CTYPE} on return or error; numeric and
collation locales are not changed. This permits non-ASCII output directories
and final summaries even when the calling process starts in the C locale.

For selected set S, FDP counts selections with nonpositive true average alpha,
divided by the number selected or one when empty. TDP divides true positive
discoveries by all positive-average-alpha funds, not the theoretical detectable
subset. Summaries include Monte Carlo means and standard errors. No figure is
drawn from frozen article CSV files.

Full runs can be computationally expensive. The examples below intentionally
request a reduced experiment; use the defaults for the paper design. The
provided output directory stores only \file{simulation_summary.txt} and an
optional final \file{fdp_tdp.pdf}. Raw replication metrics and configurations
remain available in the returned R object, not in files. Results are based on independent
innovations across dates; this does not validate arbitrary serial dependence.
}
\value{The runners return an object of class \code{fund_simulation} containing
\code{raw} replication metrics, \code{summary} means and Monte Carlo standard
errors, \code{tasks} and their seeds, \code{config}, \code{manifest},
\code{cores}, and \code{output_dir}. The plot function invisibly returns its
output path, or \code{NULL} for the current device.}
\seealso{\code{\link{simulate_fund_panel}}, \code{\link{fit_fund_methods}}}
\examples{
small <- reproduce_simulations(reps = 2, cores = 1, examples = 1,
  times = 80, distributions = "normal", signals = 0.10, N = 20,
  methods = c("TV-LS-BH", "CC-LS-BH"), plots = FALSE)
small$summary
}
