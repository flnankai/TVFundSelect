\name{reproduce_calibration}
\alias{reproduce_calibration}
\title{Fresh data-calibrated false-discovery experiments}
\description{Fits a generating model and residual pool from raw user-supplied fund
data, generates every Monte Carlo panel anew, and refits every selection method.}
\usage{
reproduce_calibration(data, replications = 1000L, cores = 25L,
  lookbacks = c(120L, 180L), train_end = 366L,
  shifts = c(0, .002, .004), generators = "TV",
  methods = five_method_names, signal_fraction = .1,
  gamma = .1, knot_candidates = 0:3, degree = 2L,
  block_length = 1L, seed = 31019000L, truth_seed = 20261019L,
  tolerance = 1e-6, max_iter = 10000L, out_dir = NULL,
  ccy_B = 9999L, ccy_seed = 202610270L)
}
\arguments{
\item{data}{Panel from \code{read_fund_panel}.}
\item{replications}{Fresh panels per case; 1000 for the article.}
\item{cores}{Requested PSOCK workers; default 25.}
\item{lookbacks}{Calibration sample lengths; 120 and 180 in the article.}
\item{train_end}{Last calibration month index. Index 366 is July 2017 for the article archive.}
\item{shifts}{Absolute positive/negative monthly average-alpha magnitudes. Zero is the all-null case.}
\item{generators}{\code{"TV"} for Table 4, optionally \code{c("CC", "TV")} for the twelve-case extended experiment.}
\item{methods}{Selection labels supported by \code{fit_fund_methods}.}
\item{signal_fraction}{Fraction assigned each positive and negative sign, rounded down. At N=504 and 0.1 this is 50 per side.}
\item{gamma}{Nominal BH level.}
\item{knot_candidates}{Restricted-residual BIC candidate interior-knot counts.}
\item{degree}{Spline degree.}
\item{block_length}{Length of circular residual blocks. Default 1 matches the independent-month article experiment.}
\item{seed}{Simulation base seed; each task uses \code{seed + T*10000 + replication}.}
\item{truth_seed}{Fixed sign-assignment base seed; each lookback uses \code{truth_seed + T}.}
\item{tolerance}{Spatial equation convergence tolerance.}
\item{max_iter}{Maximum spatial iterations.}
\item{ccy_B, ccy_seed}{CCY bootstrap draws and base seed. Each fit uses ccy_seed + T*10000 + replication; the inner bootstrap is unsigned ordinary residual resampling.}
\item{out_dir}{Optional empty output directory. Existing files are not overwritten.}
}
\details{The generator is fitted from newly loaded raw data, with complete-history
eligibility and excess-return response. Only the fitted average-alpha intercept
is removed. Centered time-varying alpha and risk exposures remain. Generating
residuals are centered once and multiplied by \code{sqrt(T/(T-rank(X)))}.
Each generated date samples a whole residual vector and a shared Rademacher sign.
The generated panel is not recentered. All methods share that panel, and shifted
cases share draws. The default TV generator and two sample lengths produce six
cases; extended CC and TV generators produce twelve. All-null TDP is undefined;
its FDR equals the probability of any selection. These conditional experiments
do not estimate the unknown FDR of observed fund choices and do not propagate
generating-model estimation uncertainty. Supplied seeds determine draws before
method fitting, independently of worker scheduling. The caller's RNG is restored.
No failed method fit is silently discarded; failures stop the computation.
Only the final calibration summary is saved; replication details remain in memory.}
\value{A list of fresh replication metrics, Monte Carlo summaries and MCSEs,
newly fitted designs, truth assignments, initial draw audits, seeds and settings.}
\seealso{\code{read_fund_panel}, \code{reproduce_empirical}}
