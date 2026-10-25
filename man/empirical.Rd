\name{read_fund_panel}
\alias{read_fund_panel}
\alias{reproduce_empirical}
\alias{portfolio_alpha}
\alias{plot_empirical}
\title{Fresh rolling fund selection and portfolio evaluation}
\description{Loads user-supplied raw return CSVs and recomputes rolling selections,
buy-and-hold portfolio returns, paired performance comparisons and spline-dimension
sensitivity. Historical selections and archived result summaries are not read.}
\usage{
read_fund_panel(fund_file, factor_file, rf_file, start_date = "1987-02-01",
  dates = NULL, zero_is_missing = TRUE, rf_in_percent = TRUE, header = FALSE)
reproduce_empirical(data, lookbacks = c(120L, 180L), holding = 6L,
  methods = c("TV-SS-BH", "TV-LS-BH"), gamma = .1,
  knot_candidates = 0:3, degree = 2L, tolerance = 1e-6,
  max_iter = 10000L, cores = 25L, benchmarks = c("FF4", "TV-FF4"),
  hac_lags = 6L, sensitivity = TRUE, sensitivity_knots = 0:3,
  out_dir = NULL)
portfolio_alpha(returns, factors, rf = rep(0, length(returns)),
  benchmark = c("FF4", "TV-FF4"), hac_lag = 6L,
  n_internal = 3L, degree = 2L, spread = FALSE)
plot_empirical(x, file)
}
\arguments{
\item{fund_file}{CSV with funds in rows and dates in columns; net returns in decimal units.}
\item{factor_file}{CSV with dates in rows and factors in columns; decimal units, without an intercept.}
\item{rf_file}{One-row or one-column CSV of risk-free returns.}
\item{start_date}{First month, used if \code{dates} is omitted.}
\item{dates}{Optional increasing vector convertible to \code{Date}.}
\item{zero_is_missing}{Whether raw zero-coded fund returns denote missing data. The mask is fixed before RF subtraction.}
\item{rf_in_percent}{Whether RF values need division by 100; the article archive uses percentages.}
\item{header}{Whether the three CSV inputs contain column headers.}
\item{data}{A panel returned by \code{read_fund_panel}.}
\item{lookbacks}{Training lengths in months. The article uses 120 and 180.}
\item{holding}{Months between formation dates; the last incomplete holding period is retained.}
\item{methods}{Selection labels supported by \code{fit_fund_methods}. The article's investment comparison uses the two TV methods.}
\item{gamma}{Nominal one-sided BH FDR level.}
\item{knot_candidates}{Candidate interior-knot counts for the restricted-residual BIC rule.}
\item{degree}{Spline polynomial degree.}
\item{tolerance}{Spatial equation convergence tolerance.}
\item{max_iter}{Maximum spatial iterations.}
\item{cores}{Requested PSOCK worker count; at most the number of tasks is started. Default 25.}
\item{benchmarks, benchmark}{Constant-factor or centered time-varying factor evaluation regression. The historical label FF4 denotes the supplied factor matrix; supply the article's four columns for exact specifications.}
\item{hac_lags, hac_lag}{Bartlett HAC lag(s), with no prewhitening and regression-df adjustment.}
\item{sensitivity}{Recalculate fixed-knot rolling fits for the TV methods.}
\item{sensitivity_knots}{Fixed interior-knot counts for Table 9 style sensitivity.}
\item{out_dir}{Optional empty output directory for newly calculated CSVs, RDS and plot. Nonempty directories are never overwritten.}
\item{returns}{Finite monthly portfolio net returns, or paired return differences if \code{spread=TRUE}.}
\item{factors}{Finite date-by-factor matrix for economic evaluation.}
\item{rf}{Monthly risk-free returns in decimal units.}
\item{n_internal}{Interior-knot count for the common TV evaluation regression.}
\item{spread}{Do not subtract RF from a paired return difference.}
\item{x}{Object returned by \code{reproduce_empirical}.}
\item{file}{Output PDF path for the newly calculated wealth/selection plot.}
}
\details{Eligibility uses training observations only. Selected funds receive equal
weights at formation and are then bought and held, so weights drift. Empty
selections earn RF. A selected fund with an unavailable holding return is moved
to RF from that month through the end of the holding block; this is a valuation
assumption, not an observed liquidation return. No transaction costs are added.
The response is always fund net return minus RF. Annual alpha and SE equal
twelve times monthly values. Two-sided economic p-values differ from the
one-sided selection p-values. Paired TV-SS minus TV-LS inference uses the matched
monthly difference directly; the single contrast within each lookback is unadjusted.
Sensitivity overlap is the mean Jaccard ratio conditional on a nonempty union.
Raw CRSP-derived fund data are deliberately not redistributed with the package.}
\value{\code{read_fund_panel} returns a validated panel with missingness, dates and
input hashes. \code{reproduce_empirical} returns freshly computed monthly returns,
window selections, fund scores, performance, benchmark alpha, paired alpha,
sensitivity and configuration tables. \code{portfolio_alpha} returns one evaluation
row. \code{plot_empirical} invisibly returns the output path.}
\examples{
set.seed(1)
f <- matrix(rnorm(80 * 4, sd = .02), 80, 4)
r <- .001 + drop(f \%*\% c(.8, .1, -.1, .2)) + rnorm(80, sd = .01)
portfolio_alpha(r, f, hac_lag = 6)
}
