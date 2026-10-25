\name{constancy_test}
\alias{constancy_test}
\alias{individual_constancy_test}
\alias{reproduce_constancy}
\title{Recompute Individual and Panel Coefficient-Constancy Diagnostics}
\description{Tests whether all alpha and factor-loading functions are constant,
using the formulas and iid whole-date residual bootstrap in Appendix B of the
October 24, 2026 manuscript. These are specification diagnostics, not tests of
positive average alpha or measurements of observed-fund FDR.}
\usage{
individual_constancy_test(y, factors, L = NULL, B = 1999L, seed = 20261020L)
constancy_test(Y, factors, L = NULL, B = 1999L, seed = 20281020L,
               cores = 25L)
reproduce_constancy(data, B = 1999L, cores = 25L, seed = 20261019L,
  panel_lengths = c(ncol(data$returns), 180L, 120L), common_length = 180L,
  min_residual_df = 30L, fund_rows = seq_len(nrow(data$returns)),
  out_dir = NULL)
}
\arguments{
\item{y}{One finite numeric fund excess-return vector.}
\item{Y}{Finite numeric N-by-T matrix of excess returns; time runs across columns.}
\item{factors}{Finite T-by-p matrix of observed factors.}
\item{L}{Spline dimension per coefficient. NULL uses 3 + floor(T\verb{^}0.2).
Quadratic splines and equally spaced interior knots are used.}
\item{B}{Number of new residual-bootstrap draws.}
\item{seed}{Nonnegative integer random seed; the caller's RNG is restored.}
\item{cores}{Maximum number of local PSOCK workers.}
\item{data}{A raw panel returned by \code{read_fund_panel()}.}
\item{panel_lengths}{Lengths of balanced panels ending at the final observation.}
\item{common_length}{Length of the common-window individual analysis.}
\item{min_residual_df}{Minimum residual degrees of freedom in the longest-history
individual sample, matching the archived sampling rule.}
\item{fund_rows}{Original input row indices to include; defaults to all funds.}
\item{out_dir}{Optional new or empty output directory; existing files are never overwritten.}
}
\details{The individual normal reference is
\eqn{Z=((T-p-1)C-(p+1)(L-1))/\sqrt{2(p+1)(L-1)}}.
The panel diagnostics use the finite-restriction and residual-correlation
normalizations printed in Appendix B. S1 uses its normal reference; S2 and the
unsquared MAX use bootstrap upper tails with inclusive exceedances and a plus-one
correction. The Cauchy statistic combines those three p-values equally.
Whole-date resampling preserves contemporaneous dependence but does not preserve
serial dependence. Fixed-design QR projection is algebraically equivalent to
refitting the constant and varying regressions; the residual correlation and all
normalizations are recomputed in every draw. No Rcpp compiler is required.
The longest uninterrupted history is selected from raw missingness only; ties
choose the earliest span. Unavailable short histories are returned explicitly.
Complete common panels retain their full calendar periods.}
\value{The two test functions return statistics, normal/bootstrap p-values, draw
statistics and configuration. \code{reproduce_constancy} returns freshly computed
individual and panel tables, sample audit, unavailable histories and bootstrap
details; optionally writes CSV/RDS output and a new p-value histogram PDF.}
\examples{
set.seed(1)
F <- matrix(rnorm(60), 60, 1)
Y <- matrix(rnorm(4 * 60), 4, 60)
test <- constancy_test(Y, F, L = 3, B = 9, cores = 1)
test$p_values
}
\seealso{\code{read_fund_panel}, \code{reproduce_paper}}
