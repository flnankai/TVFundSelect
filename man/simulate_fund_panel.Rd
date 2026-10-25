\name{simulate_fund_panel}
\alias{simulate_fund_panel}
\title{Generate the Manuscript's Time-Varying Fund Panels}
\description{Generates new observations from the three loading specifications,
centered time-varying alpha paths and four unit-covariance error models used in
Section 4 of the manuscript. This function never reads stored results.}
\usage{
simulate_fund_panel(N = 200L, TT = 200L, example = 3L,
  distribution = c("normal", "t3", "mixture", "icm"),
  signal_shift = 0.02, rho = 0.5, nonnull_fraction = 0.1,
  seed = NULL, burn_in = 50L, positive_range = c(0.02, 0.34),
  alpha_deviation_amplitude = 0.3)
}
\arguments{
  \item{N}{Number of funds; at least 10.}
  \item{TT}{Number of observed dates; at least 10. The estimation methods may
  require more observations to identify their local or spline regressions.}
  \item{example}{Loading specification 1, 2 or 3.}
  \item{distribution}{Normal, common-scale standardized Student t with three
  degrees of freedom, common-scale normal mixture, or independent standardized
  t components followed by the covariance transformation.}
  \item{signal_shift}{Shift added to each positive signal's uniform baseline.}
  \item{rho}{Toeplitz correlation parameter, strictly between minus one and one.}
  \item{nonnull_fraction}{Fraction of funds in each of the disjoint positive
  and negative sets. Each set has \code{floor(N * nonnull_fraction)} funds.}
  \item{seed}{An optional integer seed. If supplied, the caller's RNG state and
  kind are restored; if \code{NULL}, the current RNG stream is advanced normally.}
  \item{burn_in}{Discarded observations in the factor and stochastic-loading
  recursions, after their initialized state.}
  \item{positive_range}{Two increasing endpoints of the positive baseline.}
  \item{alpha_deviation_amplitude}{Maximum absolute uniform amplitude of the
  centered logistic deviation. Zero requests constant alpha paths.}
}
\details{
The positive average alpha is \eqn{\sqrt{200/T}(b_i+s)}, where the default
\eqn{b_i} is uniform on \eqn{[0.02,0.34]}. Negative average alphas are independently
uniform on \eqn{[-0.30,-0.20]}; all other average alphas are zero. Every alpha
path has its prescribed discrete time average exactly, up to numerical precision.
Factors follow the specified independent AR--GARCH recursions. In the normal
mixture, the scale equals one with probability 0.9 and three otherwise, divided
by \eqn{\sqrt{1.8}}. Every error model has covariance
\eqn{\Sigma_{ij}=\rho^{|i-j|}}. In particular, the ordinary multivariate t(3)
errors are divided by \eqn{\sqrt{3}}; the common-scale t model and independent
component model are distinct distributions.
}
\value{A list containing \code{Y} (fund-by-date returns), \code{factors}
(date-by-factor matrix), \code{beta}, \code{alpha}, \code{delta} (true average
alphas), \code{positive}, \code{negative}, \code{errors}, and \code{settings}.}
\seealso{\code{\link{reproduce_simulations}}, \code{\link{fit_fund_methods}}}
\examples{
panel <- simulate_fund_panel(N = 20, TT = 80, example = 1, seed = 123)
stopifnot(max(abs(rowMeans(panel$alpha) - panel$delta)) < 1e-12)
}
