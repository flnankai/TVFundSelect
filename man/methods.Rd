\name{fund_selection_methods}
\alias{fund_selection_methods}
\alias{tv_ss_bh}
\alias{tv_ls_bh}
\alias{cc_ss_bh}
\alias{cc_ls_bh}
\alias{ccy_bh}
\alias{fit_fund_methods}
\title{Select Funds with Positive Average Alpha}
\description{
Five signed, one-sided fund-selection procedures from the accompanying paper.
Time-varying spline, constant-coefficient, spatial-median, least-squares and
local-linear implementations are separately callable on user-supplied panels.
}
\usage{
tv_ss_bh(Y, factors, gamma = 0.1, knot_candidates = 0:3,
  degree = 2L, bic_penalty = c("paper", "standard"),
  tolerance = 1e-6, max_iter = 10000L,
  calibration = c("residual_df", "none"))

tv_ls_bh(Y, factors, gamma = 0.1, knot_candidates = 0:3,
  degree = 2L, bic_penalty = c("paper", "standard"))

cc_ss_bh(Y, factors, gamma = 0.1, tolerance = 1e-6,
  max_iter = 10000L)

cc_ls_bh(Y, factors, gamma = 0.1,
  variance_estimator = c("residual_df", "paper"))

ccy_bh(Y, factors, gamma = 0.1, bandwidth = NULL,
  keep_curves = FALSE, B = 9999L, seed = 20261027L, cores = 25L,
  batch_size = 500L, keep_bootstrap = FALSE, keep_resamples = FALSE)

fit_fund_methods(Y, factors, gamma = 0.1,
  knot_candidates = 0:3, degree = 2L, tolerance = 1e-6,
  max_iter = 10000L, ccy_bandwidth = NULL,
  methods = c("TV-SS-BH", "TV-LS-BH", "CC-SS-BH", "CC-LS-BH", "CCY-BH"),
  ccy_B = 9999L, ccy_seed = 20261027L, ccy_cores = 25L)
}
\arguments{
\item{Y}{Finite numeric \eqn{N\times T} matrix of excess returns. Funds are
rows and dates are columns. Missing observations must be handled before fitting.}
\item{factors}{Finite numeric \eqn{T\times p} matrix of observed factors,
without an intercept. A zero-column matrix is supported when identified.}
\item{gamma}{Nominal BH level, strictly between zero and one half.}
\item{knot_candidates}{Distinct nonnegative integer counts of equally spaced
interior knots. A single count fixes the dimension instead of searching.}
\item{degree}{Positive integer polynomial degree of the B-splines; the paper uses two.}
\item{bic_penalty}{The paper's restricted-residual BIC by default. The alternative
\code{"standard"} is a diagnostic and is not the paper's primary specification.}
\item{tolerance}{Joint spatial location/scale equation tolerance.}
\item{max_iter}{Positive integer iteration limit; failure to converge is an error.}
\item{calibration}{The default applies the paper's residual-df calibration once.
\code{"none"} supplies an explicitly uncalibrated diagnostic.}
\item{variance_estimator}{The default uses the unrestricted constant regression's
residual degrees of freedom. \code{"paper"} uses RSS/T and is not the paper's
primary CC-LS-BH comparator.}
\item{bandwidth, ccy_bandwidth}{Positive local-linear bandwidth on normalized time.
\code{NULL} uses \eqn{1.06\,\mathrm{sd}(t/T)T^{-1/5}} with sample SD.}
\item{keep_curves}{Retain fitted alpha curves and residuals in memory.}
\item{B, ccy_B}{Number of bootstrap draws for CCY, default 9999.}
\item{seed, ccy_seed}{Reproducible bootstrap seed; caller RNG state is restored.}
\item{cores, ccy_cores}{PSOCK workers for one CCY fit, default 25. Outer parallel workflows use one internally.}
\item{batch_size}{Bootstrap draws per in-memory batch.}
\item{keep_bootstrap, keep_resamples}{Retain bootstrap statistics or date indices in the returned object, never on disk.}
\item{methods}{A nonempty, nonduplicated subset of the five canonical labels.}
}
\details{
The target is \eqn{\delta_i=T^{-1}\sum_t\alpha_i(t/T)} and the null is
\eqn{\delta_i\leq0}. This is not a test that the alpha path is identically zero.
Only the alpha-deviation spline block is centered. Both time-varying spline
procedures choose the same design using restricted residuals without an intercept.

TV-SS-BH fits the joint spatial-median and positive diagonal scale equations to
restricted residuals. Its signed score is
\eqn{\sqrt{(T-q)\widehat\varsigma}\,\widehat\theta_i/\widehat d_i},
where \eqn{q} is nuisance rank, excluding the intercept. The returned
\code{varsigma} includes the df multiplier; \code{varsigma_uncorrected} does not.
The scale is normalized by a geometric-mean convention, which cancels in scores.
It requires at least two funds. CC-SS-BH uses the original constant-coefficient
spatial-sign normalization, without an additional time-varying df multiplier.

TV-LS-BH estimates the full regression on \code{[1,Z]} and divides unrestricted
RSS by \eqn{T-q-1}. CC-LS-BH analogously uses \eqn{T-p-1}.
The four non-CCY methods use upper standard-normal tails. CCY-BH instead uses
the ordinary residual-bootstrap tail probability described in
\code{\link{ccy_bootstrap_bh}}, with residual degrees-of-freedom correction.
Its public \code{statistics} are unstudentized average alphas; normal statistics
and standard errors are diagnostic fields only and do not enter its BH rule.
The bootstrap extension is specified by this package and the revised manuscript;
it is not attributed as an unchanged testing procedure to Cheng, Cai and Yan.
Bootstrap calibration alone does not prove BH FDR control under arbitrary dependence.
}
\value{
Each individual method returns a list including \code{method}, integer row indices
\code{selected}, discovery count \code{k}, \code{bh_threshold}, signed
\code{statistics}, \code{p_values}, \code{gamma}, and method-specific diagnostics.
Least-squares and local-linear methods return \code{delta_hat} and
\code{standard_errors}; spatial methods return \code{theta}, \code{d}, radial
moments, equation residuals and convergence information. Spline fits also return
selected knot counts, design rank and the BIC table. Selection indices always
refer to the rows of the supplied \code{Y}.
\code{fit_fund_methods} returns a named list of individual results in requested order.
}
\references{
Wang, H., Zhao, P., Feng, L. and Wang, Z. (2025).
Robust mutual fund selection with false discovery rate control.
\emph{Journal of Econometrics}, 252, 106121.

Cheng, T., Cai, B. and Yan, C. (2018).
Time-varying skills (versus luck) in U.S. active mutual funds and hedge funds.
\emph{Journal of Empirical Finance}, 49, 81--106.

Ang, A. and Kristensen, D. (2012). Testing conditional factor models.
\emph{Journal of Financial Economics}, 106, 132--156.
}
\seealso{\code{\link{make_tv_design}}, \code{\link{bh_select}}}
\examples{
set.seed(91)
f <- matrix(rnorm(120), ncol = 1)
y <- matrix(rnorm(12 * 120), 12, 120) + 0.4 * rep(f[, 1], each = 12)
fit <- tv_ls_bh(y, f, knot_candidates = 1)
fit$selected
both <- fit_fund_methods(y, f, knot_candidates = 1,
  methods = c("TV-LS-BH", "CC-LS-BH"))
vapply(both, function(x) x$k, integer(1))
}
