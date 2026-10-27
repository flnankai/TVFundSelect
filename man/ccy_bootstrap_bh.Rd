\name{ccy_bootstrap_bh}
\alias{ccy_bootstrap_bh}
\title{Residual-Bootstrap Inference for Local-Linear Average Alpha}
\description{Computes reproducible bootstrap p-values and one-sided BH selections.
The public \code{ccy_bh} entry point fixes ordinary residual resampling and
the centered smoother residual degrees-of-freedom correction.}
\usage{
ccy_bootstrap_bh(Y, factors, gamma = 0.1, bandwidth = NULL,
  B = 9999L, seed = 20261027L,
  resampling = c("residual", "signed_residual"), df_adjust = TRUE,
  batch_size = 500L, keep_bootstrap = FALSE, keep_resamples = FALSE,
  cores = 25L)
}
\arguments{
\item{Y}{Finite numeric fund-by-date matrix of excess returns.}
\item{factors}{Finite numeric date-by-factor matrix, without an intercept.}
\item{gamma}{Nominal BH level, between zero and one half.}
\item{bandwidth}{Positive bandwidth or NULL for the normal-reference bandwidth.}
\item{B}{Positive integer bootstrap count, default 9999.}
\item{seed}{Nonnegative integer seed; the caller's RNG state is restored.}
\item{resampling}{Ordinary whole-date residual sampling by default. Shared-sign sampling is a separate optional diagnostic, not the paper CCY-BH method.}
\item{df_adjust}{Apply the centered smoother residual-df correction. Always true in \code{ccy_bh}.}
\item{batch_size}{Maximum bootstrap draws per in-memory batch.}
\item{keep_bootstrap}{Retain bootstrap statistics in the returned object.}
\item{keep_resamples}{Retain common date indices, signs, residual pool and null curves in memory.}
\item{cores}{PSOCK workers for a single fit, default 25; use one within an outer parallel simulation.}
}
\details{
Let A map observed returns to the fitted alpha path and S map them to fitted
returns. With a equal to the column means of A, the statistic is the unstudentized
average alpha \eqn{\widehat\delta_i=a^\top y_i}. The bootstrap null mean is
\eqn{S y_i-\widehat\delta_i\mathbf{1}}, preserving centered time-varying alpha
and all fitted loadings. Residuals are centered once and multiplied by
\eqn{\sqrt{T/\nu}}, where
\eqn{\nu=\Vert C(I-S)\Vert_F^2} and \eqn{C=I-\mathbf{1}\mathbf{1}^\top/T}.

Every bootstrap draw samples dates independently with replacement and uses
the same sampled date vector for all funds. Factors and bandwidth remain fixed.
Resampled residuals are not recentered within draws. Applying a to each
pseudo-panel is exactly equivalent to refitting the fixed-bandwidth estimator.
The primary p-value is
\eqn{(1+\sum_b 1\{\widehat\delta_i^{*(b)}\geq\widehat\delta_i\})/(B+1)}.
It retains the estimated smoothing bias of the fitted null. No normal reference
tail or studentization is used in this primary calculation.

The minimum p-value is \eqn{1/(B+1)}; compare this resolution with the first
BH threshold \eqn{\gamma/N}. All random indices are generated on the master,
so seeds, p-values and selections do not depend on worker count or batching.
No intermediate files are written. Optional returned matrices have size limits.
The ordinary iid resampling does not preserve serial dependence, and this
calibration does not by itself establish finite-sample BH FDR control.
}
\value{A list containing average alphas, plus-one p-values, exceedance counts,
BH selections, residual degrees of freedom, bandwidth, seeds, Monte Carlo
resolution and worker information. Normal-reference and centered-bootstrap
fields are diagnostic alternatives only. Optional resampling/statistic matrices
remain in memory.}
\seealso{\code{\link{ccy_bh}}, \code{\link{reproduce_simulations}},
\code{\link{reproduce_calibration}}}
\examples{
set.seed(17)
f <- matrix(rnorm(80), ncol = 1)
y <- matrix(rnorm(8 * 80), 8, 80)
# A small software example; use B=9999 for the manuscript configuration.
fit <- ccy_bh(y, f, B = 99, cores = 1)
fit$p_values
}
