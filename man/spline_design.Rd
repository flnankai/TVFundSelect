\name{make_tv_design}
\alias{make_tv_design}
\alias{select_tv_design_bic}
\title{Construct or Select the Centered Time-Varying Spline Design}
\description{Construct the rank-reduced centered alpha-deviation block and
factor-by-spline blocks, or select the number of interior knots using the
paper's restricted-residual BIC.}
\usage{
make_tv_design(factors, n_internal = 2L, degree = 2L)
select_tv_design_bic(Y, factors, knot_candidates = 0:3,
  degree = 2L, penalty = c("paper", "standard"))
}
\arguments{
\item{Y}{Finite numeric N by T excess-return matrix.}
\item{factors}{Finite numeric T by p factor matrix, without an intercept.}
\item{n_internal}{Nonnegative integer number of equally spaced interior knots.}
\item{degree}{Positive integer spline polynomial degree.}
\item{knot_candidates}{Distinct nonnegative integer knot counts.}
\item{penalty}{Use \code{"paper"} for article reproduction;
\code{"standard"} is an additional diagnostic convention.}
}
\details{
Time is \code{seq_len(T)/T}. For a basis with L columns, a Helmert contrast
removes the redundant centered-alpha direction. Only this block is centered;
factor spline blocks are not. For full rank, nuisance rank is \eqn{(p+1)L-1}.
The augmented design including the intercept must retain full rank. Inadmissible
candidates are marked in the returned BIC table, not repaired by ridge fitting.

The paper criterion is
\eqn{\log(\lVert YM_K\rVert_F^2/(NT))+(p+1)L\log(NT)/(NT)}.
This uses the restricted regression without the average-alpha intercept.
}
\value{
\code{make_tv_design} returns \code{Z}, the original basis \code{B}, centered
basis \code{B_tilde}, reduced alpha block \code{C}, orthonormal contrast
\code{A_L}, basis dimension \code{L}, \code{rank}, QR decomposition and spline settings.
\code{select_tv_design_bic} returns the selected \code{design}, restricted
\code{residual} information, BIC, MSE, candidate table and penalty convention.
}
\seealso{\code{\link{tv_ss_bh}}, \code{\link{tv_ls_bh}}}
\examples{
set.seed(12)
f <- matrix(rnorm(80), ncol = 1)
d <- make_tv_design(f, n_internal = 1)
dim(d$Z)
}
