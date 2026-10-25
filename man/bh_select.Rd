\name{bh_select}
\alias{bh_select}
\alias{discovery_metrics}
\title{BH Selection and Simulation Discovery Metrics}
\description{Apply the Benjamini--Hochberg step-up rule to valid p-values or
evaluate simulated selections against known average-alpha truths.}
\usage{
bh_select(p_values, gamma = 0.1)
discovery_metrics(selected, delta)
}
\arguments{
\item{p_values}{Finite numeric vector in [0,1]. An empty vector gives no selections.}
\item{gamma}{Nominal level strictly between zero and one half.}
\item{selected}{Integer row indices of selected funds.}
\item{delta}{Nonempty finite numeric vector of true average alphas.}
}
\details{
False selections include both zero and negative average alphas. TDP is the
fraction recovered among all positive average-alpha funds, not a separately
designated asymptotically detectable subset. With no positive truths, TDP is
returned as zero. Real empirical data do not provide known \code{delta};
these metrics must not be interpreted as observed empirical skill or FDR.
}
\value{
\code{bh_select} returns \code{selected}, \code{k}, and \code{threshold}.
\code{discovery_metrics} returns \code{FDP}, \code{TDP}, the numbers selected,
true discoveries, false discoveries, alternatives and nulls.
}
\examples{
b <- bh_select(c(0.001, 0.02, 0.9), gamma = 0.1)
discovery_metrics(b$selected, c(0.2, 0, -0.1))
}
