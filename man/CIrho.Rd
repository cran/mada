\name{CIrho}
\alias{CIrho}
\title{
Confidence intervals for Spearman's \eqn{\rho}.
}

\concept{confidence interval}
\concept{Spearman's \eqn{\rho}}

\description{
Using Fisher's z-transformation (\code{\link{atanh}}) and the classic normal approximation confidence intervals for a vector of correlations is computed.
}
\usage{
CIrho(rho, N, level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{numeric vector, must be between -1 and 1.}
  \item{N}{integer vector, sample sizes.
}
  \item{level}{numeric, confidence level.}
}
\value{
A matrix with first column \code{rho} and two further columns with the lower and upper bound. 
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}
\examples{
CIrho(c(0.34,0.19), c(22, 48), level = 0.80)
}