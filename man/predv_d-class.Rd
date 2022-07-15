\name{predv_d-class}


\alias{predv_d-class}
\alias{summary.predv_d-class}
\alias{print.predv_d}
\alias{summary.predv_d}

\title{
Methods for the class \code{predv_d}.
}
\description{
Various methods for the output of the function \code{\link{predv_d}}.
}
\usage{
\method{print}{predv_d}(x, xlim_npv=c(0,1),xlim_ppv=c(0,1), ...)
\method{summary}{predv_d}(object, xlim_npv=c(0,1),xlim_ppv=c(0,1), ...)
}
\arguments{
  \item{x}{An object of class \code{predv_d}.}
    \item{object}{An object of class \code{predv_d}.}
  \item{xlim_npv}{limits of the x-axis for the plot on projected negative predictive values. Default is c(0,1).}
  \item{xlim_ppv}{limits of the x-axis for the plot on projected positive predictive values. Default is c(0,1).}
  \item{\dots}{arguments to be passed to methods}
}

\value{\code{summary.predv_d} returns a list of class \code{summary.predv_d}.}
\author{Bernardo Sousa-Pinto <bernardo@med.up.pt>}
\seealso{
\code{\link{predv_d}}
}
\keyword{classes}