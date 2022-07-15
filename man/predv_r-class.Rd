\name{predv_r-class}


\alias{predv_r-class}
\alias{summary.predv_r-class}
\alias{print.predv_r}
\alias{summary.predv_r}

\title{
Methods for the class \code{predv_r}.
}
\description{
Various methods for the output of the function \code{\link{predv_r}}.
}
\usage{
\method{print}{predv_r}(x, ylim_npv=c(0,1),ylim_ppv=c(0,1), ...)
\method{summary}{predv_r}(object, ylim_npv=c(0,1),ylim_ppv=c(0,1), ...)
}
\arguments{
  \item{x}{An object of class \code{predv_r}.}
    \item{object}{An object of class \code{predv_r}.}
  \item{ylim_npv}{limits of the y-axis for the plot on projected negative predictive values. Default is c(0,1).}
  \item{ylim_ppv}{limits of the y-axis for the plot on projected positive predictive values. Default is c(0,1).}
  \item{\dots}{arguments to be passed to methods}
}

\value{\code{summary.predv_r} returns a list of class \code{summary.predv_r}.}
\author{Bernardo Sousa-Pinto <bernardo@med.up.pt>}
\seealso{
\code{\link{predv_r}}
}
\keyword{classes}