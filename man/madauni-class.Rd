\name{madauni-class}


\alias{madauni-class}
\alias{summary.madauni-class}
\alias{print.madauni}
\alias{vcov.madauni}
\alias{summary.madauni}
\alias{print.summary.madauni}

\title{
Methods for the class \code{madauni}.
}
\description{
Various methods for the output of the function \code{\link{madauni}}. Also the default method \code{\link{confint}} works for this class.
}
\usage{
\method{print}{madauni}(x, digits = 3, ...)
\method{vcov}{madauni}(object, ...)
\method{summary}{madauni}(object, level = 0.95, ...)
\method{print}{summary.madauni}(x, digits = 3, ...)
}
\arguments{
  \item{x}{An object of class \code{madauni}.}
    \item{object}{An object of class \code{madauni}.}
  \item{level}{numeric, the confidence level for the confidence intervals in the summary.}
  \item{digits}{integer indicating the number of decimal places to round to.}
  \item{\dots}{arguments to be passed to methods}
}

\value{\code{summary.madauni} returns a list of class \code{summary.madauni} which is printed with \code{print.summary.madauni}.}
\author{Philipp Doebler <philipp.doebler@googlemail.com>}
\seealso{
\code{\link{madauni}}
}
\keyword{classes}% __ONLY ONE__ keyword per line
