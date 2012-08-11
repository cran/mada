\name{talpha}
\alias{talpha}
\title{
The \eqn{t_\alpha} transformation as a link function for binary GLMs.
}

\concept{transformation}
\concept{link function}

\description{
A parametric link function, i.e. a family of link functions intended for binary data.
}
\usage{
talpha(alpha, verbose = FALSE,
  splineinv = TRUE, eps = 2 * .Machine$double.eps, maxit = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{numeric, must be larger than 0 and smaller than 2.}
  \item{verbose}{logical, warn if truncation occurs when link function or inverse are used.
}
  \item{splineinv}{logical, use spline interpolation for calculation of inverse link?}
  \item{eps}{if splineinv is \code{FALSE}, a Newton-Raphson algorithm is run to calculate the inverse. The argument \code{eps} determines when to terminate this algorithm. Ignored if splineinv is \code{TRUE}.}
  \item{maxit}{maximum number of iterations for Newton-Raphson. Ignored if splineinv is \code{TRUE}.}
}
\value{
An object of class \code{"link-glm"}, see \code{\link[stats]{family}} and \code{\link[stats]{family}}. Intended for use with \code{\link[stats]{glm}}.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}
\examples{
canonical <- binomial(link = talpha(1)) # logit-link
talpha_fam <- function(alpha)binomial(link = talpha(alpha)) # talpha family
## A call to glm might look like this: glm(formula, family = talpha_fam(1.5))
}