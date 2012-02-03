\name{cochran.Q}
\alias{cochran.Q}
\title{
Cochran's Q statistic
}
\description{
Given estimates from primary studies and the weights of the single studies calculate Cochran's Q as a measure of heterogeneity.
}
\usage{
cochran.Q(x, weights)
}
\arguments{
  \item{x}{numeric, typically a vector of effect sizes like (log-)OR}
  \item{weights}{numeric, see Details}
}
\details{
In fixed effects settings the weights are often inverse proportional to the variances of the primary studies. Cochran's Q is known to have low power to detect heterogeneity.
}
\value{
A named vector of length 3. First element is \code{Q} followed by the \code{p-value} and the degrees of freedom.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}
