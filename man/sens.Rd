\name{sens}
\alias{sens}
\alias{fpr}
\alias{spec}

\title{Sensitivity, Specificity and False Positive Rate}
\description{
Calculate basic measures of diagnostic accuracy for a number of studies.
}
\usage{
sens(x)
spec(x)
fpr(x)
}
\arguments{
  \item{x}{a data frame with variables including \code{TP}, \code{FN}, \code{FP}, \code{TN}, alternatively a matrix with column names including these.}
}
\details{
These functions are the basic building blocks of many procedures to assess diagnostic accuracy. For a decent summary of  set of primary studies it is better to use \code{\link{madad}}, for graphical summaries \code{\link{crosshair}} and \code{\link{ROCellipse}} are available.
}
\value{
A numeric vector.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{madad}}, \code{\link{crosshair}}, \code{link{ROC.ellipse}}
}
\examples{
data(AuditC)
plot(fpr(AuditC), sens(AuditC), main = "AUDIT-C data on ROC space",
     ylab = "Sensitivity", xlab = "False Positive Rate")
}