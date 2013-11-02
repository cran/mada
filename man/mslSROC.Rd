\name{mslSROC}
\alias{mslSROC}

\concept{diagnostic meta-analysis}
\concept{meta-analysis}
\concept{SROC curve}

\title{
Plot the Moses-Shapiro-Littenberg SROC curve
}
\description{
The approach to SROC curve modeling is described in the paper of Moses, Shapiro and Littenberg (1993). It is considered outdated and is included in \code{mada} so that users can reproduce older results and compare different SROC curves.
}
\usage{
mslSROC(data = NULL, subset=NULL,
  TP="TP", FN="FN", FP="FP", TN="TN", 
  fpr = NULL, extrapolate = FALSE, 
  correction = 0.5, correction.control = "all",
  add = FALSE, lty = 1, lwd = 1, col = 1, ...)
}
\arguments{
  \item{data}{any object that can be converted to a data frame with integer variables for observed frequencies of true positives, false negatives, false positives and true negatives. The names of the variables  are provided by the arguments \code{TP}, \code{FN}, \code{FP} and \code{TN} (see their defaults). Alternatively the data can be a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}. If no \code{data} is specified, the function will check the \code{TP}, \code{FN}, \code{FP} and \code{TN} arguments.
}
  \item{TP}{character or integer: name for vector of integers that is a variable of \code{data} or a vector of integers. If \code{data} is not \code{NULL}, names are expected, otherwise integers are.}
  \item{FN}{character or integer: name for vector of integers that is a variable of \code{data} or a vector of integers. If \code{data} is not \code{NULL}, names are expected, otherwise integers are.}
  \item{FP}{character or integer: name for vector of integers that is a variable of \code{data} or a vector of integers. If \code{data} is not \code{NULL}, names are expected, otherwise integers are.}
  \item{TN}{character or integer: name for vector of integers that is a variable of \code{data} or a vector of integers. If \code{data} is not \code{NULL}, names are expected, otherwise integers are.}
\item{subset}{the rows of \code{data} to be used as a subset in all calculations. If \code{NULL} (the default) then the complete data is considered.}
\item{fpr}{Points between 0 and 1 on which to draw the SROC curve. Should be tightly spaced. If set to \code{NULL}, the default, it will be the vector of numbers \code{0.01, 0.02, ..., 0.99} and is truncated if the \code{extrapolate} argument is \code{FALSE}.}
\item{extrapolate}{logical, should the SROC curve be extrapolated beyond the region where false positive rates are observed?}
  \item{correction}{numeric, continuity correction applied if zero cells}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
  \item{add}{logical, should the SROC curve be added to an existing plot?}
  \item{lty}{line type, see \code{\link{lines}}.}
  \item{lwd}{line width, see \code{\link{lines}}.}
  \item{col}{color of SROC, see \code{\link{lines}}.}
  \item{...}{arguments to be passed on to plotting functions.}
}
\details{
Details are found in the paper of Moses, Shapiro and Littenberg (1993).
}
\value{
Besides plotting the SROC, an \code{\link{invisible}} list is returned which contains the parameters of the SROC.
}
\references{
Moses L.E., Shapiro D., & Littenberg B. (1993) \dQuote{Combining independent studies of a diagnostic test into a summary ROC curve: data-analytic approaches and some additional considerations.} \emph{Statistics in Medicine}, \bold{12}, 1293--1316.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{reitsma-class}}, \code{\link{talpha}}, \code{\link{SummaryPts}}
}
\examples{
## First Example
data(Dementia)
ROCellipse(Dementia)
mslSROC(Dementia, add = TRUE) # Add the MSL-SROC to this plot

## Second Example
# Make a fancy plot and look at the coefficients
msl_Dementia <- mslSROC(Dementia, col = 3, lwd = 3, lty = 3)
msl_Dementia$A2 # intercept on logit SROC space
msl_Dementia$B2 # slope on logit SROC space
}

