\name{rsSROC}
\alias{rsSROC}

\concept{diagnostic meta-analysis}
\concept{meta-analysis}
\concept{SROC curve}

\title{
Plot the Ruecker-Schumacher (2010) SROC curve
}
\description{
Assuming that a weighted Youden index is maximized in all primary studies, the Ruecker-Schumacher approach estimates individual ROC curves and then averages them.
}
\usage{
rsSROC(data = NULL, subset=NULL,
  TP="TP", FN="FN", FP="FP", TN="TN", 
  lambda = "from_bivariate",                  
  fpr = NULL, extrapolate = FALSE, plotstudies = FALSE,
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
\item{lambda}{numeric or \code{"from_bivariate"}, the weight of the weighted Youden index. Must be between 0 and 1. If set to \code{"from_bivariate"}, the \code{\link{reitsma}} function is used to calculate lambda from the data.}
\item{fpr}{Points between 0 and 1 on which to draw the SROC curve. Should be tightly spaced. If set to \code{NULL}, the default, it will be the vector of numbers \code{0.01, 0.02, ..., 0.99} and is truncated if the \code{extrapolate} argument is \code{FALSE}.}
\item{extrapolate}{logical, should the SROC curve be extrapolated beyond the region where false positive rates are observed?}
\item{plotstudies}{logical, should the ROC curves for the individual studies be added to the plot? The plot will become crowded if set to \code{TRUE}.}
  \item{correction}{numeric, continuity correction applied if zero cells}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
  \item{add}{logical, should the SROC curve be added to an existing plot?}
  \item{lty}{line type, see \code{\link{lines}}.}
  \item{lwd}{line width, see \code{\link{lines}}.}
  \item{col}{color of SROC, see \code{\link{lines}}.}
  \item{...}{arguments to be passed on to plotting functions.}
}
\details{
Details are found in the paper of Ruecker and Schumacher (2010).
}
\value{
Besides plotting the SROC, an \code{\link{invisible}} list is returned which contains the parameters of the SROC.
}
\references{
Ruecker G., & Schumacher M. (2010) \dQuote{Summary ROC curve based on a weighted Youden index for selecting an optimal cutpoint in meta-analysis of diagnostic accuracy.} \emph{Statistics in Medicine}, \bold{29}, 3069--3078.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com> 
Original code kindly supplied by G. Ruecker.
}

\seealso{
\code{\link{reitsma-class}}, \code{\link{talpha}}, \code{\link{SummaryPts}}
}
\examples{
## First Example
data(Dementia)
ROCellipse(Dementia)
rsSROC(Dementia, add = TRUE) # Add the RS-SROC to this plot

## Second Example
# Make a crowded plot and look at the coefficients
rs_Dementia <- rsSROC(Dementia, col = 3, lwd = 3, lty = 3, 
                      plotstudies = TRUE)
rs_Dementia$lambda
rs_Dementia$aa # intercepts of primary studies on logit ROC space
rs_Dementia$bb # slopes 
}

