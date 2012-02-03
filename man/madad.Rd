\name{madad}
\alias{madad}
\alias{madad-class}
\alias{print.madad}


\concept{diagnostic meta-analysis}
\concept{meta-analysis}
\concept{descriptive statistics}

\title{
Descriptive statistics for meta-analysis of diagnostic accuracy
}
\description{
Given the frequencies of true positives, false negative, false positives and true negatives from primary diagnostic studies \code{madad} calculates various summary statistics. Apart from sensitivities, specificities and false positive rates the function also calculates the diagnostic odds ratio (DOR) and the positve and negative likelihood ratios, together with their respective confidence intervals. Also two hypothesis tests are calculated: one testing the equality of the sensitivities and the same for the false positive rates. 
}
\usage{
madad(x = NULL, TP, FN, FP, TN, level = 0.95, correction = 0.5, 
  correction.control = "all", method = "wilson", yates = TRUE, ...)

\method{print}{madad}(x, digits = 3, ...)
}
\arguments{
  \item{x}{any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}.}
  \item{TP}{vector of integers, ingored if \code{X} is not \code{NULL}.}
  \item{FN}{vector of integers, ingored if \code{X} is not \code{NULL}.}
  \item{FP}{vector of integers, ingored if \code{X} is not \code{NULL}.}
  \item{TN}{vector of integers, ingored if \code{X} is not \code{NULL}.}
  \item{correction}{numeric, continuity correction applied to zero cells.}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
\item{level}{numeric, confidence level for the calculations of confidence intervals.}
\item{method}{character, method used to calculate the confidence intervals for sensitivities, specificities and false positive rates. One of \code{"wald"}, \code{"wilson"}, \code{"agresti-coull"}, \code{"jeffreys"}, \code{"modified wilson"}, \code{"modified jeffreys"}, \code{"clopper-pearson"}, \code{"arcsine"}, \code{"logit"}, \code{"witting"}}
\item{yates}{logical, should a Yates correction be used for testing the equality of sensitivities and specificities?}
\item{digits}{integer, to what decimal place is the output to be rounded?}
\item{\dots}{further arguments to be passed on the other funtions (currently none).}
}
\details{
All calculations are performed using the continuity corrected cell counts, so if there are zero cells, the sensitivities and specificities not equal to 1. This can be avoided by setting \code{correction.control} to  \code{"none"}.

The test for the equality of sensitivities and its counterpart for the specificities is based on \code{\link{prop.test}}. This function will occasionally output warnings.

}
\value{
An object of class \code{madad} which is essentially a list with the following components:
\item{sens}{A list of two components, \code{sens} (the sensitivities) and \code{sens.ci} the confidence intervals (a matrix with 2 columns).}
\item{spec}{A list of two components, \code{spec} (the specificities) and \code{spec.ci} the confidence intervals (a matrix with 2 columns).}
\item{fpr}{A list of two components, \code{fpr} (the false positive rates) and \code{fpr.ci} the confidence intervals (a matrix with 2 columns).}
\item{sens.htest}{An object of class \code{htest}.}
\item{spec.htest}{An object of class \code{htest}.}
\item{DOR}{A list of two components, \code{DOR} the diagnostic odds ratios and \code{DOR.ci} the confidence intervals (a matrix with 2 columns).}
\item{posLR}{A list of two components, \code{posLR} the positive likelihood ratios and \code{posLR.ci} the confidence intervals (a matrix with 2 columns).}
\item{negLR}{A list of two components, \code{negLR} the negative likelihood ratios and \code{negLR.ci} the confidence intervals (a matrix with 2 columns).}
\item{cor_sens_fpr}{numeric, the correlation of the sensitivities and false-positive rates.}
\item{level}{numeric}
\item{method}{character}
\item{names}{character vector, if the main argument of \code{madad} is a data frame with a variable \code{names} these names are stored here.}
\item{nobs}{integer, number of primary studies.}
\item{data}{data frame, with columns \code{TP}, \code{FN}, \code{FP} and \code{TN}.}
\item{data.name}{character, name of the main argument.}
\item{correction}{numeric}
\item{correction.control}{character}
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{madauni}}
}
\examples{
data(AuditC)
AuditC.d <- madad(AuditC)
print(AuditC.d, digits = 2) #round everything to 2 digits
}
\keyword{univar}
