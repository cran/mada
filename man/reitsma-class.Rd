\name{reitsma-class}
\alias{reitsma-class}
\alias{sroc}
\alias{mcsroc}
\alias{sroc.reitsma}
\alias{mcsroc.reitsma}
\alias{ROCellipse.reitsma}
\alias{crosshair.reitsma}
\alias{plot.reitsma}
\alias{summary.reitsma}
\alias{print.reitsma}

\concept{random effects model}
\concept{SROC curve}
\concept{ROC space}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{
Methods for \code{reitsma} objects.
}
\description{
Objects of the class \code{\link{reitsma}} are output by the function with the same name. Apart from standard methods the functions \code{sroc}, \code{mcsroc} and \code{ROCellipse} provide SROC curves and confidence regions for fits.
}
\usage{
\method{print}{reitsma}(x, ...)
\method{summary}{reitsma}(object, level = 0.95, ...)
\method{sroc}{reitsma}(fit, fpr = 1:99/100, ...)
\method{mcsroc}{reitsma}(fit, fpr = 1:99/100, replications = 10000, lambda = 100, ...)
\method{ROCellipse}{reitsma}(x, level = 0.95, add = FALSE, pch = 1, ...)
\method{crosshair}{reitsma}(x, level = 0.95, length = 0.1, pch = 1, ...)
\method{plot}{reitsma}(x, extrapolate = FALSE, plotsumm = TRUE, level = 0.95, 
     ylim = c(0,1), xlim = c(0,1), pch = 1, sroclty = 1, sroclwd = 1, ...)
}
\arguments{
  \item{x}{a \code{reitsma} object.}
  \item{object}{a \code{reitsma} object.}
\item{fit}{a \code{reitsma} object.}
  \item{level}{numeric, the level for calculations of confidence intervals (\code{summary}) or regions (\code{ROCellipse})}
  \item{fpr}{numeric, the false positives rates for which to calculate the predicted sensitivities}
  \item{replications}{integer, the number of replications for the Monte-Carlo SROC curve}
  \item{lambda}{numeric, the parameter lambda of the Monte-Carlo run, see details}
  \item{add}{logical, should the confidence region be added to the current plot? If set to \code{FALSE} a matrix of points of the ellipse is returned}
  \item{extrapolate}{logical, should the SROC curve be plotted beyond the observed false positive rates?}
  \item{plotsumm}{logical, should the summary pair of sensitivity and false positive rate together with its confidence region be plotted?}
  \item{length}{positve numeric, length of the "whiskers" of the crosshairs.}
  \item{ylim}{numeric of length 2, which section of the sensitivities to plot?}
  \item{xlim}{numeric of length 2, which section of the false positive rates to plot?}
  \item{pch}{integer, symbol for the pair of mean sensitivity and false positive rate}
  \item{sroclty}{integer, line type of the SROC curve}
  \item{sroclwd}{integer, line width of the SROC curve}  
  \item{...}{arguments to be passed on to other functions}
}
\details{
The confidence regions of \code{ROCellipse} are first calculated as ellipses on logit-ROC space, so the back-transformed regions that are output are not necessarily ellipses. The Monte-Carlo SROC curves are generated from random samples from the fitted model and a \code{\link{lowess}} smooth through them is output. Many computational details are to be found in Doebler et al. (2012).
}
\value{
\code{sroc} a  matrix ready for plotting. Each row corresponds to one point in ROC space. \code{mcsroc} returns a \code{\link{lowess}} smooth. \code{ROCellipse} returns a list, the first element being a matrix of points in ROC space that delimit the confidence region and the second is the point estimate of the pair of sensitivity and false positive rate in ROC space.
}
\references{
Doebler, P., Holling, H., Boehning, D. (2012) \dQuote{A Mixed Model Approach to Meta-Analysis of Diagnostic Studies with Binary Test Outcome.} \emph{Psychological Methods}, to appear
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{reitsma}}
}
\examples{
# load data
data(Dementia)
# fit model
fit <- reitsma(Dementia)
# calculate a confidence region but do not plot it
cr.Dementia <- ROCellipse(fit)
#calculate a SROC curve
sroc.Dementia <- sroc(fit)
# plot the confidence region in ROC space as a line
plot(cr.Dementia$ROCellipse, type = "l", xlim = c(0,1), ylim = c(0,1))
# add the point estimate of the mean
points(cr.Dementia$fprsens)
# add the SROC curve
lines(sroc.Dementia)

}
\keyword{regression}
\keyword{classes}


