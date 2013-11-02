\name{phm-class}
\alias{phm-class}
\alias{sroc.phm}
\alias{plot.phm}
\alias{summary.phm}
\alias{print.phm}



\concept{random effects model}
\concept{SROC curve}
\concept{ROC space}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{
Methods for \code{phm} objects.
}
\description{
Objects of the class \code{\link{phm}} are output by the function with the same name. Apart from standard methods the function \code{sroc} provides SROC curves and confidence bands for model fits.
}
\usage{
\method{print}{phm}(x, ...)
\method{summary}{phm}(object, level = 0.95, ...)
\method{sroc}{phm}(fit, fpr = 1:99/100, ...)
\method{plot}{phm}(x, extrapolate = FALSE, confband = TRUE, level = 0.95,
     ylim = c(0,1), xlim = c(0,1), sroclty = 1, sroclwd = 1, 
     confbandlty = 2, confbandlwd = 0.5, ...)
}
\arguments{
  \item{x}{a \code{phm} object.}
  \item{object}{a \code{phm} object.}
\item{fit}{a \code{phm} object.}
  \item{level}{numeric, the confidence level for calculations of confidence intervals (\code{summary}) or confidence bands (\code{plot}).}
  \item{fpr}{numeric, the false positives rates for which to calculate the predicted sensitivities.}
  \item{extrapolate}{logical, should the sroc curve be plotted beyond the observed false positive rates?}
  \item{confband}{logical, should confidence bands be plotted?}
  \item{ylim}{numeric of length 2, which section of the sensitivities to plot?}
  \item{xlim}{numeric of length 2, which section of the false positive rates to plot?}
  \item{sroclty}{integer, line type of the SROC curve}
  \item{sroclwd}{integer, line width of the SROC curve}  
  \item{confbandlty}{integer, line type of the SROC curve's confidence band}
  \item{confbandlwd}{integer, line width of the SROC curve's confidence band}  
  \item{...}{arguments to be passed on to other functions}
}
\details{
The SROC curve is derived from the model formula.  The confidence bands are calculated from the bounds of the confidence interval for the diagnostic accuracy parameter \eqn{\theta}. The parameter and its confidence interval are then also used to calculate the AUC and partial AUC using the formulae
\deqn{
AUC(a,b) = \int_a^bu^\theta\mathrm{d}u = \frac{1}{\theta+1}[b^{\theta+1}-a^{\theta+1}],
}
\deqn{
AUC = AUC(0,1)
}
and
\deqn{
pAUC = \frac{1}{b-a}AUC(a,b),
}
where \eqn{a} is the lower bound of the observed false positive rates and \eqn{b} the upper.
}
\value{
The \code{sroc} function returns a matrix ready for plotting. Each row corresponds to one point in ROC space.}
\references{
 Holling, H., Boehning D., Boehning, W. (2012) \dQuote{Meta-Analysis of Diagnostic Studies based upon SROC-Curves: a Mixed Model Approach using a Proportional Hazards Model.} \emph{Statistical Modelling}, \bold{12}, 347--375.
}

\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{phm}}
}
\examples{
# load data
data(AuditC)
# fit model
fit <- phm(AuditC)
#calculate a SROC curve, but do not plot it
sroc.AuditC <- sroc(fit)
# plot the SROC curve in ROC space as a line
plot(sroc.AuditC, type = "l")
# Fancy version using plot
plot(fit)
}
\keyword{regression}
\keyword{classes}
