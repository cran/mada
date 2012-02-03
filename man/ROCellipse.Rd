\name{ROCellipse}
\alias{ROCellipse}
\alias{ROCellipse.default}
\title{
Confidence Regions on ROC space
}
\description{
Plot individual confidence regions for the estimate from each primary study on ROC space or add such regions to an existing plot.
}
\usage{
\method{ROCellipse}{default}(x, correction = 0.5, level = 0.95, 
           xlim = c(0, 1), ylim = c(0, 1), method = "wilson", 
           pch = 1, add = FALSE, corr = 0, suppress = TRUE, 
           ellipsecol = "grey", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data frame with variables including \code{TP}, \code{FN}, \code{FP}, \code{TN}, alternatively a matrix with column names including these.}
  \item{correction}{numeric, continuity correction applied to zero cells.}
  \item{level}{numeric, confidence level for the calculations of confidence intervals.
}
  \item{xlim}{numeric of length 2, which portion of ROC space should be plotted? All reasonable values should be within (0,1).}
  \item{ylim}{numeric of length 2, which portion of ROC space should be plotted? All reasonable values should be within (0,1).}
\item{method}{character, method used to calculate the confidence intervals for sensitivities, specificities and false positive rates. One of \code{"wald"}, \code{"wilson"}, \code{"agresti-coull"}, \code{"jeffreys"}, \code{"modified wilson"}, \code{"modified jeffreys"}, \code{"clopper-pearson"}, \code{"arcsine"}, \code{"logit"}, \code{"witting"}}
\item{pch}{Symbol used to plot point estimates. Use \code{pch = ""} to suppress plotting point estimates.}
\item{add}{logical, should the plot be added to the current plot?}
\item{corr}{numeric or character, the correlation assumed in the calculation of the confidence ellipsoids on logit-ROC space. If set to \code{"logit"}, the correlation of the logit-transformed sensitivities and false positive rates will be used in the correlations. See details for further explanation.}
\item{suppress}{logical, should the warnings produced by the internal call to \code{madad} be suppressed? Defaults to \code{TRUE}, since only the diagnostic accuracies and their confidence intervals are used in subsequent calculations.}
\item{ellipsecol}{The color used for plotting the ellipses.}
\item{\dots}{further arguments passed on to \code{plot}.}
}
\details{
The confindence regions are ellipses on logit-ROC space, hence the name of the function. The standard deviations underlying confidence intervals for the sensitivities and false positive rates are used to determine the scale of the ellipses on logit-ROC space. These ellipses get backtransformed to ROC space and plotted. As a default no correlation is assumed on logit-ROC space.

The objects of class \code{\link{reitsma}} have their own \code{ROCellipse} method to add a confidence region for the pooled estimate, see \code{\link{reitsma-class}}.
}
\value{
Besides plotting an invisble \code{NULL} is returned.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}
\seealso{
\code{\link{crosshair}}, \code{\link{reitsma-class}}
}
\examples{
data(AuditC)
ROCellipse(AuditC)
}
\keyword{aplot}
\keyword{hplot}

