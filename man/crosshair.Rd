\name{crosshair}
\alias{crosshair}
\alias{crosshair.default}
\title{
Crosshair plot
}
\description{
Produces a crosshair plot or adds such a plot to an existing plot.
}
\usage{
\method{crosshair}{default}(x, correction = 0.5, level = 0.95, method = "wilson",
                      xlim = c(0,1), ylim = c(0,1), length = 0.1, pch = 1, 
                      add = FALSE, suppress = TRUE, ...)
}
\arguments{
  \item{x}{a data frame with variables including \code{TP}, \code{FN}, \code{FP}, \code{TN}, alternatively a matrix with column names including these.}
  \item{correction}{numeric, continuity correction applied to zero cells.}
  \item{level}{numeric, confidence level for the calculations of confidence intervals.
}
\item{method}{character, method used to calculate the confidence intervals for sensitivities, specificities and false positive rates. One of \code{"wald"}, \code{"wilson"}, \code{"agresti-coull"}, \code{"jeffreys"}, \code{"modified wilson"}, \code{"modified jeffreys"}, \code{"clopper-pearson"}, \code{"arcsine"}, \code{"logit"}, \code{"witting"}}  
\item{xlim}{part of ROC space to be plotted}
\item{ylim}{part of ROC space to be plotted}
\item{length}{length of "whiskers" of the crosshair.}
  \item{pch}{Symbol used to plot point estimates. Use \code{pch = ""} to suppress plotting point estimates.}
  \item{add}{logical, should the plot be added to the current plot?}
  \item{suppress}{logical, should the warnings produced by the internal call to \code{madad} be suppressed? Defaults to \code{TRUE}, since only the diagnostic accuracies and their confidence intervals are used in subsequent calculations.}
  \item{\dots}{further arguments passed on to \code{plot}.}
}

\details{Crosshair plots go back to Phillips et al. (2010). Note that for fits of the \code{\link{reitsma}} function a crosshair method is available to plot pooled estimate, see \code{\link{reitsma-class}}.
}

\value{
Besides plotting, the function returns an invisible \code{NULL}.
}
\references{
Phillips, B., Stewart, L.A., & Sutton, A.J. (2010). \dQuote{'Cross hairs' plots for diagnostic meta-analysis.} \emph{Research Synthesis Methods}, \bold{1}, 308--315.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{ROCellipse}}, \code{\link{reitsma-class}}

}
\examples{
data(AuditC)
crosshair(AuditC)
}

\keyword{aplot}
\keyword{hplot}
