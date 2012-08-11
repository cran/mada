\name{AUC}
\alias{AUC}
\alias{auc}
\alias{AUC.default}
\alias{AUC.reitsma}
\alias{AUC.phm}

\title{
Area under the curve (AUC)
}
\description{
Calculates the area under the curve given a function or a fitted model.
}
\usage{
\method{AUC}{default}(x, fpr = 1:99/100, ...)
\method{AUC}{phm}(x, level = 0.95, ...)
\method{AUC}{reitsma}(x, fpr = 1:99/100, ...)
}
\arguments{
  \item{x}{a function with range and domain in ROC space (default method) or an object of class \code{\link{phm}} or \code{\link{reitsma}}.}
  \item{fpr}{numeric vector, points on which the (S)ROC curve is evaluated}
  \item{level}{numeric, confidence level for the calculations of confidence intervals.
}
  \item{\dots}{further arguments, currently not used.}
}

\details{The area under the curve is calculated using the trapezoidal rule. The argument \code{fpr} is the grid on which the (S)ROC curve is evaluated. In many cases the default grid will contain points on which the SROC curve can only be calculated by extrapolation; however if only a subinterval is specified a \emph{partial AUC} is calculated and the AUC value might differ substantially.

For \code{\link{phm}} objects the AUC and its confidence interval is calculated analytically, for \code{\link{reitsma}} objects a call to the default method is invoked.
}

\value{
An object of the class \code{AUC} which is really a list with component \code{AUC} and an optional component \code{ci}, which is currently only available from the \code{AUC} method for \code{\link{phm}} ojects.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\examples{
data(AuditC)
AUC(phm(AuditC))
}
