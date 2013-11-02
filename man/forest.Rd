\name{forest}
\alias{forest}
\alias{forestmada}
\alias{forest.madad}
\alias{forest.madauni}



\title{
Forest plot for univariate measures
}
\description{
Produce a forest plot. Includes graphical summary of results if applied to output of suitable model-fitting function. \code{forest} methods for \code{\link{madad}} and \code{\link{madauni}} objects are provided.
}
\usage{
\method{forest}{madad}(x, type = "sens", log = FALSE, ...)
\method{forest}{madauni}(x, log = TRUE, ...)
forestmada(x, ci, plotci = TRUE, main = "Forest plot", xlab = NULL,
          digits = 2L,  snames = NULL, subset = NULL, pch = 15, 
          cex = 1, cipoly = NULL, polycol = NA, ...) 
}
\arguments{
\item{x}{an object for which a \code{forest} method exists or (in the case of \code{foresmada}) a vector of point estimates.}
\item{ci}{numeric matrix, each row corresponds to a confidence interval (the first column being the lower bound and the second the upper).}
\item{plotci}{logical, should the effects sizes and their confidence intervals be added to the plot (as text)?}
\item{main}{character, heading of plot.}
\item{xlab}{label of x-axis.}
\item{digits}{integer, number of digits for axis labels and confidence intervals.}
\item{snames}{character vector, study names. If \code{NULL}, generic study names are generated.}
\item{subset}{integer vector, allows to study only a subset of studies in the plot. One can also reorder the studies with the help of this argument.}
\item{pch}{integer, plotting symbol, defaults to a small square. Also see \code{\link{plot.default}}.}
\item{cex}{numeric, scaling parameter for study names and confidence intervals.}
\item{cipoly}{logical vector, which confidence interval should be plotted as a polygon? Useful for summary estimates. If set to \code{NULL}, regular confidence intervals will be used.}
\item{polycol}{color of the polygon(s), passed on to \code{\link{polygon}}. The default value of \code{NA} implies no color.}
\item{type}{character, one of \code{sens}, \code{spec}, \code{negLR}, \code{posLR} or \code{DOR}.}
\item{log}{logical, should the log-transformed values be plotted?}
\item{...}{arguments to be passed on to \code{forestmada} and further on to other plotting functions}
}
\details{
Produces a forest plot to graphically assess heterogeneity. Note that \code{forestmada} is called internally, so that the \code{...} argument can be used to pass on arguments to this function; see the examples. 
}
\value{
Returns and invisible \code{NULL}.
}

\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{madad}}, \code{\link{madauni}}
}
\examples{
data(AuditC)

## Forest plot of log DOR with random effects summary estimate
forest(madauni(AuditC))

## Forest plot of negative likelihood ratio (no log transformation)
## color of the polygon: light grey 
## draw the individual estimate as filled circles
forest(madauni(AuditC, type = "negLR"), 
       log = FALSE, polycol = "lightgrey", pch = 19)

## Paired forest plot of sensitivities and specificities
## Might look ugly if device region is too small
old.par <- par()
AuditC.d <- madad(AuditC)

plot.new()
par(fig = c(0, 0.5, 0, 1), new = TRUE)
forest(AuditC.d, type = "sens", xlab = "Sensitivity")
par(fig = c(0.5, 1, 0, 1),  new = TRUE)
forest(AuditC.d, type = "spec", xlab = "Specificity")

par(old.par)

## Including study names
## Using Letters as dummies
forest(AuditC.d, type = "spec",  xlab = "Specificity",
      snames = LETTERS[1:14])

}

\keyword{hplot}
