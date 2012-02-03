\name{madauni}
\alias{madauni}

\concept{meta-analysis}
\concept{diagnostic meta-analysis}
\concept{random effects model}
\concept{Mantel-Haenszel estimation}


\title{
Meta-Analyisis of univariate measures of diagnostic accuracy
}
\description{
The classic strategy to meta-analysis of diagnostic accuracy data is to pool a univariate measure of accuracy like the diagnostic odds ratio, the positive likelihood ratio or the negative likelihood ratio. For fixed effect estimation a Mantel-Haenszel estimator is implemented and for random effect estimation a DerSimonian-Laird estimator is available.
}
\usage{
madauni(x, type = "DOR", method = "DSL", suppress = TRUE, ...)
}
\arguments{
  \item{x}{any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}.}
  \item{type}{character, what effect size should be pooled? Either \code{"DOR"}, \code{"posLR"} or \code{"negLR"}.}
  \item{method}{character, method of estimation. Either \code{"MH"} or \code{"DSL"}.}
  \item{suppress}{logical, should warnings produced by the internal call to \code{\link{madad}} be suppressed?}
  \item{...}{further arguments to be passed on to \code{\link{madad}}, for example \code{correction.control}.}
}
\details{
The Mantel-Haenszel method performs fixed effect estimation of the effect sizes. For the DOR the variance of this estimator is calculated according to Robins et al. (1986) and for the likelihood ratios the variance is based on Greenland et al. (1985).

The DerSimonian-Laird method performs a random effects meta-analysis. For this \eqn{\tau^2}{tau^2}, the variance of the log-transformed effect size (DOR, positive or negative likelihood ratio) is calculated by the DerSimonian and Laird (1986) method. The confidence interval for \eqn{\tau^2}{tau^2} is derived by inverting the Q-Test of Viechtbauer (2007).
}
\value{An object of class \code{madauni}, for which some standard methods are available, see \code{\link{madauni-class}}}

\references{
DerSimonian, R. and Laird, N. (1986). \dQuote{Meta-analysis in clinical trials.} \emph{Controlled clinical trials}, \bold{7}, 177--188.

Greenland, S. and Robins, J.M. (1985). \dQuote{Estimation of a Common Effect Parameter from Sparse Follow-Up Data.} \emph{Biometrics}, \bold{41}, 55--68. 

Robins, J. and Greenland, S. and Breslow, N.E. (1986). \dQuote{A general estimator for the variance of the Mantel-Haenszel odds ratio.} \emph{American Journal of Epidemiology}, \bold{124}, 719--723.

Viechtbauer, W. (2007). \dQuote{Confidence intervals for the amount of heterogeneity in meta-analysis.} \emph{Statistics in Medicine}, \bold{26}, 37--52.
}

\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}
\note{
Performing univariate meta-analysis of diagnostic studies can not be recommended anymore now that bivariate methods are available, at least not if a reasonable number of primary studies is availabel. The package \code{mada} provides this functionality for exploratory purposes and for meta-analysis of a small number of studies.
}

\seealso{
\code{\link{madauni-class}}
}
\examples{
data(AuditC)
AuditC.uni <- madauni(AuditC)
summary(AuditC.uni)
}
