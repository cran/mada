\name{phm}
\alias{phm}
\alias{phm.default}

\concept{proportional hazards model}
\concept{random effects model}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}
\concept{SROC curve}

\title{
Diagnostic Meta-Analysis with the proportional hazards model approach of Boehning et.al (2012)
}
\description{
The function fits the model of Boehning" et al. (2012). The adjusted profile maximum likelihood estimator (APMLE) is implemented for homogeneity and heterogeneity of primary studies.
}
\usage{
phm(x, ...)
\method{phm}{default}(x=NULL, TP, FN, FP, TN, 
    correction = 0.5, correction.control = "all", 
    hetero = TRUE, estimator = "APMLE", l = 100, ...)
}
\arguments{
\item{x}{any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}.
}
  \item{TP}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{FN}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{FP}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{TN}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{correction}{numeric, continuity correction applied if zero cells}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
  \item{hetero}{logical, should heterogeneity of studies be assumed? Will fit model for homogeneity otherwise.}
  \item{estimator}{character, determines estimator used. Currently only \code{APMLE} is available.}
  \item{l}{interger, number of iterations for fixed point algorithm}
  \item{...}{arguments passed on to other functions (currently not used)}
}
\details{
The model of Boehning et al. (2012) assumes that the relationship between false positive rates \eqn{u} and and sensitivities \eqn{p} can be described by
\deqn{u^\theta = p,}
where \eqn{\theta} is the diagnostic accuracy parameter. If homogeneity of the studies can be assumed, \eqn{\theta} is estimated as a fixed effect. Under heterogeneity a random effect with variance \eqn{\tau^2} describes the variation of the diagnostic accuracy parameter in the population of studies. Since the error of each observed \eqn{\theta} depends only on the sample size and \eqn{\theta} the model has only one parameter in the case of homogeneity and two parameters under heterogeneity, making it suitable for diagnostic meta-analysis with low sample size. Estimation proceeds by a fixed point algorithm derived from the adjusted profile likelihood. More details on the computational approach can be found in Boehning et al. (2012).
}
\value{
An object of the class \code{phm} for which many standard methods are available. See \code{\link{phm-class}} for details.
}
\references{
Holling, H., Boehning, W., Boehning D.(2012) \dQuote{Meta-Analysis of Diagnostic Studies based upon SROC-Curves: a Mixed Model Approach using a Proportional Hazards Model.} \emph{Statistical Modelling -- an International Journal}
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>, 
Walailuck Boehning (original implementation of estimation algorithm)
}

\seealso{
\code{\link{phm-class}}
}
\examples{
data(AuditC)
(fit <- phm(AuditC))
summary(fit)
plot(fit)
}

\keyword{regression}
