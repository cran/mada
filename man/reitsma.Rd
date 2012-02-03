\name{reitsma}
\alias{reitsma}
\alias{reitsma.default}

\concept{random effects model}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{
Fit the bivariate model of Reitsma et al. (2005).
}
\description{
The function fits the bivariate model of Reitsma et al. (2005) that Harbord et al. (2007) have shown to be equivalent to the HSROC of Rutter&Gatsonis (2001). We specify the model as a linear mixed model with known variances of the random effects, similar to the computational approach by Reitsma et al. (2005). Variance components are estimated by restricted maximum likelihood (REML) as a default but ML estimation is available as well.
}
\usage{
reitsma(X, ...)
\method{reitsma}{default}(X = NULL, TP, FN, FP, TN, correction = 0.5, 
    correction.control = "all", REML = TRUE, ...)
}
\arguments{
  \item{X}{any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}.
}
  \item{TP}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{FN}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{FP}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{TN}{vector of integers, ingored if \code{X} is not \code{NULL}}
  \item{correction}{numeric, continuity correction applied if zero cells}
  \item{correction.control}{character, if set to \code{"all"} (the default) the continuity correction is added to the whole data if only one cell in one study is zero. If set to \code{"single"} the correction is only applied to rows of the data which have a zero.}
  \item{REML}{logical, should REML estimation be used? If set to \code{FALSE} ML estimation is performed.}
  \item{...}{arguments to be passed on to other functions, currently ignored}
}
\details{
Parameter estimation makes use of the fact that the mean parameter can be profiled in the likelihood. If optimization of the REML-likelihood with \code{\link{nlm}} fails, it is attempted with the \code{SANN} algorithm in \code{\link{optim}}. More details on the computational approach can be found in Doebler et al. (2012).
}
\value{
An object of the class \code{reitsma} for which many standard methods are available. See \code{\link{reitsma-class}} for details.
}
\references{
Rutter, C., & Gatsonis, C. (2001). \dQuote{A hierarchical regression approach to meta-analysis of
diagnostic test accuracy evaluations.} \emph{Statistics in Medicine}, \bold{20}, 2865--2884.

Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Harbord, R., Deeks, J., Egger, M., Whiting, P., & Sterne, J. (2007). \dQuote{A unification of
models for meta-analysis of diagnostic accuracy studies.} \emph{Biostatistics}, \bold{8}, 239--251.

Doebler, P., Holling, H., Boehning, D. (2012) \dQuote{A Mixed Model Approach to Meta-Analysis of Diagnostic Studies with Binary Test Outcome.} \emph{Psychological Methods}, to appear
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{reitsma-class}}
}
\examples{
data(Dementia)
(fit <- reitsma(Dementia))
summary(fit)
plot(fit)
}

\keyword{regression}
\keyword{multivariate}
