\name{SummaryPts}
\alias{SummaryPts}
\alias{SummaryPts.default}
\alias{SummaryPts.reitsma}
\alias{print.SummaryPts}
\alias{summary.SummaryPts}


\concept{MCMC}

\title{
Use the Zwindermann & Bossuyt (2008) MCMC procedure to generate summary points (positive and negative likelihood ratio, diagnostic odds ratio) for the Reitsma et al. (2005) bivariate model
}
\description{
Zwindermann & Bossuyt (2008) argue that likelihood ratios should not be pooled by univariate meta-analysis. They propose a sampling based approach that uses the parameters of a fit to the bivariate model (implemented in \code{\link{reitsma}}) to generate samples for observed sensitivities and false positive rates. From these samples the quantities of interest (and their confidence intervals) are estimated.
}
\usage{
SummaryPts(object, ...)
\method{SummaryPts}{default}(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                           n.iter = 10^6, FUN, ...)
\method{SummaryPts}{reitsma}(object, n.iter = 10^6, FUN = NULL, ...)
\method{print}{SummaryPts}(x, ...)
\method{summary}{SummaryPts}(object, level = 0.95, digits = 3, ...)
}

\arguments{
\item{object}{an object for which a method exists}
\item{x}{An object of class \code{SummaryPts}}  
\item{mu}{numeric of length 2, expected to be the mean parameter of a bivariate model}
    \item{Sigma}{2x2 variance covariance matrix, expected to be the matrix representing the standard error of \code{mu}  and the covariance of these two estimates}
  \item{alphasens}{numeric, alpha parameter for the sensitivities. Amounts to logit transformation if set to 1 (the default). See \code{\link{reitsma}}.}
    \item{alphafpr}{numeric, alpha parameter for the false positive rates. Amounts to logit transformation if set to 1 (the default). See \code{\link{reitsma}}.}
      \item{n.iter}{number of samples}
  \item{FUN}{A list of functions with 2 arguments (\code{sens} and \code{fpr}); if set to \code{NULL} in \code{SummaryPts.reitsma}, the positive, negative and inverse negative likelihood ratios are calculated and also the diagnostic odds ratio (DOR). See the example on how to supply other functions.}
  \item{level}{numeric, confidence level for confidence intervals}
  \item{digits}{number of significant digits to display}
  \item{...}{arguments to be passed on to other functions, currently ignored}
}
\details{
Samples are generated from a bivariate normal using \code{\link{rmvnorm}}. Note that the FUN argument
}
\value{
An object of the class \code{SummaryPts} for which \code{print} and \code{summary} methods are available.
}
\references{
Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Zwinderman, A., & Bossuyt, P. (2008). \dQuote{We should not pool diagnostic likelihood ratios
in systematic reviews.}\emph{Statistics in Medicine}, \bold{27}, 687--697.
}
\author{
Philipp Doebler <philipp.doebler@googlemail.com>
}

\seealso{
\code{\link{reitsma}}, \code{\link{talpha}}
}
\examples{
data(Dementia)
(fit <- reitsma(Dementia))
mcmc_sum <- SummaryPts(fit)
mcmc_sum #just the means
summary(mcmc_sum) # 95% CIs by default
summary(mcmc_sum, level = 0.80, digits = 5) ## more digits, smaller CIs

## Supplying other functions

# Example 1: theta parameter of proportional hazards model 
# see "phm" in mada's documentation for details on theta 
theta <- function(sens,fpr){log(sens) / log(fpr)}
theta_sum <- SummaryPts(fit, FUN = list(theta = theta))
summary(theta_sum)
# compare with phm:
summary(phm(Dementia)) # the two estimators almost agree in this example

# Example 2: Youden index
Youden <- function(sens, fpr){sens - fpr}
Youden_sum <- SummaryPts(fit, FUN = list(Youden = Youden))
summary(Youden_sum)
}

