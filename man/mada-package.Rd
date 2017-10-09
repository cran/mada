\name{mada-package}
\alias{mada-package}
\alias{mada}
\docType{package}
\title{
Meta-Analysis of diagnostic accuracy studies
mada
}
\description{
This package provides functions for diagnostic meta-analysis. Next to basic analysis and visualization the bivariate Model of Reitsma et al. (2005) that is equivalent to the HSROC of Rutter&Gatsonis (2001) can be fitted. A new approach based to diagnostic meta-analysis of Holling et al. (2012) is also available. Standard methods like summary, plot and so on are provided.
}
\details{
\tabular{ll}{
Package: \tab mada\cr
Type: \tab Package\cr
Version: \tab 0.5.8\cr
Date: \tab 2017-10-06\cr
License: \tab GPL-2\cr
}
The package provides tools for the meta-analysis of diagnostic accuracy data. For this the number true positives (TP), false negatives (FN), true negatives (TN) and false postives (FP) for each study must be known. The package can fit the bivariate model of Reitsma et al (2005), a bivariate random effects model. This model has been shown by Harbord et al. (2007) to be equivalent to the HSROC proposed by Rutter & Gatsonis (2001). We approach  this model as a linear mixed effects model to avoid the complications of non-linear mixed effects model. The main function to fit such model is \code{\link{reitsma}} and standard methods are available for the output of this function.
}
\author{
Author and Maintainer: Philipp Doebler
}

\references{
Rutter, C., & Gatsonis, C. (2001). \dQuote{A hierarchical regression approach to meta-analysis of
diagnostic test accuracy evaluations.} \emph{Statistics in Medicine}, \bold{20}, 2865--2884.

Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Harbord, R., Deeks, J., Egger, M., Whiting, P., & Sterne, J. (2007). \dQuote{A unification of
models for meta-analysis of diagnostic accuracy studies.} \emph{Biostatistics}, \bold{8}, 239--251.
}

\keyword{package}
\keyword{regression}
\seealso{
\code{\link{reitsma}}
}
