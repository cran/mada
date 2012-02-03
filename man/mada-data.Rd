\name{mada-data}
\docType{data}
\alias{mada-data}
\alias{AuditC}
\alias{Dementia}
\alias{IAQ}
\alias{SAQ}

\title{Diagnostic accuracy data}

\description{Four data frames with diagnostic accuracy data from binary test outcomes.}

\usage{
data("AuditC")
data("Dementia")
data("IAQ")
data("SAQ")
}


\format{Four data frames with 4 variables.
  \describe{
  \item{TP}{integer. number of true positives}  
  \item{FN}{integer. number of false negatives}  
  \item{FP}{integer. number of false positives}  
  \item{TN}{integer. number of true negatives}  
  }
}


\details{
The \code{AuditC} data is from Kriston et al. (2008). The \code{Dementia} from Mitchell (2009) and the \code{SAQ} and \code{IAQ} data are subsets from the data in Patrick et al. (1994)
}

\source{
Kriston, L., H{\"o}elzel, L., Weiser, A., Berner, M., & Haerter, M. (2008).\dQuote{ Meta-analysis: Are 3
Questions Enough to Detect Unhealthy Alcohol Use?} \emph{ Annals of Internal Medicine},
\bold{149}, 879--888.

Mitchell, A. (2009). \dQuote{A meta-analysis of the accuracy of the mini-mental state examination
in the detection of dementia and mild cognitive impairment.} \emph{Journal of Psychiatric
Research}, \bold{43}, 411--431.

Patrick, D., Cheadle, A., Thompson, D., Diehr, P., Koepsell, T., & Kinne, S. (1994). \dQuote{The
validity of self-reported smoking: a review and meta-analysis.} \emph{American Journal of
Public Health}, \bold{84}, 1086--1093.
}

\keyword{datasets}
