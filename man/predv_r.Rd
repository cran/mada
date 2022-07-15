\name{predv_r}
\alias{predv_r}

\concept{predictive values}
\concept{diagnostic meta-analysis}
\concept{meta-analysis}

\title{
Estimation of Distributions of Predictive Values Based on Prevalence Ranges and Pooled Sensitivities and Specificities
}
\description{
Estimation of projected summary predictive values based on a prevalence range and pooled (meta-analytical) sensitivities and specificities. A probability distribution for the negative and positive predictive values are obtained for each prevalence value within a predetermined range.
}
\usage{
predv_r(x,prop_min,prop_max,zb=TRUE,n_iter=100000,...)
}

\arguments{
  \item{x}{dataset containing data from the primary studies. It must correspond to any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}. These respectively concern the numbers of true positives, true negatives, false positives, and false negatives for each primary study)}
  \item{prop_min}{minimum prevalence value being considered. It must be stated as a proportion (i.e., as a numeric value between 0 and 1). If both prop_min and prop_max are not defined, a prevalence range based on available primary studies' data will be computed (see details).}
  \item{zb}{logical. If TRUE (default), the Zwindermann & Bossuyt approach will be used to generate samples for observed sensitivities and false positive rate (as in SummaryPts function). If FALSE, beta distributions will be obtained based on 95 percent confidence interval bounds of pooled sensitivities and specificities (while this latter approach may not take fully into account the correlation between sensitivity and false positive rate, it may lead to faster results).}
  \item{prop_max}{maximum prevalence value being considered. It must be stated as a proportion (i.e., as a numeric value between 0 and 1). If both prop_min and prop_max are not defined, a prevalence range based on available primary studies' data will be computed (see details).}
  \item{n_iter}{number of simulations being performed. Default value is 100,000.}
  \item{...}{further arguments to be passed on to \code{\link{predv_r}}.}
}

\details{
The predv_r function projects summary predictive values from (i) a prevalence range, and (ii) pooled sensitivities and specificities obtained in the context of diagnostic test accuracy meta-analysis using a bivariate random-effects model.
The bivariate random-effects model is equivalent to the hierarchical summary receiver operating characteristic model. By default, a sampling-based approach is used to generate samples for observed sensitivities and false positive rate. From these samples, and for each prevalence value within the range being considered, distributions of predictive values will be obtained based on the application of the Bayes theorem.
The prevalence range can be user-defined, by providing a value for the minimum (argument prop_min) and a value for the maximum value of that range (argument prop_max).
If both prop_min and prop_max are missing/not defined, a prevalence range based on available primary studies' data will be computed. That is, the lowest and highest frequency of patients with disease/condition across included primary studies will be considered. This may be a suboptimal option compared to user-defined arguments, particularly if good prevalence studies are available.

\bold{Guided example}

The dataset skin_tests contains results from a set of primary studies assessing the accuracy of skin tests for diagnosing penicillin allergy (they are part of the data analysed by Sousa-Pinto et al [2021]).
This dataset contains four columns, displaying - for each primary study - the number of true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN).
Let us assume that the prevalence of penicillin allergy ranges between 0.01 and 0.10 (1 and 10 percent). Pooled negative and positive predictive values can be estimated by:

\code{predv_r(x=skin_tests,prop_min=0.01,prop_max=0.15,zb=TRUE)}

The results indicate that the point estimates for the negative predictive value range between 0.88 (prevalence=0.15) and 0.99 (prevalence=0.01).
For the positive predictive value, point estimates range between 0.09 (prevalence=0.01) and 0.59 (prevalence=0.15), although uncertainty is particularly high for the latter estimate (95 percent credible interval=0.36-0.80).
Values may differ slightly from the ones just described, as we are dealing with simulation results.

If we had no information on how the prevalence range of penicillin allergy, we would opt for solely relying on data provided by included primary studies:

\code{pred_skin_tests1 <- predv_r(x=skin_tests)}

}

\value{
An object of class \code{predv_r}, for which some standard methods are available, see \code{\link{predv_r-class}}. Some of the obtainable components include:
\item{NPV}{A dataframe displaying the mean, standard-deviation (SD) and percentiles (p) for the probability distribution of negative predictive values for each prevalence value within the defined range.}
\item{PPV}{A dataframe displaying the mean, standard-deviation (SD) and percentiles (p) for the probability distribution of positive predictive values for each prevalence value within the defined range.}
}

\references{
Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Zwinderman, A., & Bossuyt, P. (2008). \dQuote{We should not pool diagnostic likelihood ratios in  systematic reviews.} \emph{Statistics in Medicine}, \bold{27}, 687--697.

Sousa-Pinto, B., Tarrio, I., Blumenthal, K.G., Azevedo, L.F., Delgado, L., & Fonseca, J.A. (2021). \dQuote{Accuracy of penicillin allergy diagnostic tests: A systematic review and meta-analysis.} \emph{Journal of Allergy and Clinical Immunology}, \bold{147}, 296--308.
}

\author{
Bernardo Sousa-Pinto <bernardo@med.up.pt>
}

\seealso{
\code{\link{reitsma}},
\code{\link{SummaryPts}},
\code{\link{predv_d}}
}

\examples{
data(skin_tests)
pred_skin_tests <- predv_r(x=skin_tests,prop_min=0.01,prop_max=0.15,zb=TRUE)
pred_skin_tests
}