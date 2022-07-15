\name{predv_d}
\alias{predv_d}

\title{
Estimation of Distributions of Predictive Values Based on Prevalence Probability Distributions and Pooled Sensitivities and Specificities
}

\description{
Estimation of projected summary predictive values based on a prevalence probability distribution and pooled (meta-analytical) sensitivities and specificities. Probability distributions for negative and positive predictive values are obtained.
}

\usage{
predv_d(x,prop_m,prop_sd,zb=TRUE,n_iter=100000,...)
}

\arguments{
  \item{x}{dataset containing data from the primary studies. It must correspond to any object that can be converted to a data frame with integer variables \code{TP}, \code{FN}, \code{FP} and \code{TN}, alternatively a matrix with column names including \code{TP}, \code{FN}, \code{FP} and \code{TN}. These respectively concern the numbers of true positives, true negatives, false positives, and false negatives for each primary study)}
  \item{prop_m}{mean value of the prevalence probability distribution. It must be stated as a proportion (i.e., as a numeric value between 0 and 1). If both prop_m and prop_sd are not defined, a probability distribution for the prevalence based on available primary studies' data will be computed (see details).}
  \item{prop_sd}{standard-deviation of the prevalence probability distribution. It must be stated as a value between 0 and 1.  If both prop_m and prop_sd are not defined, a probability distribution for the prevalence based on available primary studies' data will be computed (see details).}
  \item{zb}{logical. If TRUE (default), the Zwindermann & Bossuyt approach will be used to generate samples for observed sensitivities and false positive rate (as in SummaryPts function). If FALSE, beta distributions will be obtained based on 95 percent confidence interval bounds of pooled sensitivities and specificities (while this latter approach may not take fully into account the correlation between sensitivity and false positive rate, it may lead to faster results).}
  \item{n_iter}{number of simulations being performed. Default value is 100,000.}
  \item{...}{further arguments to be passed on to \code{\link{predv_d}}.}
}
\details{
The predv_d function projects summary predictive values distributions from (i) a prevalence probability distribution, and (ii) pooled sensitivities and specificities obtained in the context of diagnostic test accuracy meta-analysis using a bivariate random-effects model.
The bivariate random-effects model is equivalent to the hierarchical summary receiver operating characteristic model. By default, a sampling-based approach is used to generate samples for observed sensitivities and false positive rates. From these samples, and based on the prevalences probability distribution being considered, distributions of predictive values will be obtained based on the application of the Bayes theorem.
The prevalence probability distribution can be obtained by providing a value for the mean (argument prop_m) and a value for the standard-deviation (argument prop_sd).
If both prop_m and prop_sd are missing/not defined, a probability distribution for the prevalence based on available primary studies' data will be computed. That is, random-effects meta-analysis of log-transformed prevalences will be performed (using metafor) using data from included primary studies; the pooled results will then be used to obtain the probability distribution for prevalences. This may be a suboptimal option (as there may be considerable heterogeneity, diagnostic accuracy primary studies may not be the best ones to estimate the prevalence of a disease/condition...) compared to user-defined arguments, particularly if good prevalence studies exist.

\bold{Guided example}

The dataset skin_tests contains results from a set of primary studies assessing the accuracy of skin tests for diagnosing penicillin allergy (they are part of the data analysed by Sousa-Pinto et al [2021]).
This dataset contains four columns, displaying - for each primary study - the number of true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN).
Let us now assume that the prevalence of penicillin allergy can be modeled by a probability distribution, having a mean of 0.05 (5 percent) and a standard-deviation of 0.015. Distributions of negative and positive predictive values can be estimated by:

\code{predv_d(x=skin_tests,prop_m=0.05,prop_sd=0.015,zb=TRUE)}

For negative predictive values, we obtain a probability distribution defined by a mean value of 0.96 and a standard-deviation of 0.01 (95 percent credible interval=0.93-0.98).
For positive predictive values, we obtain a probability distribution defined by a mean value of 0.31 and a standard-deviation of 0.12 (95 percent credible interval=0.11-0.57).
Values may differ slightly from the ones just described, as we are dealing with simulation results.

If we had no information on how the prevalence of penicillin allergy could be modeled by a probability distribution, we would opt for solely relying on data provided by included primary studies:

\code{predv_d(x=skin_tests)}

In that case, in addition to the results, we would get an warning message stating that considerable heterogeneity was found when doing meta-analysis of prevalences. Results should be carefully interpreted.
}
\value{An object of class \code{predv_d}, for which some standard methods are available, see \code{\link{predv_d-class}}. Some of the obtainable components include:
\item{SummaryData}{A dataframe displaying the mean, standard-deviation (SD) and percentiles (p) for the probability distribution of the summary negative predictive values ("NPV" row) and positive predictive values ("PPV" row).}
\item{results_pred}{A dataframe displaying the results for all samples.}
}

\references{
Reitsma, J., Glas, A., Rutjes, A., Scholten, R., Bossuyt, P., & Zwinderman, A. (2005).
\dQuote{Bivariate analysis of sensitivity and specificity produces informative summary
measures in diagnostic reviews.} \emph{Journal of Clinical Epidemiology}, \bold{58}, 982--990.

Zwinderman, A., & Bossuyt, P. (2008). \dQuote{We should not pool diagnostic likelihood ratios in  systematic reviews.} \emph{Statistics in Medicine}, \bold{27}, 687--697.

Sousa-Pinto, B., Tarrio, I., Blumenthal, K.G., Azevedo, L.F., Delgado, L., & Fonseca, J.A. (2021). \dQuote{Accuracy of penicillin allergy diagnostic tests: A systematic review and meta-analysis.} \emph{Journal of Allergy and Clinical Immunology}, \bold{147}, 296--308.

Joseph L, Belisle P. (2017). \dQuote{Computing Beta distribution parameters.} [Internet] Accessible at: https://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/BetaParmsFromQuantiles.html
}
\author{
Bernardo Sousa-Pinto <bernardo@med.up.pt>
}

\seealso{
\code{\link{reitsma}},
\code{\link{SummaryPts}},
\code{\link{predv_r}}
}

\examples{
data(skin_tests)
pred_skin_tests <- predv_d(x=skin_tests,prop_m=0.05,prop_sd=0.015,zb=TRUE)
pred_skin_tests
}