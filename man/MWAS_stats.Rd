\name{MWAS_stats}
\alias{MWAS_stats}
\title{Metabolome-Wide Associations}
\description{
This function tests for association between individual metabolites and a disease
phenotype.}
\usage{
MWAS_stats (metabo_SE, disease_id, confounder_ids = NULL, assoc_method, mt_method = "BH",
            output = "pvalues", CV_metabo = NULL)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object. See "MWAS_SummarizedExperiment()".}
  \item{disease_id}{character vector corresponding to the ID of the response to be
  modeled.}
  \item{confounder_ids}{optional character vector corresponding to the IDs of the
  covariates to be included in the model (e.g. age or gender).}
  \item{assoc_method}{character constant indicating the association method that
  will be used. Possible values for assoc_method are: "pearson" (Pearson correlation),
  "spearman" (Spearman correlation), "kendall" (Kendall correlation), "linear"
  (linear regression) or "logistic" (logistic regression).}
  \item{mt_method}{character constant indicating the multiple-testing correction
  method that will be used. Possible values for mt_method are: "BH" (Benjamini and
  Hochberg), "bonferroni", "holm", "hochberg", "hommel", "BY" (Benjamini and Yekutieli),
  "qvalues", or "none".}
  \item{output}{character constant indicating the output of the function. If
  output = "pvalues", p-values and estimates for each metabolic variable will be
  returned as a matrix. If output = "models", detailed information about the
  statistical model fitted for each metabolic variable will be returned.}
  \item{CV_metabo}{optional numeric vector containing the coefficients of variation
  of the metabolic variables. This vector will be added as an additional column
  of the output matrix.}
}
\value{
By default, a matrix where each row contains the model coefficient estimate
and the p-value obtained for each metabolic variable. When output = "models", the
function returns a list, each list element containing detailed information about
the statistic model fitted for each metabolic variable.
}
\references{
Benjamini Y, Hochberg Y. (1995). Controlling the false discovery rate: a practical
and powerful approach to multiple testing. Journal of the Royal Statistical Society
Series B, 57, 289–300.

Benjamini Y, Yekutieli D. (2001). The control of the false discovery rate in multiple
testing under dependency. Annals of Statistics, 29, 1165–1188.

Dobson AJ. (1990). An Introduction to Generalized Linear Models. London: Chapman and Hall.

Kim, S. (2015). ppcor: An R Package for a Fast Calculation to Semi-partial Correlation
Coefficients. Communications for Statistical Applications and Methods, 22, 665-674.

Holm S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian
Journal of Statistics, 6, 65–70.

Hommel G. (1988). A stagewise rejective multiple test procedure based on a modified
Bonferroni test. Biometrika, 75, 383–386.

Hochberg Y. (1988). A sharper Bonferroni procedure for multiple tests of significance.
Biometrika, 75, 800–803.

Shaffer JP. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46,
561–576.

Storey JD. (2002). A direct approach to false discovery rates. Journal of the Royal
Statistical Society, Series B, 64: 479-498.

Venables WN, Ripley BD. (2002). Modern Applied Statistics with S. New York: Springer.

}
\examples{
## Load data
data(metabo_SE)
data(targetMetabo_SE)

## Test for association between BMI and metabolic_data
BMI_model <- MWAS_stats (metabo_SE, disease_id = "BMI", assoc_method = "spearman",
                         mt_method = "BH", output = "pvalues")

## Test for association between diabetes and target_metabolites (age-gender adjusted)
T2D_model <- MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         confounder_ids = c("Age", "Gender"),
                         assoc_method = "logistic", mt_method = "BY",
                         output = "pvalues")
}
