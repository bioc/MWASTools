\name{MWAS_bootstrapping}
\alias{MWAS_bootstrapping}
\title{MWAS bootstrap resampling}
\description{
This function generates bootstrap replicates (non-parametric resampling) of a model
testing for association between a given metabolite and a disease phenotype,
and calculates the confidence intervals of model coefficients.
}
\usage{
MWAS_bootstrapping (metabo_SE, metabolite_id, disease_id, confounder_ids = NULL,
                    assoc_method, iterations = 10000)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object generated using the function
  "MWAS_SummarizedExperiment()".}
  \item{metabolite_id}{character vector corresponding to the id of the metabolite
  to be modeled.}
  \item{disease_id}{character vector corresponding to the id of the response
  to be modeled.}
  \item{confounder_ids}{optional character vector corresponding to the ids of
  covariates to be included in the model (e.g. age or gender).}
  \item{assoc_method}{character constant indicating the association method that
  will be used. Possible values for assoc_method are: "pearson" (pearson correlation),
  "spearman" (spearman correlation), "kendall" (kendall correlation), "linear"
  (linear regression) or "logistic" (logistic regression).}
  \item{iterations}{numeric value indicating the number of bootstrap replicates}
}
\value{
A list with 3 elements, each list element reporting the following information:
i) object of class "boot"; ii) summary of the previous object; iii) 95-confidence
interval of the metabolite model coefficient. For more details, check the function
"boot()" from the "boot" package.
}
\references{
Davison AC, Hinkley, DV. (1997). Bootstrap Methods and Their Application. Cambridge
University Press.
}
\examples{
## Load data
data(targetMetabo_SE)

## Bootstrap model testing for association between diabetes (T2D) and 3OH-butyrate
MWAS_bootstrapping (targetMetabo_SE, metabolite_id = "3OH-butyrate", disease_id = "T2D",
                    assoc_method = "logistic", iterations = 1000)
}

