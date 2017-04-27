\name{MWAS_heatmap}
\alias{MWAS_heatmap}
\title{Visualize MWAS results as a multiple-phenotype heatmap}
\description{
This function allows visualizing MWAS results generated using multiple phenotypes
as a heatmap. The values of the heatmap are the individual MWAS scores: -log10
p-values (corrected for multiple-testing) adjusted for the direction of the
association. The metabolites are ordered based on hierarchical cluster analysis
of the auto-correlation metabolic matrix.
}
\usage{
MWAS_heatmap (metabo_SE, MWAS_list, alpha_th = 0.05, ncut = 3, ...)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object. See "MWAS_SummarizedExperiment()".}
  \item{MWAS_list}{list of matrices generated with the function "MWAS_stats()".
  The names of the individual matrices must correspond to the phenotype names.
  The dimensions of all matrices must be the same, and consistent with metabo_SE
  dimensions.}
  \item{alpha_th}{numeric value indicating MWAS significance threshold. Metabolites
  with p-value (corrected for multiple-testing) above alpha_th will have a MWAS
  score of 0. }
  \item{ncut}{numeric value indicating where the tree will be cut.}
  \item{...}{other arguments passed to the function "Heatmap()" from the
  ComplexHeatmap package.}
}
\references{
Gu Z, et al. (2016). Complex heatmaps reveal patterns and correlations in
multidimensional genomic data. Bioinformatics, 32, 2847-2849.
}
\value{
A heatmap showing MWAS results generated with multiple phenotypes. The function
also returns a matrix indicating the metabolic clusters.
}
\examples{
## Load data
data(targetMetabo_SE)

## Test for association between diabetes and target_metabolites
T2D_model <- MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         confounder_ids = c("Age", "Gender", "BMI"),
                         assoc_method = "logistic")

## Test for association between BMI and target_metabolites
BMI_model <- MWAS_stats (targetMetabo_SE, disease_id = "BMI",
                         confounder_ids = c("Age", "Gender", "T2D"),
                         assoc_method = "spearman")

## Generate MWAS_list: do not forget the names!
MWAS_list <- list(T2D = T2D_model, BMI = BMI_model)

## Generate heatmap
MWAS_heatmap (targetMetabo_SE, MWAS_list, alpha_th = 0.05)

}
