\name{MWAS_SummarizedExperiment}
\alias{MWAS_SummarizedExperiment}
\title{Create a SummarizeExperiment object}
\description{
This function formats the metabolic and clinical data into a SummarizedExperiment
object.}

\usage{
MWAS_SummarizedExperiment(metabo_matrix, clinical_matrix, sample_type)
}
\arguments{
  \item{metabo_matrix}{numeric matrix containing the metabolic data (e.g. NMR peak
  intensities or metabolite concentrations). The columns of the matrix must correspond
  to the metabolic variables and the rows to the samples. Column and row names must
  contain the metabolite IDs (e.g. chemical shifts for NMR data) and the sample IDs,
  respectively.}
  \item{clinical_matrix}{numeric matrix containing the clinical data (e.g. age, gender).
  The columns of the matrix must correspond to the phenotypic variables and the
  rows to the samples. Column and row names must contain the phenotype IDs and
  the sample IDs, respectively. For samples without clinical data (e.g. quality
  control (QC) samples), NA values must be used. }
  \item{sample_type}{numeric vector indicating sample type (i.e. experimental sample
  or QC sample). The vector must be coded as follows: experimental sample = 0,
  QC sample = 1. If QC samples are not available, all the elements of this vector
  must be 0.}
}
\value{
A SummarizedExperiment object.
}
\references{
Morgan M, et al. (2016). SummarizedExperiment: SummarizedExperiment container.
R package.
}
\examples{
## Load data
data(metabo_SE)

\donttest{## Get metabolic_data, clinical_data, and sample_type
library(SummarizedExperiment)
metabolic_data = t(assays(metabo_SE)$metabolic_data)
clinical_data = as.matrix(colData(metabo_SE)[, -5])
sample_type = as.vector(colData(metabo_SE)[, 5])

## Reconstruct SummarizedExperiment
data_SE = MWAS_SummarizedExperiment(metabolic_data, clinical_data,
                                    sample_type)
}
}
