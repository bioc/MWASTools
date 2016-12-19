\name{CV_filter}
\alias{CV_filter}
\title{Filter metabolic data by CV}
\description{
This function allows filtering a matrix of metabolic variables based on the coefficient
of variation (CV) of each variable across the quality control (QC) samples. See also
function "QC_CV()".
}
\usage{
CV_filter(metabo_SE, CV_metabo, CV_th = 0.30)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object. See "MWAS_SummarizedExperiment()".}
  \item{CV_metabo}{numeric vector containing the CVs of the metabolic variables.
  See function "QC_CV()".}
  \item{CV_th}{numeric value indicating the CV threshold. Only features with CV below
  CV_th will be retained in the matrix.}
}
\value{
A SummarizedExperiment object containing the CV-filtered metabolic_data.
}
\references{
Dumas ME, et al. (2006). Assessment of analytical reproducibility of 1H NMR spectroscopy
based metabonomics for large-scale epidemiological research: the INTERMAP Study.
Analytical Chemistry, 78, 2199-1208.
}
\examples{
## Load data
data(metabo_SE)

## Calculate CVs
CV_metabo  <-  QC_CV (metabo_SE)

## Filter metabolic_data by CV
metabo_CVfiltered <- CV_filter(metabo_SE, CV_metabo, CV_th = 0.30)
metabo_CVfiltered2 <- CV_filter(metabo_SE, CV_metabo, CV_th = 0.15)
}
