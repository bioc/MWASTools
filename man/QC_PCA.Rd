\name{QC_PCA}
\alias{QC_PCA}
\title{Principal Component Analysis}
\description{
This function performs PCA on a matrix of metabolic data and returns the results as
an object of class "prcomp". When quality control (QC) samples are available,
"QC_PCA()" can be used to assess the stability and reproducibility of the dataset.}
\usage{
QC_PCA (metabo_SE, scale = FALSE, center = TRUE,...)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object generated using the function
  "MWAS_SummarizedExperiment()".}
  \item{scale}{logical constant indicating whether the metabolic variables will be
  scaled to have unit variance before the analysis. For more details, check "prcomp()".}
  \item{center}{logical constant indicating whether the metabolic variables will be
  shifted to be zero-centered before the analysis. For more details, check "prcomp()".}
  \item{...}{other arguments passed to "prcomp()".}
}
\value{
A list with class "prcomp". For more details, check "prcomp()".
}
\references{
Mardia K, et al. (1979). Multivariate Analysis, London: Academic Press.
}
\examples{
## Load data
data(metabo_SE)
data(targetMetabo_SE)

\donttest{
## PCA model using all metabolic data
PCA_model <- QC_PCA (metabo_SE)
}
## PCA model using target metabolites
PCA_subset <- QC_PCA (targetMetabo_SE)
}
