\name{QC_CV}
\alias{QC_CV}
\title{Calculate coefficients of variation}
\description{
This function calculates the coefficient of variation (CV) ( |sd/mean| ) of each
metabolic feature across the quality control (QC) samples. The CV distribution is
represented in a histogram. This function can be used to assess the reproducibility
of individual metabolic features. Notice that CV = 0.30 and CV = 0.15 are the thresholds
established by the FDA guidelines for biomarker discovery and quantification, respectively.
}
\usage{
QC_CV (metabo_SE, CV_th = 0.30, plot_hist = TRUE, hist_bw = 0.005,
       hist_col = "moccasin", size_lab = 12, size_axis = 12)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object. See "MWAS_SummarizedExperiment()".}
  \item{CV_th}{numeric value indicating the CV threshold.}
  \item{plot_hist}{logical constant indicating whether a histogram showing CV
  distribution will be plotted.}
  \item{hist_bw}{numeric value indicating histogram bin width.}
  \item{hist_col}{character string indicating the color to be used to fill the histogram
  bars.}
  \item{size_lab}{numeric value indicating the font size of x- and y-axis titles.}
  \item{size_axis}{numeric value indicating the font size of x- and y-axis labels.}
}
\value{
A numeric vector containing the CV of each metabolic feature and a histogram
showing CV distribution. In the histogram, CVs above 1 are set to 1.
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
metabo_CV <-  QC_CV (metabo_SE)
metabo_CV2 <-  QC_CV (metabo_SE, hist_bw = 0.008, hist_col = "lightblue")
}

