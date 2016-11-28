\name{QC_CV_specNMR}
\alias{QC_CV_specNMR}
\title{Plot NMR spectrum colored based on CV}
\description{
This function allows plotting a reference NMR spectrum colored based on the coefficient
of variation (CV) of each NMR signal. See function "QC_CV()".
}
\usage{
QC_CV_specNMR (metabo_SE, ref_sample,  CV_th = 0.30, xlab = "ppm",
               ylab = "intensity", size_axis = 12, size_lab = 12, xlim = NULL,
               ylim = NULL, xbreaks = waiver (), xnames = waiver (),
               ybreaks = waiver (), ynames = waiver ())
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object generated using the function
  "MWAS_SummarizedExperiment()".}
  \item{ref_sample}{character vector indicating the id of the sample that will
  be used to plot the NMR spectrum.}
  \item{CV_th}{numeric value indicating the CV threshold. NMR signals with CV equal or
  above CV_th will be colored in red.}
  \item{xlab}{character vector specifying a title for the x-axis.}
  \item{ylab}{character vector specifying a title for the y-axis.}
  \item{size_axis}{numeric vector indicating the font size of x- and y-axis labels.}
  \item{size_lab}{numeric vector indicating the font size of x- and y-axis tittles.}
  \item{xlim}{numeric vector containing the minimum and maximum values of the x-axis.
  Notice that ppm is displayed in reverse scaled (e.g. xlim = c(10, 0)).}
  \item{ylim}{numeric vector containing the minimum and maximum values of the y-axis.}
  \item{xbreaks}{numeric vector indicating the positions of the breaks of the x-axis.}
  \item{xnames}{character vector (same length as xbreaks) containing the labels
  of each break of the x-axis.}
  \item{ybreaks}{numeric vector indicating the positions of the breaks of the y-axis.}
  \item{ynames}{character vector (same length as ybreaks) containing the labels
  of each break of the y-axis.}
}
\value{
An NMR spectrum plot colored based on the CV of each NMR signal.
}
\references{
Dumas ME, et al. (2006). Assessment of analytical reproducibility of 1H NMR spectroscopy
based metabonomics for large-scale epidemiological research: the INTERMAP Study.
Analytical Chemistry, 78, 2199-1208.
}
\examples{
## Load data
data(metabo_SE)

## Plot NMR spectrum colored by CV
QC_CV_specNMR (metabo_SE, ref_sample = "QC1", CV_th = 0.30)
QC_CV_specNMR (metabo_SE, ref_sample = "QC1", CV_th = 0.15)

}


