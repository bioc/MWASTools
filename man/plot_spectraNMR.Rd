\name{plot_spectraNMR}
\alias{plot_spectraNMR}
\title{Plot NMR spectra}
\description{
This function generates an NMR spectra plot, with the chemical shifts displayed
along the x-axis, and the peak intensities displayed on the y-axis.
}
\usage{
plot_spectraNMR (metabo_SE, type = "l", lty = 1, xlab ="ppm",
                 ylab ="intensity", xlim = NULL, ...)
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object. See "MWAS_SummarizedExperiment()".}
  \item{type}{character vector indicating the type of plot for each row of metabo_matrix.
  For all possible types, see "plot()".}
  \item{lty}{character vector of line types. For all possible types, see "plot()".}
  \item{xlab}{character vector specifying a title for the x-axis.}
  \item{ylab}{character vector specifying a title for the y-axis.}
  \item{xlim}{numeric vector containing the minimum and maximum values of the x axis.
  Notice that ppm is displayed in reverse scale (e.g. xlim = c(10, 0)).}
  \item{...}{other arguments passed to "matplot()".}
}
\value{
An NMR spectra plot.
}
\examples{
## Load data
data(metabo_SE)

## Plot first 2 spectra
plot_spectraNMR (metabo_SE[, 1:2])
plot_spectraNMR (metabo_SE[, 1:2], xlim = c(1.03, 0.85), main = "NMR spectra")
}

