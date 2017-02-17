\name{QC_CV_scatterplot}
\alias{QC_CV_scatterplot}
\title{Plot MS-based scatter plot colored based on CV}
\description{
This function creates a MS-based scatter plot (rt vs mz) where the metabolic
features are colored based on their coefficient of variation (CV). See "QC_CV()".
}
\usage{
QC_CV_scatterplot (rt, mz, CV_metabo, CV_th = 0.30, xlab = "rt",
                   ylab = "mz", pch = 20, marker_size = 1, xlim = NULL,
                   ylim = NULL, size_axis = 10, size_lab = 10)
}
\arguments{
  \item{rt}{numeric vector of retention time values.}
  \item{mz}{numeric vector of mz values.}
  \item{CV_metabo}{numeric vector containing the CV of each metabolic feature. The
  length of this vector should be consistent with the length of rt and mz.}
  \item{CV_th}{numeric value indicating the CV threshold. NMR signals with CV equal or
  above CV_th will be colored in red.}
  \item{xlab}{character vector specifying a title for the x-axis.}
  \item{ylab}{character vector specifying a title for the y-axis.}
  \item{pch}{value specifying the symbol used to represent each MS feature in the
  scatter plot. To see all possible symbols, check "plot()" options.}
  \item{marker_size}{numeric value indicating the size of the symbol used to represent
  each metabolic feature in the scatter plot.}
  \item{xlim}{numeric vector containing the minimum and maximum values of the x-axis.}
  \item{ylim}{numeric vector containing the minimum and maximum values of the y-axis.}
  \item{size_axis}{numeric value indicating the font size of x- and y-axis title.}
  \item{size_lab}{numeric value indicating the font size of x- and y-axis labels.}
}
\value{
A MS-based scatter plot where MS features are represented according on CV.
}
\examples{
## Load data
data(MS_data)
rt <- MS_data[, 1]
mz <- MS_data[, 2]

## Simulate CV values
CV_metabo <- runif(length(rt), 0.05, 0.31)

## MS-based scatter plot
QC_CV_scatterplot(rt, mz, CV_metabo)
QC_CV_scatterplot(rt, mz, CV_metabo, xlim = c(0, 10))
QC_CV_scatterplot(rt, mz, CV_metabo, CV_th = 0.15)
}

