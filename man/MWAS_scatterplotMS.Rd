\name{MWAS_scatterplotMS}
\alias{MWAS_scatterplotMS}
\title{Visualize MWAS results in MS-based scatter plot}
\description{
This function creates a MS-based scatter plot (rt vs mz) based on the output
from "MWAS_stats()". MS-features are filtered according to a given significance
threshold and only significant features are represented in the scatter plot. The
color of the marker representing each significant MS feature indicates the direction
of the association (i.e. downregulation or upregulation) and the size of the marker
indicates the strength of the association (i.e. -log10(p-value)).
}
\usage{
MWAS_scatterplotMS (rt, mz, MWAS_matrix, alpha_th = 0.05, xlab = "rt",
                    ylab = "mz", pch = 20, scale_color = c("cornflowerblue", "red"),
                    xlim = NULL, ylim = NULL, size_axis = 10, size_lab = 10,
                    legend_position = "bottom")
}
\arguments{
  \item{rt}{numeric vector of retention time values.}
  \item{mz}{numeric vector of mz values.}
  \item{MWAS_matrix}{numeric matrix resulting from the function "MWAS_stats()".
  The dimensions of this matrix must be consistent with the length of rt and mz}
  \item{alpha_th}{numeric value indicating the significance threshold. Only
  variables with p-value (corrected for multiple-testing) below alpha_th will
  be plotted.}
  \item{xlab}{character vector specifying a title for the x-axis.}
  \item{ylab}{character vector specifying a title for the y-axis.}
  \item{pch}{value specifying the symbol used to represent each MS feature in the
  scatter plot. To see all possible symbols, check "plot()" options.}
  \item{scale_color}{character vector corresponding to the 2-color scale that
  will be used to represent the association results. The first color of the scale
  indicates "downregulation", and the second color indicates "upregulation".}
  \item{xlim}{numeric vector containing the minimum and maximum values of the x-axis.}
  \item{ylim}{numeric vector containing the minimum and maximum values of the y-axis.}
  \item{size_axis}{numeric value indicating the font size of x- and y-axis title.}
  \item{size_lab}{numeric value indicating the font size of x- and y-axis labels.}
  \item{legend_position}{character vector indicating the position of the legend:
  "top", "bottom", "right", "left", "none".}
}
\value{
A MS-based scatter plot where MS features are represented according to MWAS_results.
}
\examples{
## Load data
data(MS_data)
rt <- MS_data[, 1]
mz <- MS_data[, 2]

## Simulate MWAS data
set.seed(100)
estimates <- runif(length(rt), -1, 1)
pvalues <- rbeta(length(estimates), 0.5, 1)
pFDR <- p.adjust(pvalues, method = "BH")
MWAS_matrix <- cbind(estimates, pvalues, pFDR)

## MS-based scatter plot
MWAS_scatterplotMS(rt, mz, MWAS_matrix)
MWAS_scatterplotMS(rt, mz, MWAS_matrix, alpha_th = 0.01)
MWAS_scatterplotMS(rt, mz, MWAS_matrix, alpha_th = 0.01,
                  scale_color = c("yellow", "blue"))
}

