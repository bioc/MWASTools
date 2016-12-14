\name{MWAS_skylineNMR}
\alias{MWAS_skylineNMR}
\title{Visualize MWAS results in an NMR-skyline plot}
\description{
This function generates a 2-panel figure showing the results from "MWAS_stats()"
applied to NMR data. The upper panel shows an NMR-skyline plot (comparable to a
GWAS-Manhattan plot), where the chemical shifts are displayed along the x-axis and
the -log10 p-values (sign-adjusted for the direction of the association) are
displayed on the y-axis. The lower panel shows an NMR spectrum colored according to
MWAS results.
}
\usage{
MWAS_skylineNMR (metabo_SE, MWAS_matrix, ref_sample, alpha_th = 0.05, output = "all",
                 xlab ="ppm", ylab1 = "sign*log(pFDR)", ylab2 = "intensity", pch = 20,
                 marker_size = 1, scale_color = c("black", "cornflowerblue", "red"),
                 size_lab = 12, size_axis = 12, xlim = NULL, ylim1 = NULL,
                 ylim2 = NULL, guide_type = "legend", xbreaks = waiver(),
                 xnames = waiver(), ybreaks1 = waiver(), ybreaks2 = waiver(),
                 ynames1 = waiver(), ynames2 = waiver())
}
\arguments{
  \item{metabo_SE}{SummarizedExperiment object generated using the function
  MWAS_SummarizedExperiment().}
  \item{MWAS_matrix}{numeric matrix resulting from the function "MWAS_stats()".}
  \item{ref_sample}{character vector indicating the id of the sample that will
  be used to plot the NMR spectrum.}
  \item{alpha_th}{numeric value indicating the significance threshold.}
  \item{output}{character constant indicating the outcome of the function ("skyline",
  "spectrum" or "all"). If outcome = "all", both the skyline and the spectrum with be plotted
  in a 2-panel plot.}
  \item{xlab}{character vector specifying a title for the x-axis.}
  \item{ylab1}{character vector specifying a title for the y-axis of the upper panel.}
  \item{ylab2}{character vector specifying a title for the y-axis of the lower panel.}
  \item{pch}{value specifying the symbol used to represent each ppm value in the
  skyline plot. To see all possible symbols, check "plot()" options.}
  \item{marker_size}{numeric value indicating the size of the symbol used to represent
  each ppm value in the skyline plot.}
  \item{scale_color}{character vector corresponding to the 3-color scale that
  will be used to represent the association results. The first color of the scale
  indicates "no change", the second color indicates "downregulation", and the third
  color indicates "upregulation".}
  \item{size_lab}{numeric value indicating the font size of x- and y-axis titles.}
  \item{size_axis}{numeric value indicating the font size of x- and y-axis labels.}
  \item{xlim}{numeric vector containing the minimum and maximum values of the x-axis.
  Notice that ppm is displayed in reverse scale (e.g. xlim = c(5, 2)).}
  \item{ylim1}{numeric vector containing the minimum and maximum values of the y-axis
  for the upper panel.}
  \item{ylim2}{numeric vector containing the minimum and maximum values of the y-axis
  for the lower panel.}
  \item{guide_type}{character constant indicating the guide ("legend" or "none") that
  will be added to the plots.}
  \item{xbreaks}{numeric vector indicating the positions of the breaks of the x-axis.}
  \item{xnames}{character vector (same length as xbreaks) containing the labels
  of each break of the x-axis.}
  \item{ybreaks1}{numeric vector indicating the positions of the breaks of the y-axis
  for the upper panel.}
  \item{ybreaks2}{numeric vector indicating the positions of the breaks of the y-axis
  for the lower panel.}
  \item{ynames1}{character vector (same length as ybreaks1) containing the labels
  of each break of the y-axis for the upper panel.}
  \item{ynames2}{character vector (same length as ybreaks2) containing the labels
  of each break of the y-axis for the lower panel.}
}
\value{
By default, a plot with 2 panels, the upper panel showing an NMR-skyline plot and
the lower panel showing an NMR spectrum colored based on MWAS results.
}
\references{
Elliott P, et al. (2015). Urinary metabolic signatures of human adiposity. Science
Translational Medicine, 7, 285ra62.
}
\examples{
## Load data
data(metabo_SE)

## Test for association between BMI and metabolic_data
BMI_model <- MWAS_stats (metabo_SE, disease_id = "BMI", assoc_method = "spearman",
                         output = "pvalues")

## Create skyline plots
MWAS_skylineNMR (metabo_SE, BMI_model, ref_sample = "QC1")
MWAS_skylineNMR (metabo_SE, BMI_model, ref_sample = "QC1", pch = "*", marker_size = 3)
}


