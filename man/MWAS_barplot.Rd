\name{MWAS_barplot}
\alias{MWAS_barplot}
\title{Visualize MWAS results in a bar plot}
\description{
This function creates a bar plot based on the output from "MWAS_stats()". This
function is designed to visualize MWAS results in the case of discrete metabolic
variables (e.g. target GC/MS metabolites).
}
\usage{
MWAS_barplot(MWAS_matrix, alpha_th = 0.05, width = NULL,
             scale_color = c("darkgray", "cornflowerblue", "firebrick1"),
             legend_labs =c("unchanged","downregulated","upregulated"),
             ylab = "sign*log(pFDR)", size_yaxis = 12, size_ylab = 12,
             size_names = 10, angle_names = 45, sort = TRUE)
}
\arguments{
  \item{MWAS_matrix}{numeric matrix resulting from the function "MWAS_stats()".}
  \item{alpha_th}{numeric value indicating the significance threshold.}
  \item{width}{numeric value indicating bar width.}
  \item{scale_color}{character vector corresponding to the 3-color scale that
  will be used to represent the association results. The first color of the scale
  indicates "no change", the second color indicates "downregulation", and the third
  color indicates "upregulation".}
  \item{legend_labs}{character vector containing the legend labels, according to
  scale_color.}
  \item{ylab}{character vector specifying a title for the y-axis.}
  \item{size_yaxis}{numeric value indicating the font size of y-axis title.}
  \item{size_ylab}{numeric value indicating the font size of y-axis labels.}
  \item{size_names}{numeric value indicating the font size of the metabolite ids
  displayed on the x-axis.}
  \item{angle_names}{numeric value indicating the angle in which the metabolite
  ids will be displayed on the x-axis.}
  \item{sort}{logical constant indicating whether the metabolites will be sorted
  based on MWAS results.}
}
\value{
A bar plot.
}
\examples{
## Load data
data(targetMetabo_SE)

## Test for association between diabetes and target_metabolites
T2D_model <- MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         confounder_ids = c("Age", "Gender", "BMI"),
                         assoc_method = "logistic")

## Bar plot
MWAS_barplot(T2D_model)
MWAS_barplot(T2D_model, width = 0.7) # change bar width
MWAS_barplot(T2D_model, width = 0.7, angle_names = 90)
}

