\name{JBA_corDistribution}
\alias{JBA_corDistribution}
\title{Setting the ct threshold for JBA}
\description{
This function compares the distribution of correlations between st adjacent variables
in a spectral region dominated by noise and a spectral region dominated by metabolic
signals. This function can be used to set the ct threshold for JBA binning.
}
\usage{
JBA_corDistribution(NMR_data, st = 4, cm = "pearson", metabo_range = c(3.50, 3.96),
                    noise_range = c(9.72, 9.99), color_scale = c("lightcoral", "honeydew3"))
}
\arguments{
  \item{NMR_data}{numeric matrix containing the NMR data (i.e. NMR peak intensities).
  The columns of the matrix must correspond to the metabolic variables (chemical shifts)
  and the rows to the samples. Column and row names must contain the metabolite IDs
  (i.e chemical shifts) and the sample IDs, respectively.}
  \item{st}{numeric value indicating the minimum bin size.}
  \item{cm}{character vector specifying the correlation method ("pearson" or "spearman").}
  \item{metabo_range}{numeric vector indicating the limits of a spectral region dominated
  by metabolic signals.}
  \item{noise_range}{numeric vector indicating the limits of a spectral region dominated by noise.}
  \item{color_scale}{character vector indicating color of the metabolic curve (first value),
  and the noise curve (second value).}
  }
\value{
A plot comparing the distribution of average correlations between st adjacent variables
in a spectral region dominated by metabolic signals (metabo_range) and in a spectral
region dominated by electronic noise (noise_range). The suggested ct value corresponds
to the correlation coefficient where the cumulative proportion of noise clusters is 1.
}
\examples{
## Not available.
}


