\name{JBA_binning}
\alias{JBA_binning}
\title{Binning of pJRES spectra}
\description{
This function performs binning of pJRES spectra using the JBA algorithm (see details).
}
\usage{
JBA_binning (NMR_data, st = 4, ct = 0.85, int = "sum", cm = "pearson",
             ef = 2, merge = TRUE, mt = 0.9)
}
\arguments{
  \item{NMR_data}{numeric matrix containing the NMR data (i.e. NMR peak intensities).
  The columns of the matrix must correspond to the metabolic variables (chemical shifts)
  and the rows to the samples. Column and row names must contain the metabolite IDs
  (i.e chemical shifts) and the sample IDs, respectively.}
  \item{st}{numeric value indicating the minimum bin size.}
  \item{ct}{numeric value indicating the correlation threshold. Bins with average correlation
  below ct will be neglected. This value can be established by comparing the distribution
  of average correlations in a spectral region dominated by electronic noise, and a spectral
  region dominated by metabolic signals. See function "JBA_corDistribution()". }
  \item{int}{character vector indicating the method used to calculate the bin intensity.
  	Possible values are: "sum", "mean", "max", "median". }
  \item{cm}{character vector specifying the correlation method ("pearson" or "spearman").}
  \item{ef}{numeric value establishing the maximum number of upfield or downfield
  variables that can be used to expand the seed. The maximum number of variables that
  can be added on each side of a given seed (size = st) is st*ef.}
  \item{merge}{character constant indicating whether highly correlated (correlation > mt)
  	adjacent bins will be merged (i.e. integrated as a single bin).}
  \item{mt}{numeric value indicating the correlation threshold used to merge adjacent bins.
  	This argument is ignored if merge = FALSE.}
}
\value{
A list containing binned NMR data and information about the bins, as indicated below:

-The first element of the list ("all_clusters") reports the average correlation
of st adjacent NMR variables along the chemical shift axis, using a sliding window of size one.

-The second element ("JBA_seeds") contains the local maxima (i.e. seeds) of the
correlation-based spectrum along the chemical shift axis.

-The fourth element ("JBA_bins_expanded") indicates the bin edges after expanding
the seeds by aggregating upfield and downfield NMR variables.

-The fourth element ("JBA_data") contains the binned NMR data.
}
\details{
JBA ("pJRES Binning Algorithm") is a new binning method designed to extend the
applicability of SRV (Blaise et al., 2009) to pJRES data. The main steps of the
JBA algorithm are described below:

1) The algorithm scans the NMR spectra (from low to high frequencies)
and calculates the average correlation of st adjacent variables, using
a sliding window of size one. This means that a given bin i starts at
the NMR variable i and finishes at NMR variable with i + (st -1).

2) The vector of average correlations can be visualized as pseudo-NMR spectrum,
displaying the average correlation values in the y-axis and the chemical shifts in
the x-axis. This correlation-based spectrum is then scanned to identify local maxima
passing the ct threshold. Each of these local maxima is used as seed that can be
expanded by progressively aggregating upfield and downfield NMR variables, as long as the
following criteria are met: (i) the average correlation of the bin remains equal or
above ct; (ii) for a given upfield variable (vi), correlation (vi, vi+1) needs to
be equal or higher than correlation (vi, vi-1); (iii) for a given downfield variable (vz),
correlation (vz, vz+1) needs to be equal or lower than correlation (vz, vz-1).

3) The intensity of each bin is calculated as the mean, median,
sum or maximum intensity of all variables within the bin. Notice that
due to misalignments/signal overlap, it is possible that a single peak is
split into several bins. These bins can be detected based on a
given correlation threshold and integrated as a single bin.
}
\references{
Blaise, et al. (2009). Statistical recoupling prior to significance testing in
nuclear magnetic resonance based metabonomics. Analytical Chemistry, 81, 6242-6251.
}
\examples{
## Not available.
}
