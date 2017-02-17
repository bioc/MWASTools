\docType{data}
\name{metabo_SE}
\alias{metabo_SE}
\title{NMR plasma metabolic profiles dataset}
\format{SummarizedExperiment}
\description{
  This SummarizedExperiment object contains the following information:

  -An assay matrix containing the 1H NMR profiles (1.60 - 0.80 ppm) of 506 plasma samples
  from the FGENTCARD cohort and 10 identical quality control (QC) samples. The QC samples
  were prepared from a representative pool of the experimental samples, and were injected
  regularly throughout the run to ensure analytical reproducibility.

  -A data.frame containing clinical information (age, gender, type II diabetes status and BMI)
  and sample class (i.e. experimental sample or QC sample) information for each sample row
  in the assay matrix.
}
\usage{data(metabo_SE)}
\value{SummarizedExperiment}

