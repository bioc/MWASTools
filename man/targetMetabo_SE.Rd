\docType{data}
\name{targetMetabo_SE}
\alias{targetMetabo_SE}
\title{Target NMR metabolites dataset}
\format{SummarizedExperiment}
\description{
  This SummarizedExperiment object contains the following information:

  -An assay matrix containing the levels of 8 target 1H NMR metabolites (lactate,
  3-hydroxy-butyrate, leucine, valine, isoleucine,acetate, alanine and 1,5-anhydroglucitol)
  across the experimental samples and the quality control (QC) samples.

  -A data.frame containing clinical information (age, gender, type II diabetes status and BMI)
  and sample class (i.e. experimental sample or QC sample) information for each sample row
  in the assay matrix.
}
\usage{data(targetMetabo_SE)}
\value{SummarizedExperiment}
