\name{MWAS_KEGG_pathways}
\alias{MWAS_KEGG_pathways}
\title{Map metabolites into KEGG pathways}
\description{
This function allows mapping the metabolites of interest detected by MWAS
analysis onto the KEGG pathways. The function also exports a network file and an
attribute file which can be imported into Cytoscape to visualize the results
as a pathway-based metabolic network.
}
\usage{
MWAS_KEGG_pathways(metabolites, MWAS_matrix = NULL, file_name = "KeggPaths")
}
\arguments{
  \item{metabolites}{character vector containing the KEGG IDs of the metabolites
  of interest detected by MWAS. The order of the metabolite IDs in this vector must
  match the order in MWAS_matrix. Compound KEGG IDs can be obtained using the
  function "MS_keggFinder()" from the MetaboSignal package.}
  \item{MWAS_matrix}{numeric matrix generated with the function "MWAS_stats()".
  It can also be a submatrix containing only the significant metabolites, generated
  with the function "MWAS_filter()". }
  \item{file_name}{character vector that allows customizing the name of the
  exported files.}
}
\references{
Rodriguez-Martinez A, et al. (2017).MetaboSignal: a network-based approach for
topological analysis of metabotype regulation via metabolic and signaling pathways.
Bioinformatics, 33, 773-775.

Shannon P, et al. (2003). Cytoscape: a software environment for integrated models of
biomolecular interaction networks. Genome Research, 13, 2498-2504.

Tenenbaum D. (2017). KEGGREST: Client-side REST access to KEGG. R package.
}
\value{
A six-column matrix indicating the KEGG pathways where each metabolite was mapped.
The results are formatted as a six-column matrix containing the following
information: metabolite KEGG ID (column 1), metabolite name (column 2), pathway
KEGG ID (column 3), pathway name (column 4), pathway class (column 5), pathway
organism (i.e. "Human"/"Not_human") (column 6).

The function also exports a network file ("KeggPaths_NetworkFile.txt") and an
attribute file ("KeggPaths_AttributeFile.txt") that can be imported into
Cytoscape to visualize the results as a network. The attribute file allows
customizing the metabolites of interest based on a score reflecting the degree
of association with the phenotype under study (i.e. log10(pvalue) adjusted for
the sign of the association).
}
\examples{
\donttest{## Test for association between diabetes and target_metabolites
T2D_model <- MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         confounder_ids = c("Age", "Gender", "BMI"),
                         assoc_method = "logistic")

## Select the metabolites of interest and get their corresponding KEGG IDs
T2D_model_subset = T2D_model[1:5, ]
kegg_metabolites = c("cpd:C00186", "cpd:C01089", "cpd:C00123", "cpd:C00183",
                     "cpd:C00407")

## Map metabolites into KEGG pathways
kegg_pathways = MWAS_KEGG_pathways(metabolites = kegg_metabolites,
                                   MWAS_matrix = T2D_model_subset)
}
}
