\name{MWAS_KEGG_shortestpaths}
\alias{MWAS_KEGG_shortestpaths}
\title{Build a shortest-path metabolic subnetwork}
\description{
This function allows calculating the shortest paths between the metabolites of
interest detected by MWAS analysis, and representing them as a network. The
function also generates a network file and an attribute file, which can be easily
imported into Cytoscape to visualize the network.
}
\usage{
MWAS_KEGG_shortestpaths(network_table, metabolites, MWAS_matrix = NULL,
                        type = "all", distance_th = "Inf", names = TRUE,
                        file_name = "KeggSP")
}
\arguments{
  \item{network_table}{four-column matrix where each row represents an edge
  between two nodes. See function "MWAS_KEGG_network()".}
  \item{metabolites}{character vector containing the KEGG IDs of the metabolites
  of interest detected by MWAS. The order of the metabolite IDs in this vector must
  match the order in MWAS_matrix. Compound KEGG IDs can be obtained using the
  function "MS_keggFinder()" from the MetaboSignal package.}
  \item{MWAS_matrix}{numeric matrix generated with the function "MWAS_stats()".
  It can also be a submatrix containing only the significant metabolites, generated
  with the function "MWAS_filter()". }
  \item{type}{character constant indicating whether all shortest paths
  (type = "all") or a single shortest path (type = "first") will be considered
  when there are several shortest paths between a given source metabolite and a
  given target metabolite.}
  \item{distance_th}{establishes a shortest path length threshold. Only shortest
  paths with length below this threshold will be included in the network.}
  \item{names}{logical scalar indicating whether the metabolite KEGG IDs
  will be transformed into common metabolite names.}
  \item{file_name}{character vector that allows customizing the name of the
  exported files.}
}
\references{
Csardi G, Nepusz T. (2006). The igraph software package for complex network
research. InterJournal, Complex Systems, 1695.

Posma JM, et al. (2014). MetaboNetworks, an interactive Matlab-based toolbox
for creating, customizing and exploring sub-networks from KEGG. Bioinformatics,
30, 893-895.

Rodriguez-Martinez A, et al. (2017).MetaboSignal: a network-based approach for
topological analysis of metabotype regulation via metabolic and signaling pathways.
Bioinformatics, 33, 773-775.

Shannon P, et al. (2003). Cytoscape: a software environment for integrated models of
biomolecular interaction networks. Genome Research, 13, 2498-2504.
}
\value{
A four-column matrix where each row represents an edge connecting two metabolites
(from metabolite in column 1 to metabolite in column 2). The reactions involved
in each metabolic conversion as well as the reaction type (i.e. reversible
or irreversible) are reported in the third and fourth columns, respectively.
This network can be visualized in R using the igraph package or similar packages.

The function also exports a network file ("KeggSP_NetworkFile.txt") and an
attribute file ("KeggSP_AttributeFile.txt"), which can be easily imported into
Cytoscape to visualize the network. The attribute file allows customizing
the metabolites of interest based on a score reflecting the degree of association
with the phenotype under study (i.e. log10(pvalue) adjusted for the sign of the
association).
}
\examples{
\donttest{## Build metabolic network: might take few minutes
metabolic_network = MWAS_KEGG_network(kegg_paths = KEGG_metabolic_paths[[1]][, 1])

## Test for association between diabetes and target_metabolites
T2D_model <- MWAS_stats (targetMetabo_SE, disease_id = "T2D",
                         confounder_ids = c("Age", "Gender", "BMI"),
                         assoc_method = "logistic")

## Select the metabolites of interest and get their corresponding KEGG IDs
T2D_model_subset = T2D_model[1:5, ]
kegg_metabolites = c("cpd:C00186", "cpd:C01089", "cpd:C00123", "cpd:C00183",
                     "cpd:C00407")

## Build shortest-path subnetwork
keggSP_subnetwork = MWAS_KEGG_shortestpaths(network_table = metabolic_network,
                                            metabolites = kegg_metabolites,
                                            MWAS_matrix = T2D_model_subset)
}
}
