\name{MWAS_KEGG_network}
\alias{MWAS_KEGG_network}
\title{Build a KEGG-based metabolic network}
\description{
This function generates a KEGG-based metabolic network connecting substrate-product
pairs. The network is formated as a four-column matrix, where each row represents
an edge connecting two metabolites (from metabolite in column 1 to metabolite in
column 2). The third column contains the identifiers (IDs) of the reactions
performing each of the metabolic conversions, while the fourth column indicates the
direction of each reaction.
}
\usage{
MWAS_KEGG_network(kegg_paths = NULL)
}
\arguments{
  \item{kegg_paths}{character vector containing the KEGG IDs of the metabolic
  pathways of interest (organism-specific). For example, the KEGG ID for the
  human "glycolysis/gluconeogenesis" pathway is"hsa00010". By default, the KEGG
  IDs contained in the dataset "KEGG_metabolic_paths" will be used.}
}
\note{
Like in the MetaboSignal package, reaction directionality has been cross-checked
and corrected (when required) based on previous literature (Duarte et al., 2007).
}
\references{
Duarte NC, et al. (2007). Global reconstruction of the human metabolic network
based on genomic and bibliomic data. Proceedings of the National Academy of Sciences,
104, 1777-1782.

Posma JM, et al. (2014). MetaboNetworks, an interactive Matlab-based toolbox
for creating, customizing and exploring sub-networks from KEGG. Bioinformatics,
30, 893-895.

Rodriguez-Martinez A, et al. (2017).MetaboSignal: a network-based approach for
topological analysis of metabotype regulation via metabolic and signaling pathways.
Bioinformatics, 33, 773-775.

Zhang JD, Wiemann S. (2009). KEGGgraph: a graph approach to KEGG PATHWAY in
R and Bioconductor. Bioinformatics, 25, 1470-1471.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\value{
A four-column matrix where each row represents an edge between two nodes.
}
\examples{
\donttest{## Define the paths that will be used to build the network
data(KEGG_metabolic_paths)
metabo_paths =  KEGG_metabolic_paths[[1]][, 1]

## Build metabolic network-table: might take few minutes
metabolic_network = MWAS_KEGG_network(kegg_paths = metabo_paths)
}
}
