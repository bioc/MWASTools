\docType{data}
\name{KEGG_metabolic_paths}
\alias{KEGG_metabolic_paths}
\title{KEGG human metabolic pathways}
\format{List}
\description{
  The first element of this list contains the KEGG identifiers (IDs) and names
  of 51 human metabolic pathways. By default, the function "MWAS_KEGG_network()"
  builds a reaction network using these KEGG IDs. The second element of the list
  is a matrix containing KEGG reactions with incorrect/inconsistent directionality.
  The directionality of these reactions has been corrected based on published
  literature. This matrix can be updated or edited by the user if required.
}
\usage{data(KEGG_metabolic_paths)}
\value{List}
