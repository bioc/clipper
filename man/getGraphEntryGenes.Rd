\name{getGraphEntryGenes.R}
\alias{getGraphEntryGenes}
\title{Extract all the possible entry point (genes with no entering edges) from graph.}
\description{It extracts the possible entry point of the graph. Entry
  points are defined as nodes with no entering edges.}
\usage{getGraphEntryGenes(graph, byCliques=FALSE, root=NULL)}
\arguments{
  \item{graph}{a \code{graphNEL} object.}
  \item{byCliques}{when TRUE it returns a list where entry point are organized by cliques.}
  \item{root}{nodes by which rip ordering is performed (as far as possible) on the variables using the maximum cardinality search algorithm.}
}
\value{a vector of gene names representing the entry point of graph.}

\references{
  Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal
  paths: an empirical gene set approach exploiting pathway
  topology. NAR. 2012 Sep.

  Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the
  topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
}

\examples{
if (require(graphite)) {
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  getGraphEntryGenes(graph)
}
}
