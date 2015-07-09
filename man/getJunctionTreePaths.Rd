\name{getJunctionTreePaths}
\alias{getJunctionTreePaths}
\title{Extract the shortest paths along the junction tree of the graph.}

\description{Find the shortest paths in the Junction tree designed with the cliques of the graph.}

\usage{getJunctionTreePaths(graph, root=NULL)}
\arguments{
  \item{graph}{a \code{graphNEL} object.}
  \item{root}{nodes by which rip ordering is performed (as far as possible) on the variables using the maximum cardinality search algorithm.}
}
\value{list of clique indices representing the shortest paths of the graph.}

\references{
  Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal
  paths: an empirical gene set approach exploiting pathway
  topology. NAR. 2012 Sep.

  Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the
  topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
}

\examples{
if (require(graphite)) {
  kegg  <- pathways("hsapiens", "kegg")
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  getJunctionTreePaths(graph)
}
}
