\name{deleteEdge}
\alias{deleteEdge}
\title{Remove an edge from \code{graphNel} object.}

\description{Remove from a \code{graphNEL} object the edge specified.}

\usage{deleteEdge(graph, from, to)}
\arguments{
  \item{graph}{a graphNEL object.}
  \item{from}{a string with the name of the node where the edge start.}
  \item{to}{a string with the name of the node where the edge end.}
}
\value{a \code{graphNEL} object.}

\examples{
if (require(graphite)) {
  kegg  <- pathways("hsapiens", "kegg")
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  head(edges(graph))
  ## We are going to remove the edge 1026-1019
  head(edges(deleteEdge(graph, "1026", "1019")))
}
}