\name{nameCliques}
\alias{nameCliques}
\title{Generate clique names from their own elements.}

\description{Starting from the sorted elements of each clique of the
  list, this function generates names fusing in a string the element names.}

\usage{nameCliques(cliques)}
\arguments{
  \item{cliques}{a list where each element is a clique.}
}
\value{vector of strings}

\examples{
toyCliques <- list(c(45,36,90), c(36,1000,35))
nameCliques(toyCliques)
}