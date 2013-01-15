\name{pathwayTest.R}
\alias{pathQ}
\title{Whole pathway test using qpipf.}
\description{Performs variance and mean test using qpipf on the whole pathway.}
\usage{pathQ(expr, classes, graph, nperm=100, alphaV=0.05, b=100, permute=TRUE)}
\arguments{
  \item{expr}{an expression  matrix or ExpressionSet with colnames for samples and
    rownames for expression features.}
  \item{classes}{vector of 1,2 indicating the classes of the samples (columns).}
  \item{graph}{a \code{graphNEL} object.}
  \item{nperm}{number of permutations. Default = 100.}
  \item{alphaV}{pvalue significance threshold for variance test to be used during
    mean test. Default = 0.05.}
  \item{b}{number of permutations for mean analysis. Default = 100.}
  \item{permute}{always performs permutations in the concentration
    matrix test. If FALSE, the test is made using the asymptotic
    distribution of the log-likelihood ratio. This option should be use
    only if samples size is >=40 per class.
  }
}
\value{a list with alphaVar (pvalue for the variance test) and alphaMean (pvalue for mean test).}

\references{
  Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal
  paths: an empirical gene set approach exploiting pathway
  topology. NAR. 2012 Sep.

  Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the
  topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
}

\examples{
if (require(graphite) & require(ALL)){
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  genes <- nodes(graph)
  data(ALL)
  all <- ALL[1:length(genes),1:20]
  classes <- c(rep(1,10), rep(2,10))
  featureNames(all@assayData)<- genes
  graph <- subGraph(genes, graph)
  pathQ(all, classes, graph, nperm=100, permute=FALSE)
}
}

\note{
  This function is based on the Gaussian Graphical Models and to use it
  in a proper way it is necessary that the graph is an Direct Acyclic
  Graph. Please check any graph in input using \code{isAcyclic} from ggm package.
}
  
