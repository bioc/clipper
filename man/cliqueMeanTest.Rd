\name{cliqueMeanTest.R}
\alias{cliqueMeanTest}
\title{Mean test for cliques.}
\description{It decomposes the graph in cliques and performs the mean test in every one.}
\usage{cliqueMeanTest(expr, classes, graph, nperm, alphaV=0.05, b=100,
root=NULL, permute=TRUE)}
\arguments{
  \item{expr}{an expression  matrix or ExpressionSet with colnames for
    samples and row name for genes.}
  \item{classes}{vector of 1,2 indicating the classes of samples (columns).}
  \item{graph}{a \code{graphNEL} object.}
  \item{nperm}{number of permutations.}
  \item{alphaV}{pvalue threshold for variance test to be used during mean test.}
  \item{b}{number of permutations for mean analysis.}
  \item{root}{nodes by which rip ordering is performed (as far as
    possible) on the variables using the maximum cardinality search
    algorithm.}
  \item{permute}{always performs permutations in the concentration
    matrix test. If FALSE, the test is made using the asymptotic
    distribution of the log-likelihood ratio. This option should be use
    only if samples size is >=40 per class.}
}
\value{a list with alphas (vector of cliques pvalues based on the variance test) and cliques (list of the cliques and related elements).}

\references{
  Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal
  paths: an empirical gene set approach exploiting pathway
  topology. NAR. 2012 Sep.

  Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the
  topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
}

\seealso{\code{\link{cliqueVarianceTest}}.}

\examples{
if (require(graphite) & require(ALL)){
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  genes <- nodes(graph)
  data(ALL)
  all <- ALL[1:length(genes),1:20]
  classes <- c(rep(1,10), rep(2,10))
  featureNames(all@assayData)<- genes
  graph <- subGraph(genes, graph)
  cliqueMeanTest(all, classes, graph, nperm=100, permute=FALSE)$alpha
}
}