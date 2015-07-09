\name{cliquePairedTest}
\alias{cliquePairedTest}
\title{Paired mean test for cliques.}
\description{It decomposes the graph in cliques and performs the paired mean test in every one.}
\usage{cliquePairedTest(expr, classes, graph, nperm, alphaV=0.05, b=100,
root=NULL, permute=TRUE, alwaysShrink=FALSE)}
\arguments{
  \item{expr}{an expression  matrix or ExpressionSet with colnames for
    samples and row name for genes.}
  \item{classes}{vector of 1,2 indicating the classes of samples
    (columns). It is assumed that class labels are ordered so that the first occurrence of class 2 is
    paired with the first occurrence of class 1 and so on.}
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
  \item{alwaysShrink}{always perform the shrinkage estimates of
    variance.}
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
  kegg  <- pathways("hsapiens", "kegg")
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  genes <- nodes(graph)
  data(ALL)
  all <- ALL[1:length(genes),1:20]
  classes <- c(rep(1,10), rep(2,10))
  featureNames(all@assayData)<- genes
  graph <- subGraph(genes, graph)
  cliquePairedTest(all, classes, graph, nperm=100, permute=FALSE)$alpha
}
}