\name{clipperAllRoots}
\alias{clipperAllRoots}
\title{Dissect the pathway to find the path with the greatest
  association with phenotype.}
\description{Basing on either variance or mean clique test, this
  function identifies the paths that are mostly related with the
  phenotype under study.}
\usage{clipperAllRoots(expr, classes, graph, method=c("variance","mean",
"both", "paired"), nperm=100, alphaV=0.05, b=100, trZero=0.001, signThr=0.05,
maxGap=1, permute=TRUE, alwaysShrink=FALSE)}
\arguments{
  \item{expr}{an expression  matrix or ExpressionSet with colnames for
    samples and row name for genes.}
  \item{classes}{vector of 1,2 indicating the classes of samples (columns).}
  \item{graph}{a \code{graphNEL} object.}
  \item{method}{the kind of test to perform on the cliques. It could be
    mean, variance, mixed (the best between variance and mean) or paired
    mean.}
  \item{nperm}{number of permutations. Default = 100.}
  \item{alphaV}{pvalue threshold for variance test to be used during mean test. Default = 0.05.}
  \item{b}{number of permutations for mean analysis. Default = 100.}
  \item{trZero}{lowest pvalue detectable. This threshold avoids that
    -log(p) goes infinite.}
  \item{signThr}{significance threshold for clique pvalues.}
  \item{maxGap}{allow up to maxGap gaps in the best path
    computation. Default = 1.}
  \item{permute}{always performs permutations in the concentration
    matrix test. If FALSE, the test is made using the asymptotic
    distribution of the log-likelihood ratio. This option should be use
    only if samples size is >=40 per class.}
  \item{alwaysShrink}{always perform the shrinkage estimates of
    variance.}
}

\details{
  The both method combines the results obtained from the mean and
  variance test. In particular it assign to the cliques the minimum of
  mean and variance p-values.
}

\value{A matrix with a row for each paths. Rownames have the form:

  roots-paths.

  Columns are organized as follows:
  \enumerate{
    \item Index of the starting clique
    \item Index of the ending clique
    \item Index of the clique where the maximum value is reached
    \item Length of the path
    \item Maximum score of the path
    \item Average score along the path
    \item Percentage of path activation
    \item Impact of the path on the entire pathway
    \item Cliques involved and significant
    \item Cliques forming the path
    \item Genes forming the significant cliques
    \item Genes forming the path
  }
}

\references{
  Martini P, Sales G, Massa MS, Chiogna M, Romualdi C. Along signal
  paths: an empirical gene set approach exploiting pathway
  topology. NAR. 2012 Sep.

  Massa MS, Chiogna M, Romualdi C. Gene set analysis exploiting the
  topology of a pathway. BMC System Biol. 2010 Sep 1;4:121.
}

\seealso{\code{\link{cliqueVarianceTest}}, \code{\link{cliqueMeanTest}},
\code{\link{getJunctionTreePaths}}}

\examples{
if (require(graphite) & require(ALL)){
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  genes <- nodes(graph)
  data(ALL)
  all <- ALL[1:length(genes),1:20]
  classes <- c(rep(1,10), rep(2,10))
  featureNames(all@assayData)<- genes
  graph <- subGraph(genes, graph)
  clipped <- clipperAllRoots(all, classes, graph, "var", trZero=0.01, permute=FALSE)
  clipped[,1:5]
}
}
