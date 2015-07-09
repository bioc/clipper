\name{easyClip}
\alias{easyClip}

\title{Easy clip analysis.}

\description{Easy clip function allows the full exploitation of Clipper
  Package features in a unique and easy to use function. Starting from
  an expression matrix and a pathway, these function extact the most
  transcriptionally altered portions of the graph.}

\usage{easyClip(expr, classes, graph, method=c("variance","mean"),
pathThr=0.05, pruneLevel=0.2, nperm=100, alphaV=0.05, b=100,
root=NULL, trZero=0.001, signThr=0.05, maxGap=1, permute=TRUE)}

\arguments{
  \item{expr}{an expression  matrix or ExpressionSet with colnames for
    samples and row name for genes.}
  \item{classes}{vector of 1,2 indicating the classes of samples (columns).}
  \item{graph}{a \code{graphNEL} object.}
  \item{method}{the kind of test to perform on the cliques. It could be
    either mean or variance.}
  \item{pathThr}{The significance threshold of the whole pathway
    test. Deafault = 0.05}
  \item{pruneLevel}{a dissimilarity threshold. NULL means no pruning.}
  \item{nperm}{number of permutations. Default = 100.}
  \item{alphaV}{pvalue threshold for variance test to be used during
    mean test. Default = 0.05.}
  \item{b}{number of permutations for mean analysis. Default = 100.}
  \item{root}{nodes by which rip ordering is performed (as far as
    possible) on the variables using the maximum cardinality search
    algorithm.}
  \item{trZero}{lowest pvalue detectable. This threshold avoids that
    -log(p) goes infinite.}
  \item{signThr}{significance threshold for clique pvalues.}
  \item{maxGap}{allow up to maxGap gaps in the best path
    computation. Default = 1.}
  \item{permute}{always performs permutations in the concentration
    matrix test. If FALSE, the test is made using the asymptotic
    distribution of the log-likelihood ratio. This option should be use
    only if samples size is >=40 per class.}
}

\value{a matrix with row as the different paths. Columns are organized
  as follwes:
  1 - Index of the starting clique
  2 - Index of the ending clique
  3 - Index of the clique where the maximum value is reached
  4 - length of the path
  5 - maximum score of the path
  6 - average score along the path
  7 - percentage of path activation
  8 - impact of the path on the entire pathway
  9 - clique involved and significant
  10 - clique forming the path
  11 - genes forming the significant cliques
  12 - genes forming the path)}

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
  kegg  <- pathways("hsapiens", "kegg")
  graph <- pathwayGraph(convertIdentifiers(kegg$'Chronic myeloid leukemia', "entrez"))
  genes <- nodes(graph)
  data(ALL)
  all <- ALL[1:length(genes),1:20]
  classes <- c(rep(1,10), rep(2,10))
  featureNames(all@assayData)<- genes
  graph <- subGraph(genes, graph)
  clipped <- easyClip(all, classes, graph, nperm=10)
  clipped[,1:5]
}
}
