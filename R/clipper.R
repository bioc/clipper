# Copyright 2012 Paolo Martini <paolo.martini@unipd.it>
#
#
# This file is part of clipper.
#
# clipper is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License
# version 3 as published by the Free Software Foundation.
#
# clipper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with clipper. If not, see <http://www.gnu.org/licenses/>.

runCoreClipper <- function(cliqueTest, pathList, trZero, thr, maxGap){  
  cliques <- cliqueTest$cliques
  clNames <- nameCliques(cliques)
  alphas <- cliqueTest$alpha
  names(alphas) <- clNames

  sapply(pathList, function(idx) {
    formatBestSubPath(alphas, idx, trZero, thr, maxGap)
  })
}

clipper <- function(expr, classes, graph, method=c("variance","mean"), nperm=100, alphaV=0.05, b=100, root=NULL, trZero=0.001, signThr=0.05, maxGap=1, permute=TRUE){
  expr <- t(getExpression(expr, classes))
  expGenes <- row.names(expr)
  genes <- nodes(graph)
  genes <- intersect(genes, expGenes)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  graph <- subGraph(genes, graph)
  method <- match.arg(method)
  fu <- switch(method,
               variance = cliqueVarianceTest,
               mean     = cliqueMeanTest)
  
  ct <- fu(expr, classes, graph, nperm, alphaV, b, root, permute)
  jtp <- getJunctionTreePaths(graph, root)
  
  if (is.null(jtp))
    return(NULL)
  
  clipped <- runCoreClipper(ct, jtp, trZero, signThr, maxGap)
  clpprNames <- c("startIdx", "endIdx", "maxIdx", "lenght", "maxScore", "aScore", "activation", "impact", "involvedCliques", "cliqueOnPath", "involvedGenes", "pathGenes")
  
  if (!is.matrix(clipped))
    clipped <- as.matrix(clipped)
  
  clipped <- t(clipped)
  colnames(clipped) <- clpprNames
  clipped <- addNames(clipped)
  clipped <- removeNArows(clipped)
  as.data.frame(clipped, stringsAsFactors=FALSE)
}
