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

clipper <- function(expr, classes, graph, method=c("variance","mean","both", "paired"), nperm=100, alphaV=0.05, b=100, root=NULL, trZero=0.001, signThr=0.05, maxGap=1, permute=TRUE){
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
               mean     = cliqueMeanTest,
               both     = cliqueMixedTest,
               paired   = cliquePairedTest)
  
  ct <- fu(expr, classes, graph, nperm, alphaV, b, root, permute)

  if (is.null(ct)){
    return(NULL)
  }
  
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

clipperAllRoots <- function(expr, classes, graph, method=c("variance","mean","both", "paired"), nperm=100, alphaV=0.05, b=100, trZero=0.001, signThr=0.05, maxGap=1, permute=TRUE){
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
               mean     = cliqueMeanTest,
               both     = cliqueMixedTest,
               paired   = cliquePairedTest)

  roots <- unique(unlist(sapply(getGraphEntryGenes(graph, TRUE), function(x) {
    if (length(x) == 0)
      return(NULL)
    x[1]
  })))

  
  allTests <- lapply(roots, function(root) {
    
    ct <- fu(expr, classes, graph, nperm, alphaV, b, root, permute)
    if (is.null(ct)){
      return(NULL)
    }
    
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
    
  })
  
  rnull <- length(allTests) + 1

  ct <- fu(expr, classes, graph, nperm, alphaV, b, root=NULL, permute)
  if (is.null(ct)){
    return(NULL)
  }
  
  jtp <- getJunctionTreePaths(graph, root=NULL)
  
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
  allTests[[rnull]] <- as.data.frame(clipped, stringsAsFactors=FALSE)

  names(allTests) <- c(roots,"null")
  
  maxScore <- 0
  maxRoot <- NULL
  for (r in names(allTests)) {
    if (NROW(allTests[[r]])!=0){
      maxRelativeScore <- max(as.numeric(allTests[[r]][,5]))
      if (maxScore < maxRelativeScore) {
        maxScore <- maxRelativeScore
        maxRoot <- r
      }
    }
  }
  if (is.null(maxRoot))
    return(NULL)
  
  row.names(allTests[[maxRoot]]) <- paste(maxRoot, row.names(allTests[[maxRoot]]), sep="-")
  allTests[[maxRoot]]
  
}

