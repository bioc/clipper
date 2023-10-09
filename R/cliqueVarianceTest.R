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

runVarianceTest <- function(expr, classes, graph, nperm, root, permute, alwaysShrink) {

  gns  <- colnames(expr)
  cliques <- extractCliquesFromDag(graph, root=root)

  if (is.null(cliques)){
    warning("No cliques available or the DAG provided can not be ripped. Please check if your input graph is a DAG.")
    return(NULL)
  }
  
  maxcliques <- max(sapply(cliques, length))
  
  ncl1 <- sum(classes==2)
  ncl2 <- sum(classes==1)

  shrink <- ncl1 < maxcliques || ncl2 < maxcliques || alwaysShrink

  cov  <- estimateCov(expr[classes==2,, drop=FALSE], expr[classes==1,, drop=FALSE], shrink)
  
  geneIxL <- lapply(cliques, function(c) match(c, gns))
  
  cpt <- sapply(geneIxL,
                function(ix) {
                  cliqueParamTest(cov, ix, ncl1, ncl2, TRUE)
                })
  
  alpha <- if (permute || shrink) {
    perms <- getperms(length(classes), nperm)$perms
    counts <- apply(perms, 1, function(x)
                    {
                      classesP <- classes[x]
                      ep1 <- expr[classesP==2,, drop=FALSE]
                      ep2 <- expr[classesP==1,, drop=FALSE]
                      covP <- estimateCov(ep1,ep2, shrink)
                      sapply(geneIxL,
                             function(ix) {
                               cliqueParamTest(covP, ix, ncl1, ncl2)
                             })
                    })
    
    if (!is.matrix(counts))
      counts <- t(as.matrix(counts))
    
    apply(counts >= cpt[1,], 1, sum)/nperm

  } else
    cpt[2,]
  
  list(cov=cov, alpha=alpha, cliques=cliques)
}

cliqueVarianceTest <- function(expr, classes, graph, nperm, alphaV=0.05, b=100, root=NULL, permute=TRUE, alwaysShrink=FALSE) {
  expr <- getExpression(expr, classes)

  genes <- graphite::nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  graph <- KEGGgraph::subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]
  
  res <- runVarianceTest(expr, classes, graph, nperm, root, permute, alwaysShrink)
  res$cov <- NULL
  res
}
