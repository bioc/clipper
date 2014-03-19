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

cliquePairedTest <- function(expr, classes, graph, nperm, alphaV=0.05, b=100, root=NULL, permute=TRUE, alwaysShrink=FALSE) {
  expr    <- getExpression(expr, classes)
  
  if (sum(classes==1) != sum(classes==2)) {
    stop("Your are working woth paired mode. The number of samples per class must be equal (and paired).")
  }
  
  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]
  
  cvt     <- runVarianceTest(expr, classes, graph, nperm, root, permute, alwaysShrink)

  if (is.null(cvt)){
    return(NULL)
  }
  
  check   <- cvt$alpha <= alphaV
  cliques <- cvt$cliques
  
  maxcliques <- max(sapply(cliques, length))
  
  ncl1 <- sum(classes==2)
  ncl2 <- sum(classes==1)
  
  alpha <- sapply(seq_along(cliques), function(i) {

    genes <- unlist(cliques[i])
    e  <- expr[, genes, drop=FALSE]
    e1 <- e[classes==2,, drop=FALSE]
    e2 <- e[classes==1,, drop=FALSE]
    
    if (length(genes) == 1)
      t.test(e1, e2, paired=TRUE)$p.value 

    else {
      e.diff <- e1 - e2
      e.bar <- colMeans(e.diff)
      e.centr <- e.diff - e.bar
      e.num <- nrow(e1)
      
      e.s <- (t(e.centr) %*% e.centr) / e.num
      t2 <- tryCatch(e.num * (t(e.bar) %*% solve(e.s) %*% e.bar), error=function(e) return(NA))
      
      p <- ncol(e1)
      np <- e.num - p
      t.value <- t2 * np / (p * (e.num-1))
      
      1-pf(t.value, p, np)
    }
  })
  
  list(alpha=alpha, cliques=cliques)
}
