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

runPathwayVar <- function(expr, classes, graph, nperm, permute) {
  e1 <- expr[classes==2,, drop=FALSE]
  e2 <- expr[classes==1,, drop=FALSE]
  
  adj <- as(moralize.graphNEL(graph), "matrix")

  cliques <- maxClique(moralize.graphNEL(graph))$maxCliques
  cliques <- lapply(cliques, function(x) match(x, nodes(graph)))
  
  maxcliques <- max(sapply(cliques, length))

  shrink <- sum(classes==1) < maxcliques || sum(classes==2) < maxcliques

  cov <- estimateCov(e1, e2, shrink)

  s1.hat <- qpIPF(cov$s1, cliques)
  s2.hat <- qpIPF(cov$s2, cliques)
  s.hat  <- qpIPF(cov$s,  cliques)
    
  k1.hat <- solve(s1.hat)
  k2.hat <- solve(s2.hat)
  k.hat  <- solve(s.hat)
  
  k1.det <- det(k1.hat)
  k2.det <- det(k2.hat)
  k.det  <- det(k.hat)

  lambda <- NROW(e1)*log(k1.det/k.det) + NROW(e2)*log(k2.det/k.det)

  arcs <- (sum(adj)/2) + NCOL(e1)
  lambdaT <- qchisq(0.95, arcs)
  
  alpha <- if(permute || shrink) {

    if (is.na(lambda))
      stop("Impossible to compute lambda. Maybe one determinant is equal to 0 or Inf.")
    
    perms <- getperms(length(classes), nperm)$perms
    counts <- apply(perms, 1, function(x)
                    {
                      classesP <- classes[x]
                      ep1 <- expr[classesP==2,, drop=FALSE]
                      ep2 <- expr[classesP==1,, drop=FALSE]
                      
                      covP <- estimateCov(ep1,ep2,shrink)
                      
                      s1.hat <- qpIPF(covP$s1, cliques)
                      s2.hat <- qpIPF(covP$s2, cliques)
                      s.hat  <- qpIPF(covP$s,  cliques)

                      k1.hat <- solve(s1.hat)
                      k2.hat <- solve(s2.hat)
                      k.hat  <- solve(s.hat)
                      
                      k1.det <- det(k1.hat)
                      k2.det <- det(k2.hat)
                      k.det  <- det(k.hat)
                      
                      NROW(e1)*log(k1.det/k.det) + NROW(e2)*log(k2.det/k.det) 
                    })

    if (all(is.na(counts)))
      stop("All permutational matrices have no calculable lambda")

    counts[is.na(counts)] <- lambda + 1
    
    sum(counts >= lambda)/nperm
  } else {
    1 - pchisq(lambda, arcs)
  }
  list(cov=cov, lambda=lambda, lambdaTeo=lambdaT, alpha=alpha, cliques=cliques)
}

pathQ <- function(expr, classes, graph, nperm=100, alphaV=0.05, b=100, permute=TRUE){
  expr <- getExpression(expr, classes)
  
  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))

  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")

  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]
  
  varTest <- runPathwayVar(expr, classes, graph, nperm, permute)
  check <- varTest$alpha <= alphaV

  if (is.na(check)){
    warning("Test on the concentration matrix is not calculable.")
    return(NA)
  }
  exp1 <- expr[classes==2,, drop=FALSE]
  exp2 <- expr[classes==1,, drop=FALSE]

  cli.moral <- varTest$cliques
  
  stat.obs  <- hoteIPF(exp1, exp2, check, cli.moral)
  stat.perm <- vector("numeric", nperm)
  
  for (i in seq_len(nperm)) {

    ind          <- sample(NROW(expr))
    exp1.perm    <- expr[ind[1:NROW(exp1)],, drop=FALSE]
    exp2.perm    <- expr[ind[(NROW(exp1)+1):NROW(expr)],, drop=FALSE]
    stat.perm[i] <- hoteIPF(exp1.perm, exp2.perm, check, cli.moral)
  }

  alpha <- sum(stat.perm >= stat.obs) / nperm
  
  list(alphaVar=varTest$alpha, alphaMean=alpha)
}
