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

hoteIPF <- function(exp1, exp2, exact, cliques, alwaysShrink=FALSE) {
  exp1.num <- nrow(exp1)
  exp2.num <- nrow(exp2)
  gene.num <- ncol(exp1)

  exp1.bar <- colMeans(exp1)
  exp2.bar <- colMeans(exp2)
  
  maxcliques <- max(sapply(cliques, length))
  shrink <- exp1.num < maxcliques || exp2.num <maxcliques || alwaysShrink
  cov <- estimateCov(exp1, exp2, shrink)

  exp1.s   <- cov$s1
  exp2.s   <- cov$s2
  
  exp.diff <- exp1.bar - exp2.bar
  
  exp1.s <- qpIPF(exp1.s, cliques) ## Ma le passo io
  exp2.s <- qpIPF(exp2.s, cliques)
  
  if (exact) {
    s  <- ((exp1.num-1)*exp1.s + (exp2.num-1)*exp2.s) / (exp1.num + exp2.num - 2)
    t2 <- ((exp1.num*exp2.num) / (exp1.num+exp2.num)) * (exp.diff %*% solve(s) %*% exp.diff)    
    as.numeric(t2)
  } else {
    s <- exp1.s/exp1.num + exp2.s/exp2.num
    tryCatch(as.numeric(exp.diff %*% solve(s) %*% exp.diff), error=function(e) return(NA))
  }
}
