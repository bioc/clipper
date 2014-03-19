# Copyright 2012-2013 Paolo Martini <paolo.martini@unipd.it>
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

hotePaired <- function(exp1, exp2, cliques, performPerm=FALSE, alwaysShrink=FALSE) {
  
  exp1.num <- nrow(exp1)
  exp.diff <- exp1 - exp2

  if (performPerm) {
    signs <- matrix(sample(c(1,-1), exp1.num * NCOL(exp1), replace=TRUE),
                    nrow=exp1.num)
    exp.diff <- exp.diff * signs
  }
  
  exp.bar <- colMeans(exp.diff)
  
  maxcliques <- max(sapply(cliques, length))
  shrink <- exp1.num < maxcliques || alwaysShrink
  exp.s <- estimateExprCov(exp.diff, shrink)
  
  exp.s <- qpIPF(exp.s, cliques)
  
  t2 <- exp1.num * (t(exp.bar) %*% solve(exp.s) %*% exp.bar)
  if (length(t2) == 1) {
    t2 <- as.numeric(t2)
  } else {
    warning("t2 longer than one.")
  }
    
  p <- ncol(exp1)
  np <- exp1.num - p

  t2 * np / (p * (exp1.num-1))
}
