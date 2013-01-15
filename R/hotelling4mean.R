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

hote4mean <- function(exp1, exp2, exact) {
  exp1.num <- nrow(exp1)
  exp2.num <- nrow(exp2)
  gene.num <- ncol(exp1)
  
  exp1.bar <- colMeans(exp1)
  exp2.bar <- colMeans(exp2)
  exp1.s   <- cov(exp1)
  exp2.s   <- cov(exp2)
  
  exp.diff <- exp1.bar - exp2.bar
  
  if (exact) {
    s  <- ((exp1.num-1)*exp1.s + (exp2.num-1)*exp2.s) / (exp1.num + exp2.num - 2)
    
    t2 <- tryCatch(((exp1.num*exp2.num) / (exp1.num+exp2.num)) * (exp.diff %*% solve(s) %*% exp.diff), error=function(e) return(NA))
    c <- exp1.num + exp2.num - gene.num - 1
    t.obs <- t2 * c / (gene.num * (exp1.num + exp2.num - 2))
    alpha.obs <- 1 - pf(t.obs, gene.num, c)
    list(alpha.obs=alpha.obs, t.obs=t.obs)
  } else {
    s <- exp1.s/exp1.num + exp2.s/exp2.num
    tryCatch(as.numeric(exp.diff %*% solve(s) %*% exp.diff), error=function(e) return(NA))
  }
}
