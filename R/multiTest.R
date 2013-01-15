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

multT <- function(exp1, exp2, perm.num) {
  exp1.num <- nrow(exp1)
  exps     <- rbind(exp1, exp2)
  exps.num <- nrow(exps)
  
  t.obs <- hote4mean(exp1, exp2, FALSE)
  
  stat.perm <- vector("numeric", perm.num)
  for (i in 1:perm.num) {
    ind          <- sample(exps.num)
    exp1.perm    <- exps[ind[1:exp1.num],, drop=FALSE]
    exp2.perm    <- exps[ind[(exp1.num+1):exps.num],, drop=FALSE]
    stat.perm[i] <- hote4mean(exp1.perm, exp2.perm, FALSE)
  }
  
  alpha.obs <- sum(stat.perm >= t.obs) / perm.num
  list(alpha.obs=alpha.obs, t.obs=t.obs)
}
