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

cliqueParamTest <- function(cov, geneIx, n1, n2, alpha=FALSE) {
  s1D <- det(cov$s1[geneIx, geneIx, drop=FALSE])
  s2D <- det(cov$s2[geneIx, geneIx, drop=FALSE])
  sD  <- det(cov$s[geneIx, geneIx, drop=FALSE])

  if (all(s1D > 0,s2D > 0,sD > 0))
    l <- n1*log(sD/s1D) + n2*log(sD/s2D)
  else
    l <- NA
  
  if (alpha) {
    g <- length(geneIx)
    a <- 1 - pchisq(l, g*(g+1)/2)
    c(l, a)
  } else
    l
}
