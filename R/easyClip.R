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


easyClip <- function (expr, classes, graph, method=c("variance","mean"), pathThr=0.05, pruneLevel=0.2, nperm=100, alphaV=0.05, b=100, root=NULL, trZero=0.001, signThr=0.05, maxGap=1, permute=TRUE) {

  method <- match.arg(method)
  fu <- switch(method,
               variance = "variance",
               mean     = "mean")
  
  pTest <- pathQ(expr, classes, graph, nperm, alphaV, b, permute)
  
  if ((method == "variance") & (pTest$alphaVar > pathThr))
    return(NULL)
  
  if ((method == "mean") & (pTest$alphaMean > pathThr))
    return(NULL)
  
  clipped <- clipper(expr=expr, classes=classes, graph=graph, method=method, nperm=nperm, alphaV=alphaV, b=b, root=root, trZero=trZero, signThr=signThr, maxGap=maxGap, permute=permute)
  
  prunePaths(clipped, pruneLevel)
}
