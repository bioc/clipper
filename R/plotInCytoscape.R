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

plotInCytoscape <- function(graph, path, color="#6699FF", main="graph", layout="force-directed"){
    if(requireNamespace("RCy3", quitely=TRUE)) {
           path <- as.character(path)
           genesInvolved <- retrieveInvolvedGenes(path)
           
           g <- markMultiple(graph)
           g <- RCy3::initEdgeAttribute(g, "edgeType", "char", "undefined")
           g <- RCy3::initEdgeAttribute(g, "weight", "numeric", 1)
           
           cy <- RCy3::CytoscapeConnection()

           possible.layout.names <- RCy3::getLayoutNames(cy)
           if (!layout %in% possible.layout.names) {
               warning(paste0("Layout \"",layout,"\" not found. Grid layout setted. Please chose one of the following. \n", paste(possible.layout.names, collapse=", ")))
               layout="grid"
           }
           
           if (main %in% as.character(RCy3::getWindowList(cy)))
               RCy3::deleteWindow(cy, main)
           
           w <- RCy3::CytoscapeWindow(main, g)
           RCy3::deleteNodeAttribute(w,"node.fillColor")
           RCy3::displayGraph(w)
  
           RCy3::setNodeColorDirect(w, genesInvolved, color)
           RCy3::layoutNetwork(w, layout.name=layout)
           RCy3::redraw(w)
       } else {
           stop("Package RCy3 not installed. Please install.")
       }
}

markMultiple <- function(g) {
  d <- edgeData(g)
  if (length(d) == 0)
    return(g)

  ns <- names(d)
  for (i in 1:length(d)) {
    tp <- d[[i]]$edgeType
    if (length(grep(";", tp, fixed=T)) > 0) {
      nodes <- unlist(strsplit(ns[[i]], "|", fixed=T))
      edgeData(g, nodes[1], nodes[2], "edgeType") <- "multiple"
    }
  }

  return(g)
}


retrieveInvolvedGenes <- function(path, sep=";"){
  if (length(path) != 12)
    stop("Some path information are missing. Path must be of length 12.")
  unlist(strsplit(path[11], sep))
}

