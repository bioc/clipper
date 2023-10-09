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

plotInCytoscape <- function(graph, path, color="#6699FF", main="graph"){
    if(requireNamespace("RCy3", quitely=TRUE)) {
        try(RCy3::deleteNetwork(main), silent = TRUE)
        path <- as.character(path)
        genesInvolved <- retrieveInvolvedGenes(path)
        
        g <- markMultiple(graph)
        suid <- RCy3::createNetworkFromGraph(g, main)
        
        genes <- graphite::nodes(g)
        nodesAttrib <- data.frame(id = genes, label = sub("^[^:]*:", "", genes),
                                  type = sub(":.*", "", genes), stringsAsFactors = FALSE)
        RCy3::loadTableData(nodesAttrib, "id", "node")
        RCy3::setNodeLabelMapping("label")
        RCy3::setNodeShapeMapping("type", nodesAttrib$type, nodesAttrib$shape)
        invisible(list(graph = g, suid = suid))
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

