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

mmmoralize <- function(graph) {
    g <- igraph::igraph.from.graphNEL(graph)
    m <- igraph::as_adjacency_matrix(g, sparse=F)
    m <- gRbase::moralizeMAT(m)
    g <- igraph::graph_from_adjacency_matrix(m, mode="directed")
    g
}
