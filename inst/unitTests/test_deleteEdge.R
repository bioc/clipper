# Copyright 2012 Paolo Martini <paolo.cavei@gmail.com>

library(gRbase)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"))
graph.out <- dag(c("me"), c("ve"), c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"))
graph <- igraph::igraph.to.graphNEL(graph)
graph.out <- igraph::igraph.to.graphNEL(graph.out)

checkSameGraph <- function(g1, g2) {
  identical(as(g1, "matrix"),as(g1, "matrix"))
}

test_deleteEdge <- function(){
  checkTrue(checkSameGraph(deleteEdge(graph, "ve", "me"), graph.out))
  checkException(deleteEdge(graph, "me", "ve"))
}
