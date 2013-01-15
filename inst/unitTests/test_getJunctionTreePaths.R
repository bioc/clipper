# Copyright 2012 Paolo Martini <paolo.cavei@gmail.com>

library(gRbase)
graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"))

test_getJunctionTreePaths <- function(){
  checkEquals(list(c(1,2)),getJunctionTreePaths(graph))
}
