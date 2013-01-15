library(gRbase)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"))

test_getGraphEntryGenes <- function(){
  checkEquals(getGraphEntryGenes(graph), "st")
}
