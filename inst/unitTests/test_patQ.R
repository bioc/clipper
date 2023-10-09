library(gRbase)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"))
graph <- igraph::igraph.to.graphNEL(graph)

set.seed(1234)

exp1 <- matrix(0,5,20)
for (i in 1:5){
    exp1[i,] <- rnorm(20,i+5,1)
}
row.names(exp1)<-graphite::nodes(graph)

exp2 <-matrix(0,5,20)
for (i in 1:5){
    exp2[i,] <- rnorm(20,10-i,2)
}
row.names(exp2)<-graphite::nodes(graph)

colnames(exp1)<-LETTERS[1:20]
colnames(exp2)<-letters[1:20]
classes<-c(rep(2,20), rep(1,20))

exp <- cbind(exp1, exp2)

test_pathQ <- function(){
  set.seed(1234)
  pq <- pathQ(exp, classes, graph, nperm=100)
  checkEqualsNumeric(2, length(pq))
  checkEqualsNumeric(0, pq$alphaVar)
  checkEqualsNumeric(0, pq$alphaMean)
  row.names(exp) <- letters[1:5]
  checkException(pathQ(exp, classes, graph, nperm=100))
}
