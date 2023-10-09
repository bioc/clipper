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


test_cliqueVarianceTest <- function(){
  set.seed(1234)
  t <- cliqueVarianceTest(exp, classes, graph, 100)
  checkEqualsNumeric(t$alpha, c(0,0))
  checkEqualsNumeric(length(t),2)
  checkEqualsNumeric(length(t$cliques),2)
  row.names(exp) <- letters[1:5]
  checkException(cliqueVarianceTest(exp, classes, graph, 100))
}
