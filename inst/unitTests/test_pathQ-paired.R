library(gRbase)
library(MASS)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"))

set.seed(1234)

sigma <- matrix(0,20,20)
diag(sigma)<-1

exp1 <- mvrnorm(5, rep(5,20), sigma)
row.names(exp1)<-nodes(graph)

exp2 <- mvrnorm(5, rep(5,20), sigma)
row.names(exp2)<-nodes(graph)

colnames(exp1)<-LETTERS[1:20]
colnames(exp2)<-letters[1:20]
classes<-c(rep(2,20), rep(1,20))

exp <- cbind(exp1, exp2)

test_pathQ <- function(){
  set.seed(1234)
  pq <- pathQ(exp, classes, graph, nperm=100, paired=TRUE)
  checkEqualsNumeric(2, length(pq))
  checkTrue(is.null(pq$alphaVar))
  checkTrue(pq$alphaMean > 0.05)
  checkException(cliquePairedTest(exp, classes[-1], graph, 100))
  row.names(exp) <- letters[1:5]
  checkException(pathQ(exp, classes, graph, nperm=100))
}

