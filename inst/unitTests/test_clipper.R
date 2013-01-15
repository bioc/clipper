library(gRbase)

graph <- dag(c("me","ve"),c("me","al"),c("ve","al"),c("al","an"),c("al","st"),c("an","st"),c("a","b","c"), c("c","d"))

set.seed(1234)

exp1 <- matrix(0,9,20)
for (i in 1:5){
    exp1[i,] <- rnorm(20,i+5,1)
}
for (i in 1:4){
    exp1[5+i,] <- rnorm(20,i+5,1)
}

row.names(exp1)<-nodes(graph)

exp2 <-matrix(0,9,20)
for (i in 1:5){
    exp2[i,] <- rnorm(20,10-i,2)
}
for (i in 1:4){
    exp2[5+i,] <- rnorm(20,i+5,1)
}

row.names(exp2)<-nodes(graph)

colnames(exp1)<-LETTERS[1:20]
colnames(exp2)<-letters[1:20]
classes<-c(rep(2,20), rep(1,20))

exp <- cbind(exp1, exp2)

test_clipper <- function(){
  set.seed(1432)
  clipped <- clipper(exp, classes, method="var",graph, 100)
  checkEqualsNumeric(-log(0.001)*2, as.numeric(clipped[1,5]))
  checkEqualsNumeric(c(1,12), dim(clipped))
}
