
toySummary <- matrix(c(1,1,5,3,5,2,5,3,8.2,3,2,1,0.3,0.1,2,1,"1;2;3;4;5","1;2;3","1;2;3;4;5","1;2;3","1;2;3;4;5","1;2;3","1;2;3;4;5","1;2;3"),2,12)
row.names(toySummary) <- c("1;5","1;3")

test_prunePaths <- function(){
  pp <- prunePaths(toySummary, thr=0.1)
  checkEqualsNumeric(c(1,12), dim(pp))
  checkEquals("1;5", row.names(pp))
}
  
