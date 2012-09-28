ranking <-
function(results) {

  M <- results$M

  D <- results$D

  prob <- results$prob

  ## Changing NA by 0
  M[is.na(M)] <- 0
  D[is.na(D)] <- 0
  prob[is.na(prob)] <- 0


  ## Ranking

#   ranking1 <- M*prob
# 
#   ranking2 <- sign(M)*prob
# 
#   ranking3 <- M*D
# 
#   ranking4 <- M*D*prob

  ranking5 <- sqrt(M*M + D*D)*sign(M)


  ## Ranking results
  #list(ranking1, ranking2, ranking3, ranking4, ranking5)
  theranking <- data.frame(rownames(results), ranking5)
  rownames(theranking) <- NULL
  colnames(theranking) <- c("ID", "statistic")
  
  theranking

}

