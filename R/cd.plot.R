##### Plot to compare count distributions for two experimental conditions


### Generating data

cd.dat <- function (input, columns = c(1:2)) {
  
  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")
    
  
    
  if (length(columns) < 2)
    stop("ERROR! You must indicate at least to samples to compare.\n")
  
  if (!is.null(assayData(input)$exprs)) {
    if (ncol( assayData(input)$exprs) < 2)
      stop("ERROR! Your input object should have at least two different samples.\n")
    
    cond1 <- assayData(input)$exprs[,columns[1]]
    cond2 <- assayData(input)$exprs[,columns[2]]
  } else {
    if (ncol( assayData(input)$counts) < 2)
      stop("ERROR! Your input object should have at least two different samples.\n")
    
    cond1 <- assayData(input)$counts[,columns[1]]
    cond2 <- assayData(input)$counts[,columns[2]]
  }
  
  
  suma <- cond1 + cond2
  
  detect <- which(suma > 0)
  
  cond1.0 <- cond1[detect]
  cond2.0 <- cond2[detect]
  
  nu <- length(cond1.0) # number of detected features
  
  qq <- 1:min(nu,100)
  
  cum1 <- cumsum(sort(cond1.0, decreasing = TRUE))/sum(cond1.0)
  cum2 <- cumsum(sort(cond2.0, decreasing = TRUE))/sum(cond2.0)
  
  yy1 <- cum1[round(nu*qq/min(nu,100), 0)]*100
  yy2 <- cum2[round(nu*qq/min(nu,100), 0)]*100
  
  data2plot = data.frame("%detected_features" = qq,
                         "%cumulative_reads_1" = yy1,
                         "%cumulative_reads_2" = yy2)
  
  if (!is.null(assayData(input)$exprs)) {
    
    colnames(data2plot)[2:3] = colnames(assayData(input)$exprs)[columns]
        
  } else {
    
    colnames(data2plot)[2:3] = colnames(assayData(input)$counts)[columns]
    
  }
    
  data2plot
   
}


###########################################################################
###########################################################################
###########################################################################




### Generating plot

cd.plot <- function (data2plot,...) {
  
  plot(data2plot[,1], data2plot[,2], xlab = "% detected features", ylab = "% cumulative reads",
       type = "l", col = 2, main = "Count cumulative distribution")

  lines(data2plot[,1], data2plot[,3], col = 4)

  legend("bottom", legend = colnames(data2plot[-1]), text.col = c(2,4), bty = "n",
         lty = 1, lwd = 2, col = c(2,4), cex = 1.5)

}
