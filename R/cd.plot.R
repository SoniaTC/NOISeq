##### Plot to compare count distributions for two or more samples


### Generating data

cd.dat <- function (input) {
  
  if (inherits(input,"eSet") == FALSE)
    stop("ERROR: The input data must be an eSet object.\n")
    
  
  if (!is.null(assayData(input)$exprs)) {
    if (ncol( assayData(input)$exprs) < 2)
      stop("ERROR: The input object should have at least two samples.\n")
  
    datos <- assayData(input)$exprs
    
  } else {
    if (ncol( assayData(input)$counts) < 2)
      stop("ERROR: The input object should have at least two samples.\n")
    
    datos <- assayData(input)$counts
    
  }
  
    
  datos <- datos[which(rowSums(datos) > 0),]
  
  nu <- nrow(datos) # number of detected features
  
  qq <- 1:nu
    
  data2plot = data.frame("%features" = 100*qq/nu)
  
  for (i in 1:ncol(datos)) {
    
    acumu <- 100*cumsum(sort(datos[,i], decreasing = TRUE))/sum(datos[,i])
          
    data2plot = cbind(data2plot, acumu)    
  }
  
  
  colnames(data2plot)[-1] = colnames(datos)
  
  
  
  #### Diagnostic test
  
  KSpval = mostres = NULL
  
  for (i in 1:(ncol(datos)-1)) {
    for (j in (i+1):ncol(datos)) {      
      mostres = c(mostres, paste(colnames(datos)[c(i,j)], collapse = "_"))
      #KSpval = c(KSpval, ks.test(datos[,i], datos[,j], alternative = "two.sided")$"p.value")
      KSpval = c(KSpval, suppressWarnings(ks.test(datos[,i], datos[,j], alternative = "two.sided"))$"p.value")
    }
  }
  
  KSpval = p.adjust(KSpval, method = "fdr")
  print("Summary of FDR adjusted p-values:")
  print(summary(KSpval))
    
  if (sum(KSpval < 0.05) > 0) {
    print("Diagnostic test: FAILED. Normalization is required to correct this bias.")
    print("According to Kolmogorov-Smirnov tests, at least a pair of samples have significantly different distributions")
    print("Minimum adjusted p-value was: "); print(min(KSpval, na.rm = TRUE))
  }  else {
    print("Diagnostic test: PASSED.")
  }
  
  
  #### Results
  
  list("data2plot" = data2plot, 
       "DiagnosticTest" = data.frame("ComparedSamples" = mostres, "KSpvalue" = KSpval))
   
}


###########################################################################
###########################################################################
###########################################################################




### Generating plot

cd.plot <- function (dat, samples = NULL,...) {
  
  dat = dat$data2plot
  
  if (is.null(samples)) samples <- 1:(ncol(dat)-1)
  if(length(samples) > 12) stop("Please select 12 samples or less to be plotted.")
  
  if (is.numeric(samples)) { samples = colnames(dat)[samples+1] }
      
  
  plot(dat[,1], dat[,samples[1]], xlab = "% features", ylab = "% reads",
       type = "l", col = miscolores[1],...)
  
  for (i in 2:length(samples)) {
    
    lines(dat[,1], dat[,samples[i]], col = miscolores[i])
    
  }  

  legend("bottom", legend = samples, text.col = miscolores[1:length(samples)], bty = "n",
         lty = 1, lwd = 2, col = miscolores[1:length(samples)])

}
