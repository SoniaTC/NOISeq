
##### Plot to compare count distributions for two or more samples


### Generating data

cd.dat <- function (input, norm = FALSE, refColumn = 1) {
  
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
  
  
  ceros = which(rowSums(datos) == 0)
  hayceros = (length(ceros) > 0)
  
  if (hayceros) {
    print(paste("Warning:", length(ceros), 
                "features with 0 counts in all samples are to be removed for this analysis."))
    datos = datos[-ceros,]
  }  
  
  
  ## scaling data and/or changing 0 to k
  if (norm) {
    datos = sinceros(datos, k = NULL)
  } else { 
    datos = rpkm(datos, long = 1000, lc = 1, k = 0.5)
  }
  
  
  
  ## to plot
  data2plot = log2(datos / datos[,refColumn])
  
  if (is.numeric(refColumn)) refColumn = colnames(datos)[refColumn]
  print(paste("Reference sample is:", refColumn))
  
  
  #### Diagnostic test
    
  MsinRef = as.matrix(data2plot[,-match(refColumn, colnames(data2plot))])
  colnames(MsinRef) = colnames(data2plot)[-match(refColumn, colnames(data2plot))]
  alpha = 0.05
  alpha = alpha/ncol(MsinRef)
  nperm = 10^3
    
  bootmed = sapply(1:nperm, function(k) {
    permut = sample(1:nrow(MsinRef), replace = TRUE, nrow(MsinRef))
    permut = as.matrix(MsinRef[permut,])
    permut = apply(permut, 2, median)    
    permut
  })
  
  bootmed = t(apply(bootmed, 1, quantile, probs = round(c(alpha/2, 1 - alpha/2), 4)))
  diagno = apply(bootmed, 1, 
                 function (x) { 
                   ddd =  (x[1] <= 0) * (0 <= x[2])
                   if (ddd == 1) { ddd = "PASSED" } else { ddd = "FAILED"}
                   ddd
                   })
  bootmed = cbind(bootmed, diagno)
    
  rownames(bootmed) = colnames(MsinRef)
  colnames(bootmed)[3] = "Diagnostic Test"
  print("Confidence intervals for median of M:")
  print(bootmed)
    
  if ("FAILED" %in% bootmed[,3]) {
    print("Diagnostic test: FAILED. Normalization is required to correct this bias.")
  }  else {
    print("Diagnostic test: PASSED.")
  }
  
  
  #### Results
  
  list("data2plot" = data2plot, "refColumn" = refColumn, "DiagnosticTest" = bootmed)
   
}


###########################################################################
###########################################################################
###########################################################################




### Generating plot

cd.plot <- function (dat, samples = NULL,...) {
  
  refColumn = dat$refColumn
  dat = dat$data2plot
  
  if (is.null(samples)) samples <- 1:ncol(dat)
  
  if (is.numeric(samples)) { samples = colnames(dat)[samples] }
  
  samples = setdiff(samples, refColumn)
  
  if(length(samples) > 12) stop("Please select 12 samples or less to be plotted (excluding reference).")
      
  dat = dat[,samples]
  
  dat.dens = apply(dat, 2, density, adjust = 1.5)
  limY = c(0,max(sapply(dat.dens, function (x) max(x$y, na.rm = TRUE))))
  
  plot(dat.dens[[1]], xlab = "M = log2(sample/refsample)", ylab = "Density", lwd = 2, ylim = limY,
       type = "l", col = miscolores[1], main = paste("Reference sample:", refColumn), ...)
  abline(v = median(dat[,1], na.rm = TRUE), col = miscolores[1], lty = 2)
  
  for (i in 2:length(samples)) {
    lines(dat.dens[[i]], col = miscolores[i], lwd = 2)
    abline(v = median(dat[,i], na.rm = TRUE), col = miscolores[i], lty = i+1)
  }

  legend("topleft", legend = samples, text.col = miscolores[1:length(samples)], bty = "n",
         lty = 1, lwd = 2, col = miscolores[1:length(samples)])

}
