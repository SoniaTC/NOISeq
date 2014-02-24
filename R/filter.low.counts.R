##########################################################################################



##***********************************************************##

## Coefficient of Variation

CV = function(data) { 100 * sd(data, na.rm = TRUE) / mean(data, na.rm = TRUE) }


##***********************************************************##



## Filtering out genes with low counts

filtered.data = function(dataset, factor, norm = TRUE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1) {
  
  dataset0 = dataset[rowSums(dataset) > 0,]
  dataset = dataset0
  
  if ((method == 3) && (norm)) {
    if (is.null(depth)) {
      stop("ERROR: Sequencing depth for each column in dataset must be provided.\n")
    }      
    dataset = t(t(dataset0) / (colSums(dataset0)/depth))  # estimate counts from normalized data
  }         
  
  if ((method < 3) && (!norm)) {
    dataset = 10^6 * t(t(dataset0) / colSums(dataset0))
  }
    
  grupos = unique(factor)
  cumple = NULL
  
  print("Filtering out low count features...")
  
  for (gg in grupos) {    
    datos = as.matrix(dataset[, grep(gg, factor)])
    
    if (method == 1) {
      
      if (ncol(datos) == 1) {
        cumplecond = (datos > cpm)
      } else {
        cumplecond = (apply(datos, 1, CV) < cv.cutoff)*(rowMeans(datos) > cpm)
        cumplecond[which(is.na(cumplecond) == TRUE)] = 0 
      }
      
      cumple = cbind(cumple, cumplecond)
    } 
    
    if (method == 2) {
      if (ncol(datos) == 1) stop("ERROR: At least 2 replicates per condition are required to apply this method.")
      mytest = apply(datos, 1, 
                     function (x) { 
                       suppressWarnings(wilcox.test(x, alternative = "greater", conf.int=FALSE, mu = 0))$"p.value" })
      cumple = cbind(cumple, 1*(mytest < 0.05))  
    }    
    
    if (method == 3) {           
      p0 = cpm / 10^6                
      mytest = apply(datos, 1, 
                         function (x) suppressWarnings(prop.test(sum(x), n=sum(datos), p = p0,
                                                alternative = "greater"))$"p.value")       
      #miproptest = p.adjust(miproptest, method = "fdr")
      cumple = cbind(cumple, 1*(mytest < 0.05))           
    }
  }
  
  cumple = which(rowSums(as.matrix(cumple)) >= 1) 
  
  print(paste(length(cumple), "features are to be kept for differential expression analysis with filtering method", method))
  
  dataset0[cumple,]
  
}




##***********************************************************##

