#### GC CONTENT PLOT


## Data for GC plot 

GC.dat <- function (input, factor = NULL)  {
  
  # This plot shows the mean expression for each GC content bin, globally or for each biotype (if available).
  
  # datos: Count data matrix. Each column is a different biological sample.

  if (inherits(input,"eSet") == FALSE)    
    stop("Error. You must give an eSet object\n")

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  ceros = which(rowSums(datos) == 0)
  
  if (length(ceros) > 0) {
    print(paste("Warning:", length(ceros), 
                "features with 0 counts in all samples are to be removed for this analysis."))
    datos = datos[-ceros,]
  }
      
  nsam <- NCOL(datos)
  if (nsam == 1) datos <- as.matrix(datos)
    
  
  # Per condition
  if (is.null(factor)) {  # per sample  
    print("GC content bias detection is to be computed for:")
    print(colnames(datos))
    
  } else {  # per condition
    mifactor = pData(input)[,factor]
    niveles = levels(mifactor)
    print("GC content bias detection is to be computed for:")
    print(niveles)
    datos = sapply(niveles, 
                   function (k) {
                     rowMeans(t(10^6*t(datos[,grep(k, mifactor)])/colSums(as.matrix(datos[,grep(k, mifactor)]))))
                   })
    colnames(datos) = niveles
  }
  
  
  
  # GC content
  if (any(!is.na(featureData(input)$GC)) == FALSE)
    stop ("GC content was not provided.\nPlease run addData() function to add 
          this information\n")
  
  GC <- as.numeric(as.character(featureData(input)$GC))
  if (length(ceros) > 0) GC = GC[-ceros] 
  
  
  
  # Biotypes
  if (!is.null(featureData(input)$Biotype)) {  # read biotypes if they are provided
    infobio <- as.character(featureData(input)$Biotype)
    if (length(ceros) > 0) infobio = infobio[-ceros] 
    biotypes <- unique(infobio)
    names(biotypes) <- biotypes 
    # which genes belong to each biotype
    biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
    names(biog) = biotypes
    bionum <- c(NROW(datos), sapply(biog, length))
    names(bionum) <- c("global", names(biotypes))   
    
  } else { infobio = NULL; biotypes = NULL; bionum = NULL }  
  
  
  
  ## Calculations for plot
  
  GCexpr = vector("list", length = 1 + length(biotypes))
  names(GCexpr) = c("global", names(biotypes))
  
  numXbin = 200  
  
  for (i in 1:length(GCexpr))  {
    
    if (i == 1) {  # GLOBAL
      
      numdatos = length(GC)      
      numbins = floor(numdatos / numXbin)
      misbins = quantile(GC, probs = seq(0,1,1/numbins), na.rm = TRUE)
      
      if (length(misbins) != length(unique(misbins))) {
        repes = names(table(misbins))[which(table(misbins) > 1)]
        for (rr in repes) {
          cuantos = length(which(misbins == rr))
          cuales = which(misbins == rr)
          sumo = (misbins[cuales[1]+cuantos] - misbins[cuales[1]])/cuantos
          for (j in cuales[-1]) misbins[j] = misbins[j-1] + sumo            
        }
      }
      
      miclasi = cut(GC, breaks = misbins, labels = FALSE)
      misbins = sapply(1:numbins, function (i) mean(misbins[i:(i+1)]))
      miclasi = misbins[miclasi]
      GCexpr[[i]] = aggregate(datos, by = list("GCbin" = miclasi), mean, trim = 0.025)      
      
    } else {  # PER BIOTYPE
      
      datos2 = datos[biog[[i-1]],]
      GC2 = GC[biog[[i-1]]]
      
      if (bionum[i] >= numXbin*10) {  # more than numXbin*10 genes in the biotype
        
        numdatos = length(GC2)      
        numbins = floor(numdatos / numXbin)
        misbins = quantile(GC2, probs = seq(0,1,1/numbins))
        
        if (length(misbins) != length(unique(misbins))) {
          repes = names(table(misbins))[which(table(misbins) > 1)]
          for (rr in repes) {
            cuantos = length(which(misbins == rr))
            cuales = which(misbins == rr)
            sumo = (misbins[cuales[1]+cuantos] - misbins[cuales[1]])/cuantos
            for (j in cuales[-1]) misbins[j] = misbins[j-1] + sumo            
          }
        }
        
        miclasi = cut(GC2, breaks = misbins, labels = FALSE)
        misbins = sapply(1:numbins, function (i) mean(misbins[i:(i+1)]))
        miclasi = misbins[miclasi]        
        GCexpr[[i]] = aggregate(datos2, by = list("GCbin" = miclasi), mean, trim = 0.025)        
        
      } else {   # less than numXbin*10 genes in the biotype
        
        GCexpr[[i]] = cbind(GC2, datos2)      
        
      }      
      
    }
  }
  
  
  
  ## SPLINES REGRESSION MODEL   
  #library(splines)
  
  datos = GCexpr[[1]]
  GCcont = datos[,1]  
  knots =  c(rep(GCcont[1],3), seq(GCcont[1], GCcont[length(GCcont)-1], length.out=round(length(GCcont)/10, 0)), 
             rep(GCcont[length(GCcont)], 4))
  bx = splineDesign (knots, GCcont, outer.ok = TRUE)
  
  
  mismodelos = vector("list", length = ncol(datos)-1)
  names(mismodelos) = colnames(datos)[-1]
  
  for (i in 2:ncol(datos)) {
    
    print(colnames(datos)[i])
    
    mismodelos[[i-1]] = lm(datos[,i] ~ bx)
        
    print(summary(mismodelos[[i-1]]))
    
  }
  

  
  
  ## Results
  
  list("data2plot" = GCexpr, "RegressionModels" = mismodelos)
}





#**************************************************************************#
#**************************************************************************#
#**************************************************************************#





## PLOT: Median expression for each length bin

GC.plot <- function (dat, samples = NULL, toplot = "global", toreport = FALSE,...)  {
  
  datos = dat[["data2plot"]]
  mismodelos = dat[["RegressionModels"]]
  

  if (is.null(samples)) samples <- 1:(ncol(datos[[1]])-1)
  if(length(samples) > 12) stop("Please select 12 samples or less to be plotted.")
  
  if (is.numeric(samples)) { samples = colnames(datos[[1]])[samples+1] }
  
  if (is.numeric(toplot)) {
    if (toplot == 1) { toplot = "global"} else { toplot = names(toplot)[toplot + 1] }
  }
    
  
  if ((toplot == "global") && (length(samples) <= 2)) {   ### DIAGNOSTIC PLOTS
    
    if((!toreport) && (length(samples) == 2)) par(mfrow = c(1,2))
    
    for (i in 1:length(samples)) {
      matplot(datos[[1]][,1], cbind(datos[[1]][,samples[i]], mismodelos[[samples[i]]]$fit), 
               type="pl", main=samples[i], pch=20, lty=1, lwd = 2, col = c(1,4),
               ylab = "Mean expression", xlab = "GC content bins", ylim = c(0,max(datos[[1]][,samples[i]])),...)   
      text(max(datos[[1]][,1]), 0.2*max(datos[[1]][,samples[i]]), col = 4, adj = 1,
           paste("R2 = ", 100*round(summary(mismodelos[[samples[i]]])$"r.squared",4), "%", sep = ""))
      laF = summary(mismodelos[[samples[i]]])$"fstatistic"
      text(max(datos[[1]][,1]), 0.1*max(datos[[1]][,samples[i]]), col = 4, adj = 1,
           paste("p-value:", signif(pf(laF[1], df1 = laF[2], df2 = laF[3], lower.tail = FALSE),2)))
    } 

    
  } else {   ### DESCRIPTIVE PLOTS
    
    matplot(datos[[toplot]][,1], datos[[toplot]][,samples], xlab = "GC content bins", ylab = "Mean expression", 
            type = "l", main = toupper(toplot), col = miscolores, lwd = 2, 
            ylim = range(datos[[toplot]][,-1]),lty = 1,...)
    
    legend("bottomright", samples, col = miscolores[1:length(samples)], lwd = 2, bty = "n")    
    
  }
  
  if((!toreport) && (length(samples) == 2)) layout(1)
  
   
}
