## Computing PCA 

PCA.dat <- function (input, norm = FALSE, logtransf = FALSE)  {

  # input: input object
  # norm: TRUE if data are already normalized, FALSE if not.
  # logtransf: TRUE if data are already log-transformed, FALSE if not.

  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")
  

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  myfactors = pData(input)
  
  
  if (!norm) datos = rpkm(datos)
  
  if (!logtransf) datos = log2(datos+1)
  
  resultat = PCA.GENES(t(datos))
  

  ## results
  resultat <- list("result" = resultat,
                   "factors" = myfactors,                   
                   "norm" = norm,
                   "logtransf" = logtransf)    
  
  resultat
}






#***************************************************************************#
#***************************************************************************#


## PLOT: PCA plot (global or per biotype)

PCA.plot <- function (dat, samples = c(1,2), plottype = "scores", factor = NULL) {

  # dat: Data coming from PCA.dat function
  # samples: Principal components to be plotted. If NULL, PC1 and PC2 (default) will be plotted.
  # toplot: Name of biotype (including "global") to be plotted.
  # plottype: One of "scores" or "loadings"
  # factor: Name of the factor to be used to color the PCA score plot. If NULL, the first one is chosen.

  
  ## Preparing data
  if (is.null(samples)) samples = 1:2     

  
  if (plottype == "loadings") { 
    
    data2plot = dat$result
        
    rango = diff(range(data2plot$loadings[,samples]))
    
    plot(data2plot$loadings[,samples], col = 1, pch = ".",
         xlab = paste("PC", samples[1], round(data2plot$var.exp[samples[1],1]*100,0), "%"),
         ylab = paste("PC", samples[2], round(data2plot$var.exp[samples[2],1]*100,0), "%"),
         main = "Loadings",
         xlim = range(data2plot$loadings[,samples]) + 0.02*rango*c(-1,1),
         ylim = range(data2plot$loadings[,samples]) + 0.02*rango*c(-1,1))    
  }
  
    
  else if (plottype == "scores") {
    
    data2plot = dat$result
    
    if (is.null(factor)) factor = 1
    myfactor = as.character(dat$factors[,factor])
        
    condis = unique(myfactor)
    
    mypch = c(17:15, 18, 8, 1, 2)  # 7
    mycolors = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)] # 12
    parapintar = data.frame("col" = rep(mycolors, 7), "pch" = rep(mypch, 12), stringsAsFactors = FALSE)         
    parapintar = parapintar[1:length(condis),]
    rownames(parapintar) = condis
    
    pch = parapintar[myfactor,"pch"]
    col = parapintar[myfactor,"col"]       
        
    rango = diff(range(data2plot$scores[,samples]))
    
    plot(data2plot$scores[,samples], col = "white", 
         xlab = paste("PC", samples[1], round(data2plot$var.exp[samples[1],1]*100,0), "%"),
         ylab = paste("PC", samples[2], round(data2plot$var.exp[samples[2],1]*100,0), "%"),
         main = "Scores",
         xlim = range(data2plot$scores[,samples]) + 0.02*rango*c(-1,1),
         ylim = range(data2plot$scores[,samples]) + 0.02*rango*c(-1,1))    
    
    points(data2plot$scores[,samples[1]], data2plot$scores[,samples[2]],
           pch = pch, col = col, cex = 1.3) 
    
    legend("topleft", condis, pch = parapintar[,"pch"], col = parapintar[,"col"], bty = "n")        
  } 
  
  
}


