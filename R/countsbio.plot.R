##Counts for detected genes Plot according to BIOTYPES (boxplots)

countsbio.dat <- function (input, biotypes = NULL, factor = NULL)  {

  # input: input object
  # biotypes: List containing groups of biotypes to be studied
  # factor: If not NULL, it should contain the conditions to be studied and
  #         calculation will be done based on the mean of replicates of each condition.
  # quartiles: Return a summary matrix with quantile information for each biotype.



  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")
  

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  
  ceros = which(rowSums(datos) == 0)
  hayceros = (length(ceros) > 0)
    
  if (hayceros) {
    print(paste("Warning:", length(ceros), 
                "features with 0 counts in all samples are to be removed for this analysis."))
    datos0 = datos[-ceros,]
  } else { datos0 = datos}
    
  
  nsam <- NCOL(datos)
  
  if (nsam == 1) {
    datos <- as.matrix(datos)
    datos0 <- as.matrix(datos0)
  } 

  
  # Per condition
  if (is.null(factor)) {  # per sample  
    print("Count distributions are to be computed for:")
    print(colnames(datos))
    
  } else {  # per condition
    mifactor = pData(input)[,factor]
    niveles = levels(mifactor)
    print("Counts per million distributions are to be computed for:")
    print(niveles)
    datos = sapply(niveles, 
                   function (k) { 
                     rowMeans(t(t(datos[,grep(k, mifactor)])/colSums(datos[,grep(k, mifactor)])))
                     })
    datos0 = sapply(niveles, 
                    function (k) { 
                      rowMeans(t(t(datos0[,grep(k, mifactor)])/colSums(datos0[,grep(k, mifactor)])))
                    })
    colnames(datos) = colnames(datos0) = niveles
  }
  
  
  
  # Biotypes
  if (!is.null(featureData(input)$Biotype)) {  # read biotypes if they are provided
    if (hayceros) { 
      infobio0 <- as.character(featureData(input)$Biotype)[-ceros] 
    } else { infobio0 =  as.character(featureData(input)$Biotype) }
    infobio <- as.character(featureData(input)$Biotype)
  } else { infobio = NULL }  
  
  if (!is.null(infobio)) {    
    if(is.null(biotypes)) {
      biotypes <- unique(infobio)
      names(biotypes) <- biotypes    
    }   
    # which genes belong to each biotype
    biog <- lapply(biotypes, function(x) { which(is.element(infobio0, x)) })
    names(biog) = biotypes
    bionum <- c(NROW(datos0), sapply(biog, length))
    names(bionum) <- c("global", names(biotypes)) 
    bio0 = which(bionum == 0)
    if (length(bio0) > 0) bionum = bionum[-bio0]
        
  } else { biotypes = NULL; bionum = NULL }   
  
  
  # Create the summary matrix information
  resumen = vector("list", length = length(bionum))
  names(resumen) = names(bionum)
  cuentas = c(0,1,2,5,10)
      
  for (i in 1:length(resumen)) { 
    if (i == 1) { 
      datosR = datos 
    } else {       
      if(!is.null(infobio)) {
        datosR = datos[which(infobio == names(resumen)[i]),]
        if (class(datosR) !=  "matrix") { datosR = t(as.matrix(datosR)) } 
      }
    }
    
        
    numeros = t(sapply(cuentas, function(x) { apply(datosR, 2, 
                                                  function (y) { length(which(y <= x)) }) }))
    numeros = rbind(numeros, nrow(datosR)-numeros[5,], nrow(datosR), 
                    round(colSums(datos)/10^6, 1))
    resumen[[i]] = data.frame(c(0, "<=1", "<=2", "<=5", "<=10", ">10", "total", "depth"),
                              numeros)
    colnames(resumen[[i]])[1] = names(resumen)[i]
                              
  }  
  

  ## results
  cosas <- list("result" = datos0, 
                 "bionum" = bionum,
                 "biotypes" = infobio0,
                 "summary" = resumen)    
  
  cosas
}






#***************************************************************************#
#***************************************************************************#




## PLOT: Mean length for detected genes Plot according to BIOTYPES

countsbio.plot <- function (dat, toplot = "global", samples = NULL, plottype = c("barplot", "boxplot"), ...) {

  # dat: Data coming from countsbio.dat function
  # samples: Samples to be plotted. If NULL, all samples are plotted.
  # toplot: Name of biotype (including "global") to be plotted.

  ## Preparing data
  if (is.null(samples)) { samples <- 1:NCOL(dat$result) }
  if(is.numeric(toplot)) toplot = names(dat$summary)[toplot]
  if (is.numeric(samples)) samples = colnames(dat$result)[samples]
  
  if (plottype == "barplot") {
    
    if (!exists("ylab")) ylab = ""
    
    datos = dat$summary[[toplot]]
    datos = as.matrix(datos[,samples])
    rownames(datos) = dat$summary[[toplot]][,1]
    
    
    par(mar = c(6,4,4,2))
        
    barplot(as.numeric(datos["total",]), las = 2,
            cex.axis = 0.8, border = NA, ylab = "", main = "", ...)
    barplot(as.numeric(datos["<=10",]), ylab = "",  las = 2,
            cex.axis = 0.8, border = NA, add = TRUE,
            main = "", col = miscolores[1], ...)
    barplot(as.numeric(datos["<=5",]), ylab = "",  las = 2,
            cex.axis = 0.8, border = NA, add = TRUE,
            main = "", col = miscolores[2], ...)
    barplot(as.numeric(datos["<=2",]), ylab = "",  las = 2,
            cex.axis = 0.8, border = NA, add = TRUE,
            main = "", col = miscolores[3], ...)
    barplot(as.numeric(datos["<=1",]), ylab = "",  las = 2,
            cex.axis = 0.8, border = NA, add = TRUE,
            main = "", col = miscolores[4], ...)
    bp = barplot(as.numeric(datos["0",]), ylab = ylab,
                 names.arg = colnames(datos), las = 2, cex.names = 0.8,
                 cex.axis = 0.8, border = NA, add = TRUE,
                 main = paste(toupper(toplot), " (", datos["total",1], ")", sep = ""),
                 col = miscolores[6], ...)
    if (length(samples) <= 20) {
      mtext(side = 3, text = datos["depth",], adj = 0.5, at = bp, cex = 0.7)
    } else {
      mtext(side = 3, text = datos["depth",], at = bp, cex = 0.7, las = 2)
    }
    legend(x = length(samples)/3, y = datos["total",1]*0.99, c(rownames(datos)[1:5],"All"), 
           fill = c(miscolores[c(6,4:1)],"grey"), title = "#Features with expression value:",
           bg = "white", ncol = 3, cex = 0.8)
    
    par(mar = c(5, 4, 4, 4) + 0.1) 
  }
  
  
  
  
  if (plottype == "boxplot") {
    
    conteos <- as.matrix(dat$result[,samples])
    if (is.numeric(samples)) colnames(conteos) = colnames(dat$result)[samples]
    else colnames(conteos) = samples
    num <- dat$bionum[toplot]
    infobio = dat$biotypes
    
    if (num == 0) stop("Error: No data available. Please, change toplot parameter.")
    
    if (!exists("ylim")) ylim = range(na.omit(log2(1+conteos)))
    if (!exists("ylab")) ylab = "Expression values"
    
    ## Plots  
    
    if (length(samples) == 1) {  # only 1 sample is to be plotted (per biotypes if available)
      
      escala = logscaling(conteos, base = 2)
      
      if (is.null(infobio)) {
        boxplot(escala$data, col = miscolores[1], ylim = ylim, ylab = ylab,
                main = "", yaxt = "n", ...)
      } else {
        par(mar = c(10, 4, 4, 2))  
        boxplot(escala$data ~ infobio, col = miscolores, ylim = ylim, ylab = ylab,
                main = colnames(conteos), las = 2, cex.axis = 0.8, cex.lab = 0.9, yaxt = "n", ...)
        cuantos = dat$bionum[-1]
        cuantos = cuantos[sort(names(cuantos))]
        mtext(cuantos, 3, at = 1:length(cuantos), cex = 0.6, las = 2)
      }    
      
      
    } else {   # more than 1 sample is to be plotted
      if (toplot != "global") conteos = conteos[which(infobio == toplot),]          
      
      escala = logscaling(conteos, base = 2)
      main <- paste(toupper(toplot), " (", num, ")", sep = "")
      
      par(mar = c(6, 4, 2, 2))  
      boxplot(escala$data, col = miscolores, ylim = ylim, ylab = ylab,
              main = main, las = 2, cex.lab = 0.9, cex.axis = 0.8, yaxt = "n", ...)       
    }
    
    axis(side = 2, at = escala$at, labels = escala$labels)
    
    par(mar = c(5, 4, 4, 4) + 0.1)        
    
  } 
  
  
}


