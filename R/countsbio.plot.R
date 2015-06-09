##Counts for detected genes Plot according to BIOTYPES (boxplots)

countsbio.dat <- function (input, biotypes = NULL, factor = NULL, norm = FALSE)  {

  # input: input object
  # biotypes: List containing groups of biotypes to be studied
  # factor: If not NULL, it should contain the conditions to be studied and
  #         calculation will be done based on the mean of replicates of each condition.


  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")
  

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  
  depth = round(colSums(datos)/10^6,1); names(depth) = colnames(datos)
  
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
    
    if (norm) {
      datos = sapply(niveles, 
                     function (k) { 
                       rowMeans(as.matrix(datos[,grep(k, mifactor)]))
                     })
      datos0 = sapply(niveles, 
                      function (k) { 
                        rowMeans(as.matrix(datos0[,grep(k, mifactor)]))
                      })
      
    } else {
      datos = sapply(niveles, 
                     function (k) { 
                       10^6 * rowMeans(t(t(datos[,grep(k, mifactor)])/colSums(as.matrix(datos[,grep(k, mifactor)]))))
                     })
      datos0 = sapply(niveles, 
                      function (k) { 
                        10^6 * rowMeans(t(t(datos0[,grep(k, mifactor)])/colSums(as.matrix(datos0[,grep(k, mifactor)]))))
                      })      
    }
    colnames(datos) = colnames(datos0) = niveles
    depth = sapply(niveles, function (k) paste(range(depth[grep(k, mifactor)]), collapse = "-"))
  }
  
  
  
  # Biotypes
  if (!is.null(featureData(input)$Biotype)) {  # read biotypes if they are provided
    if (hayceros) { 
      infobio0 <- as.character(featureData(input)$Biotype)[-ceros] 
    } else { infobio0 =  as.character(featureData(input)$Biotype) }
    infobio <- as.character(featureData(input)$Biotype)
  } else { infobio0 = NULL; infobio = NULL }  
  
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
  if (is.null(bionum)) {
    resumen = vector("list", length = 1)
    names(resumen) = "global"
  } else {
    resumen = vector("list", length = length(bionum))
    names(resumen) = names(bionum)
  }
  
  cuentas = c(0,1,2,5,10)
  
  if (is.null(factor)) {
    if (norm) {
      datosCPM = datos
    } else {
      datosCPM = 10^6 * t(t(datos)/colSums(as.matrix(datos)))
    }
    
  } else { datosCPM = datos }
      
  for (i in 1:length(resumen)) { 
    
    if (i == 1) { 
      datosR = datosCPM       
    } else {       
      if(!is.null(infobio)) {
        datosR = datosCPM[which(infobio == names(resumen)[i]),]
        if (class(datosR) !=  "matrix") { datosR = t(as.matrix(datosR)) } 
      }
    }
    
    
    nfeatures = nrow(datosR)
    datosR = datosR[which(rowSums(datosR) > 0),] 
    if (class(datosR) !=  "matrix") { datosR = t(as.matrix(datosR)) }
    
    myglobal = NULL
    mypersample = NULL
    
    for (kk in 1:length(cuentas)) {
      mypersample = rbind(mypersample, apply(datosR, 2, function (x) { length(which(x > cuentas[kk])) }))
      myglobal = c(myglobal, sum(apply(datosR, 1, function (x) { max(x) > cuentas[kk] })))
    }    
    
    mypersample = round(100*mypersample/nfeatures, 1)
    mypersample = rbind(mypersample, depth)
    rownames(mypersample) = 1:nrow(mypersample)
    myglobal = c(round(100*myglobal/nfeatures, 1), nfeatures)
        
    resumen[[i]] = data.frame(c(paste("CPM >", cuentas), "depth"), mypersample, "total" = myglobal)
    colnames(resumen[[i]])[1] = names(resumen)[i]
    colnames(resumen[[i]])[2:(ncol(resumen[[i]])-1)] = colnames(datosR)
                              
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

countsbio.plot <- function (dat, samples = c(1,2), toplot = "global", 
                            plottype = c("barplot", "boxplot"), toreport = FALSE,...) {

  # dat: Data coming from countsbio.dat function
  # samples: Samples to be plotted. If NULL, all samples are plotted.
  # toplot: Name of biotype (including "global") to be plotted.

  mypar = par(no.readonly = TRUE)
  
  ## Preparing data
  if (is.null(samples)) { 
    if (NCOL(dat$result) == 1) {
      samples = 1
    } else {
      samples <- 1:NCOL(dat$result) 
    }    
  }
  if(is.numeric(toplot)) toplot = names(dat$summary)[toplot]
  if (is.numeric(samples) && !is.null(colnames(dat$result))) samples = colnames(dat$result)[samples]
  
  if (plottype == "barplot") {
    
    if ((exists("ylab") && !is.character(ylab)) || !exists("ylab")) ylab = ""
    
    datos = dat$summary[[toplot]]
    mytotal = as.numeric(datos[,"total"])
    datos = as.matrix(datos[,samples])
    rownames(datos) = as.character(dat$summary[[toplot]][,1])
    
    
    par(mar = c(6,4,4,2))
    
    barplot(as.numeric(datos[1,]), col = miscolores[1], las = 2, main = "", ylab = "", density = 70,
            ylim = c(0,100), cex.axis = 0.8, names.arg = "",...)  
    for (i in 2:(length(mytotal)-2)) {
      barplot(as.numeric(datos[i,]), col = miscolores[i], las = 2, main = "", ylab = "", add = TRUE, 
              density = 70, ylim = c(0,100), cex.axis = 0.8, names.arg = "",...)    
    }
    bp = barplot(as.numeric(datos[(length(mytotal)-1),]), col = miscolores[(length(mytotal)-1)], las = 2, 
                 main = paste(toupper(toplot), " (", mytotal[length(mytotal)], ")", sep = ""), 
                 ylab = "Sensitivity (%)", add = TRUE, names.arg = colnames(datos), cex.axis = 0.8,
                 density = 70, ylim = c(0,100), cex.names = 0.8,...) 
    for (j in 1:(length(mytotal)-1)) abline(h = mytotal[j], col = miscolores[j], lwd = 2)
    if (length(samples) <= 10) {
      mtext(side = 3, text = datos["depth",], adj = 0.5, at = bp, cex = 0.8)
    } else {
      mtext(side = 3, text = datos["depth",], at = bp, cex = 0.7, las = 2)
    }    
    legend("top", rownames(datos)[-length(mytotal)], fill = miscolores, density = 70, bty = "n", ncol = 3)
      
    par(mar = c(5, 4, 4, 4) + 0.1) 
  }
  
  
  
  
  if (plottype == "boxplot") {
    
    conteos <- as.matrix(dat$result[,samples])
    if (is.numeric(samples)) colnames(conteos) = colnames(dat$result)[samples]
    else colnames(conteos) = samples
    num <- dat$bionum[toplot]
    if (is.null(num)) {
      if (toplot == "global") {
        num = nrow(conteos)
      } else {
        num = 0
      }
    } 
    infobio = dat$biotypes
    
    if (num == 0 && toplot != "global") stop("Error: No data available. Please, change toplot parameter.")
    
    #if (!exists("ylim")) ylim = range(na.omit(log2(1+conteos)))
    if ((exists("ylab") && !is.character(ylab)) || !exists("ylab")) ylab = "Expression values"
    
    ## Plots  
    
    if (length(samples) == 1) {  # only 1 sample is to be plotted (per biotypes if available)
      
      escala = logscaling(conteos, base = 2)
      
      if (is.null(infobio)) {
        boxplot(escala$data, col = miscolores[1], ylab = ylab, #ylim = ylim,
                main = "", yaxt = "n", ...)
      } else {
        par(mar = c(10, 4, 4, 2))  
        boxplot(escala$data ~ infobio, col = miscolores, ylab = ylab, #ylim = ylim,
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
      boxplot(escala$data, col = miscolores, ylab = ylab, #ylim = ylim,
              main = main, las = 2, cex.lab = 0.9, cex.axis = 0.8, yaxt = "n", ...)       
    }
    
    axis(side = 2, at = escala$at, labels = escala$labels)
        
  }
  
  if (!toreport) par(mypar)  
  
}


