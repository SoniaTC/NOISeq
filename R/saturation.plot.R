#### SATURATION PLOTS

## Data for saturation plot (with or without biotypes)

saturation.dat <- function (input, k = 0, biotypes = NULL, ndepth = 6) {

  # input: input object. 
  # k: A feature is considered to be detected if the corresponding number of counts is > k.
  # biotypes: List containing groups of biotypes to be studied. 
  #           If biotypes = NULL, all biotypes are plotted independently.
  # ndepth: Number of different depths to be plotted.


  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  if (!is.null(featureData(input)$Biotype)) {  # read biotypes if they are provided
    infobio <- featureData(input)$Biotype    
  } else { infobio = NULL }
  
    
  nsam <- NCOL(datos)

  if (!is.null(infobio)) {    
    if(is.null(biotypes)) {
      biotypes <- unique(infobio)
      names(biotypes) <- biotypes         
    }   
  } else { biotypes = NULL }   

  
  
  satura <- vector("list", length = length(biotypes)+1)
  names(satura) <- c("global", names(biotypes))
  
  ndepth1 = ceiling(ndepth/2)
  
  datos = round(datos, 0)
#   datos0 = datos
#   datos0[datos0 == 0] = 0.05
  datos0 = datos + 0.2
    
  
  # Random subsamples for each sample
  submuestras <- seq.depth <- vector("list", length = nsam)
  names(submuestras) <- names(seq.depth) <- colnames(datos)

  for (n in 1:nsam) {  # simulating subsamples for each sample

    total <- sum(datos[,n]) # total counts in sample n

    varias <- vector("list", length = ndepth+1) # simulation for each depth and real depth    

    for (i in 1:(ndepth1)) {    # simulating depths < real depth

      muestra <- rmultinom(10, size = round(total/(ndepth1+1),0), prob = datos[,n])
      
      if (i == 1) { 
        varias[[i]] <- muestra 
      } else { varias[[i]] <- varias[[i-1]] + muestra }
    }

    varias[[ndepth1+1]] <- as.matrix(datos[,n])
    
    for (i in (ndepth1+2):(ndepth+1)) {    # simulating depths < real depth
      
      muestra <- rmultinom(10, size = round(total/(ndepth1+1),0), prob = datos0[,n])
      
      if (i == ndepth1+2) { 
        varias[[i]] <- matrix(varias[[i-1]], ncol = 10, nrow = nrow(varias[[i-1]])) + muestra 
      } else { varias[[i]] <- varias[[i-1]] + muestra }
    }    

    submuestras[[n]] <- varias

    seq.depth[[n]] <- c(round(total/(ndepth1+1),0)*(1:ndepth1), total, 
                        round(total/(ndepth1+1),0)*((ndepth1+2):(ndepth+1)))
  }


  
  # Global saturation

  satura[[1]] <- vector("list", length = nsam)
  names(satura[[1]]) <- colnames(datos)

  for (n in 1:nsam) { # for each sample

    satura[[1]][[n]] <- sapply(submuestras[[n]],
                               function(x) { mean(apply(x, 2, noceros, k = k)) })                               
  }
  

  # Per biotypes
  if (!is.null(infobio)) { # if biotypes available

    biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
    names(biog) = names(biotypes)

    for (j in 2:length(satura))  { # for each biotype

      satura[[j]] <- vector("list", length = nsam)
      names(satura[[j]]) <- colnames(datos)

      for (n in 1:nsam) { # for each sample

        conbio <- lapply(submuestras[[n]], function(x) { as.matrix(x[biog[[j-1]],]) })

        satura[[j]][[n]] <- sapply(conbio, function(x) { mean(apply(x, 2, noceros, k = k)) })
      }
    }
  } else { biog <- NULL }



  # computing detection increasing per million reads
  newdet <- vector("list", length = 1+length(biotypes))
  names(newdet) <- c("global", names(biotypes))
  
  for (j in 1:length(newdet))  {
    
    newdet[[j]] <- vector("list", length = nsam)
    names(newdet[[j]]) <- colnames(datos)
    
    for (n in 1:nsam) {      
      puntos <- data.frame("x" = seq.depth[[n]], "y" = satura[[j]][[n]])      
      pendi <- NULL
      
      for(i in 2:nrow(puntos)) {        
        pendi <- c(pendi,
                   (puntos$y[i]-puntos$y[i-1])/(puntos$x[i]-puntos$x[i-1]))
      }
      
      newdet[[j]][[n]] <- c(NA,pendi*1000000)
      
    }
  }
  



  bionum <- c(NROW(datos), sapply(biog, length))
  names(bionum) <- c("global", names(biog))
  
  
  # Results at real sequencing depth
  real = vector("list", length = length(satura))
  names(real) = names(satura)
  realdepth = sapply(seq.depth, function (x) x[ndepth1+1])/10^6
  
  for (i in 1:length(real)) {
    real[[i]] = data.frame("depth" = realdepth,
                           "detec" = sapply(satura[[i]], function (x) x[ndepth1+1]))
    rownames(real[[i]]) = colnames(datos)
  }
  

  # Results

  satura <- list("saturation" = satura, "bionum" = bionum,
                 "depth" = seq.depth, "newdet" = newdet, "real" = real)
  
  satura

}






##**************************************************************************#
##**************************************************************************#
##**************************************************************************#



#### Saturation plot

saturation.plot <- function(satdat, samples = NULL, toplot = 1, 
                            yrightlim = NULL, toreport = FALSE, yleftlim = NULL, ...) {

  # satdat: Data coming from saturation.dat function
  # samples: Samples to be plotted. If NULL, all samples are plotted (Maximum = 12).
  # toplot: Number or name of biotype (including "global") to be plotted.
  # colL, colR, mybg: A vector with as many colors as different samples to be plotted.
  #                   If NULL, default colors are used.
  
  mypar = par(no.readonly = TRUE)
  
  # Parameters
  lwdL = 2
  lwdR = 10
  xlab = "Sequencing depth (million reads)"
  ylabL = "Number of detected features"
  ylabR = "New detections per million reads"
  cex.main = cex.lab = cex.axis = 1
  cex = 0.8
  
  
  if (is.null(samples)) {
    samples <- 1:length(sat)
  }  


  # Preparing data
  sat <- satdat$saturation[[toplot]]
  depth <- satdat$depth
  num <- satdat$bionum[[toplot]]
  nuevo <- satdat$newdet[[toplot]]
  real = satdat$real[[toplot]][samples,]
    
  
  if (is.numeric(toplot)) {
      main <- paste(toupper(names(satdat[[1]])[toplot]), " (", num, ")", sep = "")
  } else {
      main <- paste(toupper(toplot), " (", num, ")", sep = "")
  }

  
  legend = names(satdat$saturation[[1]])[samples]
  if (toreport) legend = samples

  
  # colors
  miscolores <- colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)] 
  
  if (length(samples) > 2) {
      colL <- miscolores
  }  else {
      colL <- miscolores[c(4,2)]
  }  

  colR <- miscolores[c(12,11)] 

  mybg <- colL

  

  # xlim for plot
  xlim <- range(unlist(depth[samples])/10^6)
  
  # yleftlim
  if (is.null(yleftlim)) {
    yleftlim <- range(unlist(sat[samples]))
  } else { yleftlim = yleftlim }


  # Percentage of detections at real depth
  percen <- round(100*real[,"detec"]/num, 1)


  # Drawing new detections bars    
  if (length(samples) <= 2) {
    bars <- TRUE
  } else {
    bars <-  FALSE
  } 


  
  ## PLOTS
  
  if(!bars) {  # PLOT for detections without bars for new detections

    plot(depth[[samples[1]]]/10^6, sat[[samples[1]]], pch = 21, col = colL[1], #bg = mybg[1],
         lwd = lwdL, ylim = yleftlim,
         xlim = xlim, main = main, type = "b", xlab = xlab,
         ylab = ylabL, cex.main = cex.main, cex.lab = cex.lab,
         cex.axis = cex.axis, ...)

    if (length(samples) > 1) {  # for more than 1 sample

      j <- 2

      for (i in samples[-1]) {
        lines(depth[[i]]/10^6, sat[[i]], pch = 21, col = colL[j], #bg = mybg[j], 
              lwd = lwdL, type = "b")
        j <- j+1
      }
    }
    
    points(real, pch = 21, col = colL[1:length(samples)], bg = mybg[1:length(samples)])

    legend("bottom", legend = paste(legend, ": ", percen, "% detected", sep = ""),
           pch = 21, pt.bg = mybg, text.col = colL, bty = "n", ncol = 2,
           lwd = lwdL, col = colL, cex = cex)

    

  } else {   # PLOT for detections and new detections

    # yrightlim for plot.y2
    if (is.null(yrightlim)) {
      yrightlim <- c(0, max(10,max(na.omit(unlist(nuevo[samples])))))
    }
    
    
        
    if (!toreport) nf <- layout(matrix(c(1,2),2,1,byrow=TRUE),heights=c(0.8,0.2))
    
      
    par(mar = c(5, 4, 4, 4) + 0.1)


    # PLOT with 2 axis
    plot.y2(x = depth[[samples[1]]]/10^6, yright = nuevo[[samples[1]]],
            yleft = sat[[samples[1]]], type = c("h", "o"),
            lwd = c(lwdR, lwdL), xlab = xlab, xlim = xlim,
            yrightlim = yrightlim, yleftlim = yleftlim,
            yylab = c(ylabR, ylabL), pch = c(1,21), col = c(colR[1],colL[1]),
            main = main, x2 = depth[[samples[2]]]/10^6,
            yright2 = nuevo[[samples[2]]], yleft2 = sat[[samples[2]]],
            col2 = c(colR[2],colL[2]), cex.main = cex.main, #bg = mybg,
            cex.lab = cex.lab, cex.axis = cex.axis, cex = cex, ...)
    
    points(real, pch = 21, col = colL, bg = mybg)

    par(mar = c(0,0,0,0))
    plot(0,axes=FALSE,type="n")
    
    #HEADERS
    rect(0.7,-0.7, 1.3, 1, col = "grey90", border = "grey90") 
    text(0.93, 0.7,"Left axis", font = 3,cex=1.2)
    text(1.07, 0.68, "Right axis", font = 3, cex = 1.2)
    text(1.22,0.7,"%detected", font = 3, cex = 1.2)
    
    # The rest of the legend arguments
    text(0.72,0.15,legend[1], font = 2, adj=0)
    points(0.93, 0.15, lty = 1, pch = 21,
           col = colL[1], bg = mybg[1])
    points(1.07, 0.15, pch = "-",
           col = colR[1], cex = lwdR)
    text(1.24, 0.15, percen[1],adj=1)

    if (length(samples) == 2) {

      text(x = 0.72,-0.25, legend[2], font = 2, adj=0)
      points(0.93,-0.25, lty = 1, pch = 21,
             col = colL[2], bg = mybg[2])
      points(1.07, -0.25, pch = "-",
             col = colR[2], cex = lwdR)
      text(1.24, -0.25, percen[2],adj=1)
    }
    # Reset with the default values    
    if (!toreport) par(mypar); layout(1) 
    
  }
}
