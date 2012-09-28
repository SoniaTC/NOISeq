#### SATURATION PLOTS

## Data for saturation plot (with or without biotypes)

saturation.dat <- function (input, k = 0, biotypes = NULL,
                            ndepth = 5, newdetections = TRUE) {

  # input: input object. 
  # k: A feature is considered to be detected if the corresponding number of counts is > k.
  # biotypes: List containing groups of biotypes to be studied. 
  #           If biotypes = NULL, all biotypes are plotted independently.
  # ndepth: Number of different depths to be plotted.
  # newdetections: If TRUE, a second Y-axis is drawn for new detectections per million reads.

  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")

  if (any(!is.na(featureData(input)$Biotype)) == FALSE)
    stop ("No biological classification provided.\nPlease run addData() function to add 
          this information\n")

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  infobio <- featureData(input)$Biotype
  
  nsam <- NCOL(datos)

  if (!is.null(infobio) & is.null(biotypes)) {
    biotypes <- unique(infobio)
    names(biotypes) <- biotypes
  }    

  satura <- vector("list", length = length(biotypes)+1)
  names(satura) <- c("global", names(biotypes))
  
  
  # Random subsamples for each sample
  submuestras <- seq.depth <- vector("list", length = nsam)
  names(submuestras) <- names(seq.depth) <- colnames(datos)

  for (n in 1:nsam) {  # simulating subsamples for each sample

    total <- sum(datos[,n]) # total counts in sample n

    varias <- vector("list", length = ndepth) # simulation for each depth

    for (i in 1:(ndepth-1)) {    # (100% is calculated apart)

      muestra <- rmultinom(10, size = round(i*total/ndepth,0), prob = datos[,n])

      varias[[i]] <- muestra
    }

    varias[[ndepth]] <- as.matrix(datos[,n])

    submuestras[[n]] <- varias

    seq.depth[[n]] <- sapply(varias, function(x) { mean(colSums(x)) })
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



  if (newdetections) {  # computing detection increasing per million reads

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

        pendimil <- c(puntos$y[1]/puntos$x[1], pendi)*1000000

        newdet[[j]][[n]] <- pendimil

      }
    }
  } else { newdet <- NULL }


  bionum <- c(NROW(datos), sapply(biog, length))
  names(bionum) <- c("global", names(biog))


  # Results

  satura <- list("saturation" = satura, "bionum" = bionum,
                 "depth" = seq.depth, "newdet" = newdet)
  
  satura

}



##**************************************************************************#
##**************************************************************************#
##**************************************************************************#



#### Saturation plot

saturation.plot <- function(satdat, toplot = 1, samples = NULL,                            
                            ylim = NULL, yrightlim = NULL) {

  # satdat: Data coming from saturation.dat function
  # samples: Samples to be plotted. If NULL, all samples are plotted (Maximum = 12).
  # toplot: Number or name of biotype (including "global") to be plotted.
  # newdetections: If TRUE, a second Y-axis is drawn for new detections per million reads.
  # colL, colR, mybg: A vector with as many colors as different samples to be plotted.
  #                   If NULL, default colors are used.
  
  
  # Parameters
  lwdL = 2
  lwdR = 10
  xlab = "Sequencing depth (million reads)"
  ylabL = "Number of detected features"
  ylabR = "New detections per million reads"
  cex.main = cex.lab = cex.axis = cex = 1
  


  # Preparing data
  sat <- satdat$saturation[[toplot]]
  depth <- satdat$depth
  num <- satdat$bionum[[toplot]]
  new <- satdat$newdet[[toplot]]
  
  if (is.null(new) || length(samples)>2) {newdetections = FALSE} else {newdetections = TRUE}
  
  
  if (is.numeric(toplot)) {
      main <- paste(names(satdat[[1]])[toplot], " (", num, ")", sep = "")
  } else {
      main <- paste(toplot, " (", num, ")", sep = "")
  }
    


  if (is.null(samples)) {
    samples <- 1:length(sat)
  }
  
  legend = names(satdat$saturation[[1]])[samples]

  
  # colors
  miscolores <- colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)] 
  
  if (length(samples) > 2) {
      colL <- miscolores
  }  else {
      colL <- miscolores[c(4,2)]
  }  

  colR <- miscolores[c(12,11)] 

  mybg <- colL
    
  
  # yleftlim for plot and plot.y2
  if (is.null(ylim)) {
    yleftlim <- range(unlist(sat[samples]))
  } else { yleftlim = ylim }

  # xlim for plot
  xlim <- range(unlist(depth[samples])/10^6)
  

  # Percentage of detections at maximum depth
  percen <- sapply(sat[samples], function(x) { round(100*max(x)/num, 1) })


  # Drawing new detections bars?
  if (newdetections) {

    if (length(samples) <= 2) {
      bars <- TRUE
    } else {
      bars <-  FALSE
      print("WARNING: Too many samples to plot new detection bars. Maximum number of samples = 2.")
    }

  } else { bars <- FALSE }


  
  ## PLOTS
  
  if(!bars) {  # PLOT for detections without bars for new detections

    plot(depth[[samples[1]]]/10^6, sat[[samples[1]]], pch = 21, col = colL[1],
         ylim = yleftlim, lwd = lwdL,
         xlim = xlim, main = main, type = "o", xlab = xlab,
         ylab = ylabL, cex.main = cex.main, cex.lab = cex.lab,
         cex.axis = cex.axis, bg = mybg[1])

    if (length(samples) > 1) {  # for more than 1 sample

      j <- 2

      for (i in samples[-1]) {
        lines(depth[[i]]/10^6, sat[[i]], pch = 21, bg = mybg[j], col = colL[j],
              lwd = lwdL, type = "o")
        j <- j+1
      }
    }

    legend("top", legend = paste(legend, ": ", percen, "% detected", sep = ""),
           pch = 21, pt.bg = mybg, text.col = colL, bty = "n",
           lwd = lwdL, col = colL, cex = cex)

    

  } else {   # PLOT for detections and new detections

    # yrightlim for plot.y2
    if (is.null(yrightlim)) {
      yrightlim <- c(0, max(10,max(na.omit(unlist(new[samples])))))
    }
    
    
    nf <- layout(matrix(c(1,2),2,1,byrow=TRUE),heights=c(0.8,0.2))
    par(mar = c(5, 4, 4, 4) + 0.1)


    # PLOT with 2 axis
    plot.y2(x = depth[[samples[1]]]/10^6, yright = new[[samples[1]]],
            yleft = sat[[samples[1]]], type = c("h", "o"),
            lwd = c(lwdR, lwdL), xlab = xlab, xlim = xlim,
            yrightlim = yrightlim, yleftlim = yleftlim,
            yylab = c(ylabR, ylabL), pch = c(1,21), col = c(colR[1],colL[1]),
            main = main, x2 = depth[[samples[2]]]/10^6,
            yright2 = new[[samples[2]]], yleft2 = sat[[samples[2]]],
            col2 = c(colR[2],colL[2]), cex.main = cex.main, bg = mybg,
            cex.lab = cex.lab, cex.axis = cex.axis, cex = cex)

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
    par(mar = c(5, 4, 4, 4) + 0.1)
    layout(1)
    
  }
}
