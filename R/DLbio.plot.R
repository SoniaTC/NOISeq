#### GENE LENGTH PLOTS


## Data for gene length plot (with biotypes)

DLbio.dat <- function (input, k = 0, biotypes = NULL, ndepth = 5)  {
  
  # This plot shows the mean length for detected genes for each biotype.
  
  # datos: Count data matrix. Each column is a different biological sample.
  # k: A feature is considered to be detected if the corresponding number of counts is > k.
  # infobio: Vector containing biotype for each gene in "datos".
  # biotypes: List containing groups of biotypes to be studied.
  # ndepth: Number of different depths to be plotted.
  # long: Feature length for each feature in datos.

  if (inherits(input,"eSet") == FALSE)    
    stop("Error. You must give an eSet object\n")

  if (any(!is.na(featureData(input)$Biotype)) == FALSE)
    stop ("No biological classification provided.\nPlease run addData() function to add 
          this information\n")

  if (any(!is.na(featureData(input)$Length)) == FALSE)
    stop ("Feature length was not provided.\nPlease run addData() function to add 
          this information\n")
  
  long <- as.numeric(as.character(featureData(input)$Length))
  infobio <- featureData(input)$Biotype

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs
  else
    datos <- assayData(input)$counts
  
  nsam <- NCOL(datos)

  if (is.null(biotypes) == TRUE) {

    biotypes <- as.list(unique(infobio))
    names(biotypes) <- unique(infobio)

  }
  

  satura <- vector("list", length = length(biotypes))
  names(satura) <- names(biotypes)

  # Random subsamples for each sample
  submuestras <- seq.depth <- vector("list", length = nsam)
  names(submuestras) <- names(seq.depth) <- colnames(datos)

  for (n in 1:nsam) {

    total <- sum(datos[,n])

    varias <- vector("list", length = ndepth)

    for (i in 1:(ndepth-1)) {

      muestra <- rmultinom(10, size = round(i*total/ndepth, 0), prob = datos[,n])

      varias[[i]] <- muestra
    }

    varias[[ndepth]] <- as.matrix(datos[,n])

    submuestras[[n]] <- varias

    seq.depth[[n]] <- sapply(varias, function (x) { mean(colSums(x)) })
  }
  

  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })

  
  for (j in 1:length(satura)) { # for each biotype

    satura[[j]] <- vector("list", length = nsam)
    names(satura[[j]]) <- colnames(datos)

    longbio <- long[biog[[j]]]

    for (n in 1:nsam) { # for each biological sample

      conbio <- lapply(submuestras[[n]],
                       function (x) { as.matrix(x[biog[[j]],]) })
      
      conbio.0 <- lapply(conbio, function(x) {
        apply(x, 2, function(y) {
          median(longbio[noceros(y, k = k, num = FALSE)], na.rm = TRUE) })})                       

      satura[[j]][[n]] <- sapply(conbio.0, mean, na.rm = TRUE)
      
    }
  }

 

  ## computing biotype median length
  #length.biotype <- sapply(biog, function (x) { median(na.omit(long[x])) })
  length.biotype <- NULL
  
  total <- rowSums(datos)
  totno0 <- noceros(total, num = FALSE, k = 0)

  long.mayor150 <- which(long > 150)
  long.menor150 <- which(long <= 150)

  for (i in 1:length(biog)) {
    det.menor150 <- int.mult(list(biog[[i]], totno0, long.menor150))
    mayor150 <- intersect(long.mayor150, biog[[i]])
    estos <- union(det.menor150, mayor150)
    longestos <- na.omit(long[estos])
    length.biotype <- c(length.biotype, median(longestos))
  }

  names(length.biotype) <- names(biotypes)
  
  satura <- list("result" = satura, 
                 "bionum" = sapply(biog, length),
                 "depth" = seq.depth, 
                 "biolength" = length.biotype)
  satura
}





#**************************************************************************#
#**************************************************************************#
#**************************************************************************#





## PLOT: Mean length for detected genes Plot according to BIOTYPES

DLbio.plot <- function (dat, samples = NULL, toplot = "protein_coding",
                        ylim = NULL,...)  {

  if (is.null(samples)) {
    samples <- 1:length(dat$depth)
  }
  
  legend = names(dat$result[[1]])[samples]


  # Preparing data
  sat <- dat$result[[toplot]][samples]
  depth <- dat$depth[samples]
  biolong <- dat$biolength[[toplot]]
  num <- dat$bionum[[toplot]]

  main <- paste(names(dat$result[toplot]), " (", num, ")", sep = "")
  
  # colors
  mycol <- colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
      
  

  # ylim for plot
  if (is.null(ylim)) {

    ylim <- range(c(na.omit(unlist(sat)), biolong))
    ylim <- ylim + 0.1 * diff(ylim) * c(-1,1)    
  }  

  # xlim for plot
  xlim <- range(unlist(depth))/10^6
  
  
  # PLOT
  if(is.na(biolong)) {
    
    plot(1:5, 1:5, type = "n", axes = FALSE, main = main, xlab = "", ylab = "", ...)
    text(3, 4, "Biotype not found in the dataset", adj = 0.5, cex = 1, font = 2)
    
  } else {

    if( diff(ylim) < 100 ) {   # correcting ylim
      ylim <- mean(ylim) + 50*c(-1,1)
      ylim[1] <- max(ylim[1], 0)
    }      

    plot(depth[[1]]/10^6, sat[[1]], pch = 19, col = mycol[1],
         ylim = ylim, xlim = xlim, main = main, type = "o", 
         xlab = "Sequencing depth (million reads)", 
         ylab = "Median length of detected features",
         cex.main = 1, cex.lab = 1, cex.axis = 1,...)

    if (length(samples) > 1) {

      for (i in 2:length(samples)) {

        lines(depth[[i]]/10^6, sat[[i]], pch = 19, col = mycol[i], type = "o",...)

      }
    }

    abline(h = biolong, lty = 2, col = "grey")

    text(mean(xlim), biolong + 0.02 * diff(ylim),
         "median global length", col = "grey", cex = 1)

    legend("top", legend = legend, text.col = mycol[1:length(samples)], bty = "n",
           lty = 1, lwd = 2, col = mycol[1:length(samples)], cex = 1, ncol = 2)
  }
}
