##Counts for detected genes Plot according to BIOTYPES (boxplots)

countsbio.dat <- function (input, cols = NULL, k = 0, biotypes = NULL,
                           ndepth = 5, quartiles = FALSE)  {

  # input: input object
  # cols: Data columns to be used in the analysis
  # k: Only genes with counts > k will be included in the analysis
  # biotypes: List containing groups of biotypes to be studied
  # ndepth: Number of different depths to be plotted. 
  # quartiles: Return a summary matrix with quantile information for each biotype

  # if no cols defined, all the samples are taken

  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")

  if (any(!is.na(featureData(input)@data$Biotype)) == FALSE)
    stop ("No biological classification provided.\nPlease run addData() function to add 
          this information\n")

  if (is.null(cols)) {
    if (!is.null(assayData(input)$exprs))
      cols <- c(1:NCOL(assayData(input)$exprs))
    else
      cols <- c(1:NCOL(assayData(input)$counts))
  }

  if (!is.null(assayData(input)$exprs))
    datos <- assayData(input)$exprs[,cols]
  else
    datos <- assayData(input)$counts[,cols]
  
  if (length(cols) == 1)
    datos <- as.matrix(datos)

  infobio <- featureData(input)@data$Biotype

  nsam <- NCOL(datos)

  if (is.null(biotypes)) {
    biotypes <- unique(infobio)
    names(biotypes) <- biotypes
  }
  
  # which genes belong to each biotype
  biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
  
  satura <- vector("list", length = length(biotypes)+1)
  names(satura) <- c("global", names(biotypes))
  
  # Random subsamples for each sample at each depth
  submuestras <- seq.depth <- vector("list", length = nsam)
  names(submuestras) <- names(seq.depth) <- colnames(datos)

  for (n in 1:nsam) {

    total <- sum(datos[,n])

    varias <- vector("list", length = ndepth) 
        
    for(i in 1:(ndepth-1)) {  

      varias[[i]] <- rmultinom(10, size = round(total*i/ndepth, 0),
                               prob = datos[,n])
    }

    varias[[ndepth]] <- as.matrix(datos[,n])

    submuestras[[n]] <- lapply(varias, rowMeans)

    seq.depth[[n]] <- sapply(varias, function (x) { mean(colSums(x)) })
  }
  
  ## selecting detected genes counts for each biotype
  for (j in 2:length(satura))  {
    satura[[j]] <- vector("list", length = nsam)
    names(satura[[j]]) <- colnames(datos)

    for (n in 1:nsam) {      

      # selecting genes in bioclass j for each sample
      conbio <- lapply(submuestras[[n]],
                       function(x) { x[biog[[j-1]]] })

      satura[[j]][[n]] <- vector("list", length = ndepth)
      names(satura[[j]][[n]]) <- paste("depth", 1:ndepth, sep = "")      

      for (i in 1:length(conbio)) {            
      # selecting the genes with counts > k
        noK <- noceros(conbio[[i]], k = k, num = FALSE)

        if (is.null(noK)) {
          satura[[j]][[n]][[i]] <- NA
        } else {
          satura[[j]][[n]][[i]] <- conbio[[i]][noK]
        }
      }      
    }
  }


  ## Global
  satura[[1]] <- vector("list", length = nsam)
  names(satura[[1]]) <- colnames(datos)

  for (n in 1:nsam) {

    satura[[1]][[n]] <- vector("list", length = length(biotypes))
    names(satura[[1]][[n]]) <- names(biotypes)

    for (j in 1:length(biotypes)) {

      satura[[1]][[n]][[j]] <- satura[[j+1]][[n]][[ndepth]]

    }
  }


  bionum <- c(NROW(datos), sapply(biog, length))
  names(bionum) <- names(satura)

  # Create the summary matrix information
  if (quartiles == TRUE) {
    
    annot <- names(biotypes)

    mat <- vector("list", length = length(cols))
    names(mat) <- names(satura[[1]])

    for (h in 1:length(mat)) {

      mat[[h]] <- matrix(nrow = length(biotypes), ncol = 3)

      for ( i in 1:length(biotypes)) {
        mat[[h]][i,] = quantile(satura[[1]][[h]][[i]],na.rm=TRUE)[2:4]
      }

      rownames(mat[[h]]) <- annot
      colnames(mat[[h]]) <- c("1st quartile","Median","3rd quartile")
    }
  }

  ## results
  if (quartiles == TRUE)
    satura <- list("result" = satura, 
		  "bionum" = bionum,
		  "depth" = seq.depth,
		  "quart" = mat)                 
  else
    satura <- list("result" = satura, 
		"bionum" = bionum,
		"depth" = seq.depth)                 

  satura
}




#***************************************************************************#
#***************************************************************************#




miscolores <- colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]



#***************************************************************************#
#***************************************************************************#




## PLOT: Mean length for detected genes Plot according to BIOTYPES

countsbio.plot <- function (dat, toplot = 1, samples = NULL, ylim = NULL,...) {

  # dat: Data coming from countsbio.dat function
  # samples: Samples to be plotted. If NULL, two first samples are plotted.
  # toplot: Number or name of biotype (including "global") to be plotted.

  
  ## Preparing data
  
  legend = names(dat$result[[1]])[samples]
  
  if (is.null(samples)) { samples <- 1:2 }
  sat <- dat$result[[toplot]][samples]
  depth <- dat$depth[samples]
  num <- dat$bionum[[toplot]]

  if (num == 0) {

    print("Error: No data available. Please, change toplot parameter.")
     
  } else {

    # ylim for plot
    if (is.null(ylim)) {
      ylim <- range(na.omit(unlist(sat)))   
    } 
    
    
    # main
    if (is.numeric(toplot)) {
      
      elbiotipo = names(dat$result)[toplot]
      
    } else { elbiotipo = toplot }
    
    main <- paste(names(sat), "-", elbiotipo, " (", num, ")", sep = "")
    
      
      
    
    # BOXPLOTSdim
    
    if (length(samples) == 1) {  # only 1 sample

      if (toplot == 1 | toplot == "global") { # global, only 1 sample
     
        boxplot(sat[[1]], col = miscolores[12], ylim = ylim, main = main,
                xlab = NULL, ylab = "counts of detected features", cex.main = 1,
                cex.lab = 1, cex.axis = 1, las = 0,...)

        cuantos <- sapply(lapply(sat[[1]], na.omit), length)

        mtext(cuantos, 3, at = 1:length(sat[[1]]), cex = 0.6, las = 0)
      }
      

      if (toplot != 1 & toplot != "global") { # biotype, only 1 sample       

        m <- boxplot(sat[[1]], col = miscolores[12], ylim = ylim, main = main,
                xlab = NULL, ylab = "counts of detected features", cex.main = 1,
                cex.lab = 1, cex.axis = 1,
                names = round(depth[[1]]/10^6,2), las = 0,...)

        cuantos <- sapply(lapply(sat[[1]],na.omit), length)

        mtext(cuantos, 3, at = 1:length(sat[[1]]), cex = 0.6, las = 0)
      }
      

    } else {   # 2 samples

      par (mfrow = c(1,2))
      mycolor <- miscolores[c(12,11)]

      if (toplot == 1 | toplot == "global") { # global, 2 samples
     
        boxplot(sat[[1]], col = mycolor[1], ylim = ylim,
                main = main[1], xlab = NULL, ylab = "counts of detected features", cex.main = 1,
                cex.lab = 1, cex.axis = 1, las = 0,...)

        cuantos <- sapply(lapply(sat[[1]], na.omit), length)

        mtext(cuantos, 3, at = 1:length(sat[[1]]), cex = 0.6, las = 0)

        boxplot(sat[[2]], col = mycolor[2], ylim = ylim,
                main = main[2], xlab = NULL, ylab = "counts of detected features", cex.main = 1,
                cex.lab = 1, cex.axis = 1, las = 0,...)

        cuantos <- sapply(lapply(sat[[2]], na.omit), length)

        mtext(cuantos, 3, at = 1:length(sat[[2]]), cex = 0.6, las = 0)
      }
      
      
      if (toplot != 1 & toplot != "global") { # biotype, 2 samples

        boxplot(sat[[1]], col = mycolor[1], ylim = ylim,
                main = main[1], xlab = NULL, ylab = "counts of detected features", cex.main = 1,
                cex.lab = 1, cex.axis = 1,
                names = round(depth[[1]]/10^6,2), las = 0,...)

        cuantos <- sapply(lapply(sat[[1]], na.omit), length)

        mtext(cuantos, 3, at = 1:length(sat[[1]]), cex = 0.6, las = 0)

        boxplot(sat[[2]], col = mycolor[2], ylim = ylim,
                main = main[2], xlab = NULL, ylab = "counts of detected features", cex.main = 1,
                cex.lab = 1, cex.axis = 1, names = round(depth[[2]]/10^6,2), las = 0,...)

        cuantos <- sapply(lapply(sat[[2]], na.omit), length)

        mtext(cuantos, 3, at = 1:length(sat[[2]]), cex = 0.6, las = 0)    
      }
    }
  }
  par(mfrow=c(1,1))

}
