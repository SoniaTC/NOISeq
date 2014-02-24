allMD <- function (input, factor, conditions, k = 0.5, replicates, norm = "rpkm",
            pnr = 0.2, nss = 5, v = 0.02, lc = 0)

# input:  Set of data of type Input

# conditions: Levels of the factor to be compared (when the factor has more than 2 levels)
  
# k:      When counts = 0, 0 will be changed to k. By default, k = 0.5.
  
# norm:   Normalization method. It can be one of "rpkm" (default), "uqua"
#         (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).

# pnr:    Percentage of total reads (seq.depth) for each simulated sample.
#         Only needed when noise = "simul". By default, pnr = 1.

# nss:    Number of simulated samples (>= 2). By default, nss = 5.
#         If nss = 0, real samples are used to compute noise.

# v:      Variability in sample total reads used to simulate samples.
#         By default, v = 0.02. Sample total reads is computed as a
#         random value from a uniform distribution in the interval
#         [(pnr-v)*sum(counts), (pnr+v)*sum(counts)]

# lc:     Length correction in done by dividing expression by length^lc.
#         By default, lc = 1. 

{

#   n1 <- ncol(as.matrix(datos1))
#   n2 <- ncol(as.matrix(datos2))
#   g.sinL <- names(which(is.na(long)))
  
  # Check if the factor introduced is already defined
  # If the factor introduced is defined and has more than 2 conditions, it will check if the conditions specified are defined too
  condition_fac = FALSE
  condition_lev = FALSE
  
  datos1 <- datos2 <- matrix()

  for (i in colnames(pData(input))) {
    if (factor == i) {
      condition_fac = TRUE
      if (length(levels(pData(input)[,i])) == 2) {

        if (!is.null(assayData(input)$exprs)) {
          datos1 <- assayData(input)$exprs[,which(pData(input)[,i] ==levels(pData(input)[,i])[1])]
          datos2 <- assayData(input)$exprs[,which(pData(input)[,i] ==levels(pData(input)[,i])[2])]
        } else {
          datos1 <- assayData(input)$counts[,which(pData(input)[,i] ==levels(pData(input)[,i])[1])]
          datos2 <- assayData(input)$counts[,which(pData(input)[,i] ==levels(pData(input)[,i])[2])]
        }

        # Define the comparison string
        comparison <- paste(levels(pData(input)[,i])[1], levels(pData(input)[,i])[2], sep=" - ")

        condition_lev = TRUE
      }

      else {
        if (is.null(conditions))
          stop("Error. You must specify which conditions you wish to compare when the factor has two or more conditions.\n")
        if (length(conditions) != 2)
          stop("Error. The argument conditions must contain the 2 conditions you wish to compare.")
        
        l <- conditions %in% pData(input)[,i]
        # If they are defined, they will be TRUE
        if (l[1] == TRUE && l[2] == TRUE) {

          if (!is.null(assayData(input)$exprs)) {
            datos1 <- assayData(input)$exprs[,which(pData(input)[,i] == conditions[1])]
            datos2 <- assayData(input)$exprs[,which(pData(input)[,i] == conditions[2])]
          } else {
            datos1 <- assayData(input)$counts[,which(pData(input)[,i] == conditions[1])]
            datos2 <- assayData(input)$counts[,which(pData(input)[,i] == conditions[2])]
          }                    
          # Define the comparison string
          comparison <- paste(conditions[1],conditions[2], sep=" - ")
          
          condition_lev = TRUE
        }
      }        
    }
  }  
      
  if (condition_fac == FALSE)
    stop("The factor you have written does not correspond with any of the ones you have defined.")
  
  if (condition_lev == FALSE)
    stop("The conditions you have written don't exist in the factor specified.\n")
  
  # Correction to make it work when there are simulated samples
  if (replicates == "no")
    replicates = "technical"
#  if (description(input)@samples[[1]] == "no")
#    description(input)@samples[[1]] = "technical"
 
  n1 <- ncol(as.matrix(datos1))
  n2 <- ncol(as.matrix(datos2))
  
  if (norm == "n") {      # no normalization
    datos1 <- round(datos1, 100)
    datos2 <- round(datos2, 100)
  }


  if (is.null(k)) {
      m1 <- min(datos1[noceros(datos1, num = FALSE)], na.rm = TRUE)
      m2 <- min(datos2[noceros(datos2, num = FALSE)], na.rm = TRUE)
      mm <- min(m1, m2)
      k <- mm/2    
  } 

  
    # Total counts for each gene:
  suma1 <- rowSums(as.matrix(datos1))
  suma2 <- rowSums(as.matrix(datos2))
  

    # All genes
  todos <- rownames(as.matrix(datos1))

    # Genes with counts in any condition
  concounts <- names(which(suma1+suma2 > 0))

  long <- 1000
  g.sinL <- NULL
  
  if (!is.null(featureData(input)@data$Length)) {
    g.sinL <- names(which(is.na(featureData(input)@data$Length)))
    if (any(!is.na(featureData(input)@data$Length)) == TRUE) 
      long <- featureData(input)@data[concounts, "Length"]
  }
  
  
  if (replicates == "technical") {  ### technical replicates
    suma1 <- suma1[concounts]
    suma2 <- suma2[concounts]    

    #-------------------------------------------------------------------------#
    # Normalization of counts for each condition (aggregating replicates)

    if (norm == "rpkm") {      # RPKM
      suma1.norm <- rpkm(suma1, long = long, k = k, lc = lc)
      suma2.norm <- rpkm(suma2, long = long, k = k, lc = lc)
    }

    
    if (norm == "uqua") {
      suma.norm <- uqua(cbind(suma1, suma2), long = long, lc = lc, k = k)
      suma1.norm <- as.matrix(suma.norm[ ,1])
      suma2.norm <- as.matrix(suma.norm[ ,2])
    }

    
    if (norm == "tmm") {
      suma.norm <- tmm(as.matrix(cbind(suma1, suma2)), long = long,
                       lc = lc, k = k)      
      suma1.norm <- as.matrix(suma.norm[ ,1])
      suma2.norm <- as.matrix(suma.norm[ ,2])
    }
    
  }
   

    #-------------------------------------------------------------------------#

    ## Noise distribution

  if ((n1+n2)>2) {   # with real samples

    datitos <- cbind(datos1, datos2)
    datitos <- datitos[concounts,]

    gens.sin0 <- setdiff(concounts, g.sinL)

    if (norm == "n") {       # no normalization
      datitos.0 <- sinceros(datitos, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "rpkm") {      # RPKM
      datitos.0 <- rpkm(datitos, long = long, k = k, lc = lc)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "uqua") {      # Upper Quartile
      datitos.0 <- uqua(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "tmm") {
      datitos.0 <- tmm(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    datos1.norm <- datitos.norm[ ,1:n1]
    datos2.norm <- datitos.norm[ ,(n1+1):(n1+n2)]

    if (n1 > 1) {
      MD1 <- MD(dat = datos1.norm)
    } else { MD1 <- NULL }

    if (n2 > 1) {
      MD2 <- MD(dat = datos2.norm)
    } else { MD2 <- NULL }

    
  } else {                 # with simulated samples

    if (nss == 0) {
      nss <- 5
    }

    datos.sim <- sim.samples(counts1 = sinceros(suma1, k = k),
                             counts2 = sinceros(suma2, k = k),
                             pnr = pnr, nss = nss, v = v)

    nn <- sapply(datos.sim, ncol)

    dat.sim.norm <- vector("list", length = 2)

    datitos <- cbind(datos.sim[[1]], datos.sim[[2]])

    sumita <- rowSums(datitos)
    g.sin0 <- names(which(sumita > 0))
    gens.sin0 <- setdiff(g.sin0, g.sinL)

    if (norm == "n") {       # no normalization
      datitos.0 <- sinceros(datitos, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "rpkm") {      # RPKM
      datitos.0 <- rpkm(datitos, long = long, k = k, lc = lc)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "uqua") {      # Upper Quartile
      datitos.0 <- uqua(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    if (norm == "tmm") {
      datitos.0 <- tmm(datitos, long = long, lc = lc, k = k)
      datitos.norm <- datitos.0[gens.sin0, ]
    }

    dat.sim.norm[[1]] <- datitos.norm[ ,1:nn[1]]
    dat.sim.norm[[2]] <- datitos.norm[ ,(nn[1]+1):sum(nn)]

    MD1 <- MD(dat = dat.sim.norm[[1]])
    MD2 <- MD(dat = dat.sim.norm[[2]])

  }

  Mr <- c(as.numeric(MD1$M), as.numeric(MD2$M))
  Dr <- c(as.numeric(MD1$D), as.numeric(MD2$D))

  
  
    #-------------------------------------------------------------------------#

    ## M and D for different experimental conditions

  if (replicates == "technical" & norm != "n") {
    
    MDs <- MD(dat = cbind(suma1.norm, suma2.norm))
    
    lev1 <- suma1.norm[,1]
    lev1 <- lev1[todos]
        
    lev2 <- suma2.norm[,1]
    lev2 <- lev2[todos]

  } else {

    if (norm == "n" & (n1+n1) == 2) {
      datos1.norm <- sinceros(as.matrix(datos1)[concounts,], k = k)
      datos2.norm <- sinceros(as.matrix(datos2)[concounts,], k = k)
    }

    resum1.norm <- rowMeans(as.matrix(datos1.norm))
    resum2.norm <- rowMeans(as.matrix(datos2.norm))

    lev1 <- resum1.norm[todos]
    lev2 <- resum2.norm[todos]
    
    MDs <- MD(dat = cbind(resum1.norm, resum2.norm))

  }
  

      ## Completing M and D
  
  names(lev1) <- names(lev2) <- todos
  
  Ms <- as.numeric(MDs$M)
  names(Ms) <- rownames(MDs$M)
  Ms <- Ms[todos]
  names(Ms) <- todos

  Ds <- as.numeric(MDs$D)
  names(Ds) <- rownames(MDs$D)
  Ds <- Ds[todos]
  names(Ds) <- todos

  ## Results
  list("k" = k, "comp" = comparison, "Level1" = lev1, "Level2" = lev2, "Ms" = Ms, "Ds" = Ds, "Mn" = Mr, "Dn" = Dr)
 
}






#######################################################################
#######################################################################



allMDbio = function (input, factor, conditions, k = 0.5, norm = "rpkm", lc = 1,
                      r = 10, a0per = 0.9, nclust = 15, filter = 1, depth = NULL,
                     cv.cutoff = 0, cpm = 1)
  
  # input:  Set of data of type Input
  
  # conditions: Levels of the factor to be compared (when the factor has more than 2 levels)
  
  # k:      When counts = 0, 0 will be changed to k. By default, k = 0.5.
  
  # norm:   Normalization method. It can be one of "rpkm" (default), "uqua"
  #         (upper quartile), "tmm" (trimmed mean of M) or "n" (no normalization).

# lc:     Length correction in done by dividing expression by length^lc.
#         By default, lc = 1. 

# r:      Number of permutations to compute null distribution (r=10).
# a0per:  Percentile of S to compute a0. If NULL, a0 = 0. (a0per = 0.9)

{
  
  
  
  # Check if the factor introduced is already defined
  # If the factor introduced is defined and has more than 2 conditions, 
  # it will check if the conditions specified are defined too
  condition_fac = FALSE
  condition_lev = FALSE
  
  datos1 <- datos2 <- matrix()
  
  for (i in colnames(pData(input))) {
    if (factor == i) {
      condition_fac = TRUE
      if (length(levels(pData(input)[,i])) == 2) {
        
        if (!is.null(assayData(input)$exprs)) {
          datos1 <- assayData(input)$exprs[,which(pData(input)[,i] ==levels(pData(input)[,i])[1])]
          datos2 <- assayData(input)$exprs[,which(pData(input)[,i] ==levels(pData(input)[,i])[2])]
        } else {
          datos1 <- assayData(input)$counts[,which(pData(input)[,i] ==levels(pData(input)[,i])[1])]
          datos2 <- assayData(input)$counts[,which(pData(input)[,i] ==levels(pData(input)[,i])[2])]
        }
        
        # Define the comparison string
        comparison <- paste(levels(pData(input)[,i])[1], levels(pData(input)[,i])[2], sep=" - ")
        
        condition_lev = TRUE
        
        if (!((ncol(datos1) > 1) && (ncol(datos2) > 1))) 
          stop("Error. NOISeqBIO needs at least 2 biological replicates per condition.\n")
      }
      
      else {
        if (is.null(conditions))
          stop("Error. You must specify which conditions you wish to compare when the factor has two or more conditions.\n")
        if (length(conditions) != 2)
          stop("Error. The argument conditions must contain the 2 conditions you wish to compare.")
        
        l <- conditions %in% pData(input)[,i]
        # If they are defined, they will be TRUE
        if (l[1] == TRUE && l[2] == TRUE) {
          
          if (!is.null(assayData(input)$exprs)) {
            datos1 <- assayData(input)$exprs[,which(pData(input)[,i] == conditions[1])]
            datos2 <- assayData(input)$exprs[,which(pData(input)[,i] == conditions[2])]
          } else {
            datos1 <- assayData(input)$counts[,which(pData(input)[,i] == conditions[1])]
            datos2 <- assayData(input)$counts[,which(pData(input)[,i] == conditions[2])]
          }                    
          # Define the comparison string
          comparison <- paste(conditions[1],conditions[2], sep=" - ")
          
          condition_lev = TRUE
        }
      }        
    }
  }  
  
  if (condition_fac == FALSE)
    stop("The factor specified does not correspond with any of the ones you have defined.")
  
  if (condition_lev == FALSE)
    stop("The conditions specified don't exist for the factor specified.\n")
  
  
  ##-------------------------------------------------------------------------##
  
  
  # Number of observations within each condition 
  n1 <- ncol(as.matrix(datos1))
  n2 <- ncol(as.matrix(datos2))
  if (max(n1,n2) == 1) 
    stop("There is only one replicate per condition. Please, use NOISeq instead of NOISeqBIO.\n")
  
  # Rounding off data
  if (norm == "n") {      # no normalization
    datos1 <- round(datos1, 10)
    datos2 <- round(datos2, 10)
  }
  
  
  # Computing k
  if (is.null(k)) {
    m1 <- min(datos1[noceros(datos1, num = FALSE)], na.rm = TRUE)
    m2 <- min(datos2[noceros(datos2, num = FALSE)], na.rm = TRUE)
    k <- min(m1, m2)/2          
  } 
  
  
  # Total counts for each gene:
  suma1 <- rowSums(as.matrix(datos1))
  suma2 <- rowSums(as.matrix(datos2))
  
  # Genes with counts in any condition
  concounts <- names(which(suma1+suma2 > 0))
  
  # All genes
  todos <- rownames(as.matrix(datos1))
  
  
  # Gene length
  long <- 1000
  g.sinL <- NULL  # genes with no length defined
  
  if (!is.null(featureData(input)@data$Length)) {
    g.sinL <- names(which(is.na(featureData(input)@data$Length)))
    if (any(!is.na(featureData(input)@data$Length)) == TRUE) 
      long <- featureData(input)@data[concounts, "Length"]
  }
  
  
  # Genes with counts and with length
  gens.sin0 <- setdiff(concounts, g.sinL)
  
  
  # cond1 and cond2 in the same matrix
  datitos <- cbind(datos1, datos2)  
  datitos <- datitos[concounts,]  # selecting only genes with counts
  
  # Sequencing depth when filtering method = 3
  if (filter == 3 && is.null(depth)) depth = colSums(datitos)
  
  
  #-------------------------------------------------------------------------#
  #-------------------------------------------------------------------------#
  
  
  
  ## Normalization  
  
  if (norm == "n") {       # no normalization
    datitos.0 <- sinceros(datitos, k = k)
    datitos.norm <- datitos.0[gens.sin0, ]
  }
  
  if (norm == "rpkm") {      # RPKM
    datitos.0 <- rpkm(datitos, long = long, k = k, lc = lc)
    datitos.norm <- datitos.0[gens.sin0, ]
  }
  
  if (norm == "uqua") {      # Upper Quartile
    datitos.0 <- uqua(datitos, long = long, lc = lc, k = k)
    datitos.norm <- datitos.0[gens.sin0, ] 
  }
  
  if (norm == "tmm") {
    datitos.0 <- tmm(datitos, long = long, lc = lc, k = k)
    datitos.norm <- datitos.0[gens.sin0, ]       
  }
  
  
  
  #-------------------------------------------------------------------------#
  
  ## Filtering out low count features
  if (filter != 0) {
    datos.filt = filtered.data(dataset = datitos.norm, factor = c(rep("cond1", n1), rep("cond2", n2)), 
                               norm = TRUE, depth = depth, method = filter, cv.cutoff = cv.cutoff, cpm = cpm)    
  } else {
    datos.filt = datitos.norm
  }
  
  datos1.filt <- datos.filt[ ,1:n1]
  datos2.filt <- datos.filt[ ,(n1+1):(n1+n2)]   
  
  
  #-------------------------------------------------------------------------#
  
  
  
  ## Noise distribution
  
  Zr = NULL
  
  if (n1+n2 <= 8) {  # sharing information within clusters
    
    Zr = share.info(mydata = datos.filt, n1 = n1, n2 = n2, r = r, nclust = nclust)
    
  } else {   # r permutations
    
    for (i in 1:r) {
      
      print(paste("r =", i))
      
      mipermu = sample(1:(n1+n2))
            
      mipermu = datos.filt[,mipermu]
      
      mean1 = rowMeans(mipermu[,1:n1])
      mean2 = rowMeans(mipermu[,(n1+1):(n1+n2)])        
      
      sd1 = apply(mipermu[,1:n1], 1, sd)
      sd2 = apply(mipermu[,(n1+1):(n1+n2)], 1, sd)
      
      myparam = list("n" = c(n1,n2), "sd" = cbind(sd1,sd2))
      
      MDperm <- MDbio(dat = cbind(mean1, mean2), param = myparam, a0per = a0per)
      
      Zr = cbind(Zr, myDfunction(mydif = MDperm$D, myrat = MDperm$M, stat = 1, coef = 0.5))
       
    }
  }
  
  
  
  #-------------------------------------------------------------------------#
  
  
  ## Z-score for different experimental conditions (SIGNAL)
  
  mean1 = rowMeans(as.matrix(datos1.filt))
  mean2 = rowMeans(as.matrix(datos2.filt))        
  
  sd1 = apply(as.matrix(datos1.filt), 1, sd)
  sd2 = apply(as.matrix(datos2.filt), 1, sd)
  
  myparam = list("n" = c(n1,n2), "sd" = cbind(sd1,sd2))
  
  MDs <- MDbio(dat = cbind(mean1, mean2), param = myparam, a0per = a0per)
  
  
  Zs = myDfunction(mydif = MDs$D, myrat = MDs$M, stat = 1, coef = 0.5)
    
  
  
    
  
  #-------------------------------------------------------------------------#
  
  
  
  ## Completing M and D (in signal)
  
  lev1 <- mean1[todos]
  lev2 <- mean2[todos]
  names(lev1) <- names(lev2) <- todos
  
  Zs <- as.numeric(Zs)
  names(Zs) <- rownames(MDs$M)
  
  Zs <- Zs[todos]
  names(Zs) <- todos
  
  
  ## Computing Zn
  Zn = as.numeric(Zr)
  
  
  
  #-------------------------------------------------------------------------#
  
  ## Results
  list("k" = k, "comp" = comparison, "Level1" = lev1, "Level2" = lev2, "Zs" = Zs, "Zn" = Zn)
  
}




##############################################################################
##############################################################################




## Function to summarize difference and ratio information (D and D0)

myDfunction <- function (mydif, myrat, stat, coef) {
  
  if (stat == 1) { # linear combination of difference and ratio
    
    myDvalues = coef*mydif + (1-coef)*myrat
    
  }
  
  if (stat == 2) { # distance to origin from (ratio, difference)
    
    myDvalues = sign(mydif) * sqrt((mydif)^2 + (myrat)^2)  
    
  }
  
  myDvalues
  
}




#######################################################################







