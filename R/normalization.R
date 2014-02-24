#######################

tmm = function (datos, long = 1000, lc = 0, k = 0, refColumn = 1,
                logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                Acutoff = -1e10) {
  
  # lc: Length correction. Expression is divided by long^lc. lc can be any real number.
  
  L <- (long/1000)^lc
  datos = datos/L
  
  total <- colSums(as.matrix(datos))
  
  datos0 <- sinceros(datos, k)
  
  if (ncol(as.matrix(datos)) > 1) {
    
    fk <- .calcNormFactors(as.matrix(datos), refColumn = refColumn,
                           logratioTrim = logratioTrim, sumTrim = sumTrim,
                           doWeighting = doWeighting, Acutoff = Acutoff)
          
    fk = fk * total
    
    datos.norm <- t(t(datos0)/fk)*10^6
    
  } else {
    
    datos.norm <- datos0/L
    
  }
  
  na.omit(datos.norm)  
  
}


#######################


rpkm <- function (datos, long = 1000, lc = 1, k = 0) {

  total <- colSums(as.matrix(datos))

  datos0 <- sinceros(datos, k)                   

  datos.norm <- (t(t(datos0)/total)*10^6)/((long/1000)^lc)
  
  na.omit(datos.norm)   
}


##################################


uqua <- function (datos, long = 1000, lc = 0, k = 0) {
  
# lc: Length correction. Expression is divided by long^lc. lc can be any real number.

  L <- (long/1000)^lc
  datos = datos/L

  datos0 <- sinceros(datos, k)

  if (ncol(as.matrix(datos)) > 1) {
    sumatot <- rowSums(datos)
    supertot <- sum(sumatot)
    counts0 <- which(sumatot == 0)

    if (length(counts0) > 0) {
      datitos <- datos[-counts0,]
    } else {
      datitos <- datos
    }

    q3 <- apply(datitos, 2, quantile, probs = 0.75)
    d <- q3*supertot/sum(q3)

    datos.norm <- t(t(datos0)/d)*10^6

  } else {

    datos.norm <- datos0/L

  }
  
  na.omit(datos.norm)  
  
}



##################################

## Taken from the edgeR package with minor modifications

.calcNormFactors <- function(object, method=c("TMM","quantile"), refColumn=NULL,
                            logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, 
                            quantile=0.75) {
                            
  method <- match.arg(method)
                            
  if( is.matrix(object) ) {
    if(is.null(refColumn))
      refColumn <- 1
    data <- object
    libsize <- colSums(data)
  } else {
    stop("calcNormFactors() only operates on 'matrix' objects")
  }

  f <- switch(method,
              TMM = apply(data,2,.calcFactorWeighted,ref=data[,refColumn], 
                          logratioTrim=logratioTrim, sumTrim=sumTrim, doWeighting=doWeighting, 
                          Acutoff=Acutoff),              
              quantile = .calcFactorQuantile(data, libsize, q=quantile))
              
  f <- f/exp(mean(log(f)))

  return(f)
  
}

.calcFactorQuantile <- function (data, lib.size, q=0.75) {
    y <- t(t(data)/lib.size)
    f <- apply(y,2,function(x) quantile(x,p=q))
    f/exp(mean(log(f)))
}


.calcFactorWeighted <- function(obs, ref, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10) {

  if( all(obs==ref) )
    return(1)

  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  nO <- sum(obs)
  nR <- sum(ref)
  logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
  absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
  
  # remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  
  # taken from the original mean() function
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  
  #keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  # a fix from leonardo ivan almonacid cardenas, since rank() can return
  # non-integer values when there are a lot of ties
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)
  
  if (doWeighting) 
    2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
  else
    2^( mean(logR[keep], na.rm=TRUE) )
}
