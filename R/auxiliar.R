#################################################################

busca <-
function (x, S) {
  which(S[,1] == x[1] & S[,2] == x[2])
}


#################################################################



int.mult <- function(lista, todos = NULL) {
  
  if(is.null(todos)) {
    todos <- unlist(lista)
  }

  comunes <- todos

  for(i in 1:length(lista)) {
    comunes <- intersect(comunes, lista[[i]])
  }

  comunes
}



#################################################################


n.menor <-
function (x, S1, S2) {
 
  length(which(S1 <= x[1] &  S2 <= x[2]))

}



#################################################################



noceros <-
function (x, num = TRUE, k = 0) {
  
  nn <- length(which(x > k))
  
  if (num) {
    nn
    
  } else {
    if(nn > 0) { which(x > k) } else { NULL }
  }
}



#################################################################


sinceros <-
function (datos, k) {
  datos = as.matrix(datos)
  datos0 <- as.matrix(datos)

  if (is.null(k)) {

    mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])

    kc <- mini0/2

    datos0[datos0 == 0] <- kc

  } else {

    datos0[datos0 == 0] <- k

  }
  
  datos0
}



#################################################################


#### Simulating samples


sim.samples <-
  function(counts1, counts2 = NULL, pnr = 1, nss = 5, v = 0.02) {
    seqdep <- c(sum(counts1), sum(counts2))
    num.reads1 <- (pnr + c(-v,v))*seqdep[1]
    
    muestras <- vector("list")
    
    muestras$c1 <- NULL
    for (s in 1:nss) {
      tama <- round(runif(1, num.reads1[1], num.reads1[2]), 0)
      muestras$c1 <- cbind(muestras$c1,
                           rmultinom(1, size = tama, prob = counts1))
    }
    
    if(!is.null(counts2)) {
      num.reads2 <- (pnr + c(-v,v))*seqdep[2]
      muestras$c2 <- NULL
      for (s in 1:nss) {
        tama <- round(runif(1, num.reads2[1], num.reads2[2]), 0)
        muestras$c2 <- cbind(muestras$c2,
                             rmultinom(1, size = tama, prob = counts2))
      }
    }
    
    muestras
  }



#################################################################




ranking <-
  function(results) {
    
    M <- results$M
    
    D <- results$D
    
    prob <- results$prob
    
    ## Changing NA by 0
    M[is.na(M)] <- 0
    D[is.na(D)] <- 0
    prob[is.na(prob)] <- 0
    
    
    ## Ranking
    
    #   ranking1 <- M*prob
    # 
    #   ranking2 <- sign(M)*prob
    # 
    #   ranking3 <- M*D
    # 
    #   ranking4 <- M*D*prob
    
    ranking5 <- sqrt(M*M + D*D)*sign(M)
    
    
    ## Ranking results
    #list(ranking1, ranking2, ranking3, ranking4, ranking5)
    theranking <- data.frame(rownames(results), ranking5)
    rownames(theranking) <- NULL
    colnames(theranking) <- c("ID", "statistic")
    
    theranking
    
  }




#################################################################



#############################################################################
##############   Plot with 2 different Y axis (left and right)   ############
#############################################################################


# By Ajay Shah (taken from [R] Plot 2 time series with different y axes (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html) 

# Modified by: Sonia Tarazona

### PARAMETERS (default):
# x: data to be drawn on X-axis
# yright: data to be drawn on Y right axis
# yleft: data to be drawn on Y left axis
# yrightlim (range(yright, na.rm = TRUE)): ylim for rigth Y-axis
# yleftlim (range(yleft, na.rm = TRUE)): ylim for left Y-axis
# xlab (NULL): Label for X-axis
# yylab (c("","")): Labels for right and left Y-axis
# pch (c(1,2)): Type of symbol for rigth and left data
# col (c(1,2)): Color for rigth and left data
# linky (TRUE): If TRUE, points are connected by lines.
# smooth (0): Friedman's super smoothing
# lwds (1): Line width for smoothed line
# length (10): Number of tick-marks to be drawn on axis
# ...: Other graphical parameters to be added by user (such as main, font, etc.)
###



plot.y2 <- function(x, yright, yleft, yrightlim = range(yright, na.rm = TRUE),
                    yleftlim = range(yleft, na.rm = TRUE),
                    xlim = range(x, na.rm = TRUE),
                    xlab = NULL, yylab = c("",""), lwd = c(2,2),
                    pch = c(1,2), col = c(1,2), type = c("o","o"),
                    linky = TRUE, smooth = 0, bg = c("white","white"),
                    lwds = 1, length = 10, ...,
                    x2 = NULL, yright2 = NULL, yleft2 = NULL, col2 = c(3,4))
{
  #par(mar = c(5,2,4,2), oma = c(0,3,0,3))

  ## Plotting RIGHT axis data

  plot(x, yright, axes = FALSE, ylab = "", xlab = xlab, ylim = yrightlim,
       xlim = xlim, pch = pch[1], type = type[1], lwd = lwd[1],
       col = col[1], ...)
  
  axis(4, pretty(yrightlim, length), col = 1, col.axis = 1)

  if (is.null(yright2) == FALSE) {
    points(x2, yright2, type = type[1], pch = pch[1], lwd = lwd[1], col = col2[1], ...)
  }
  
  #if (linky) lines(x, yright, col = col[1], ...)
  
  if (smooth != 0) lines(supsmu(x, yright, span = smooth), col = col[1], lwd = lwds, ...)
  
  if(yylab[1]=="") {
    mtext(deparse(substitute(yright)), side = 4, outer = FALSE, line = 2,
          col = 1,...)
  } else {
    mtext(yylab[1], side = 4, outer = FALSE, line = 2, col = 1, ...)
  }
  

  par(new = T)

  ## Plotting LEFT axis data
  
  plot(x, yleft, axes = FALSE, ylab = "" , xlab = xlab, ylim = yleftlim, 
       xlim = xlim, bg = bg[1],
       pch = pch[2], type = type[2], lwd = lwd[2], col = col[2], ...)
  
  box()
  
  axis(2, pretty(yleftlim, length), col = 1, col.axis = 1)

  if (is.null(yleft2) == FALSE) {
    points(x2, yleft2, type = type[2], pch = pch[2], bg = bg[2],
           lwd = lwd[2], col = col2[2], ...)
  }
  

  #if (linky) lines(x, yleft, col = col[2], ...)
  
  if (smooth != 0) lines(supsmu(x, yleft, span = smooth), col = col[2], lwd=lwds, ...)
  
  if(yylab[2] == "") {
    mtext(deparse(substitute(yleft)), side = 2, outer = FALSE, line = 2, col = 1, ...)
  } else {
    mtext(yylab[2], side = 2, outer = FALSE, line = 2, col = 1, ...)
  }
  
  
  ## X-axis
  axis(1, at = pretty(xlim, length))
  
   
}








#################################################################


## Log-scale for plots


logscaling = function (data, base = 2, k = 1) {
  
  # IDEA
  # plot(data,...,yaxt = "n")
  # axis(side = 2, at = donde, labels = etiquetas)
  
  logmaximo = round(max(data, na.rm = TRUE), 0)
  numceros = nchar(logmaximo)-1
    
  etiquetas = c(0, 10^(1:numceros))
  
  donde = log(etiquetas + k, base = base)
  
  data = log(data + k, base = base)
  
  list("data" = data, "at" = donde, "labels" = etiquetas)
    
}





##***************************************************************************##
##***************************************************************************##




miscolores <- colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)]
