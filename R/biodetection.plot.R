biodetection.dat <- function(input, factor = NULL, k = 0) {

  if (inherits(input,"eSet") == FALSE)
    stop("Error. The input data must be an eSet object\n")

  if (any(!is.na(featureData(input)@data$Biotype)) == FALSE)
    stop ("No biological classification was provided.\nPlease run addData() function to add 
          this information\n")

  if (!is.null(assayData(input)$exprs)) {
    dat <- as.matrix(assayData(input)$exprs)
    mysamples = colnames(assayData(input)$exprs)
  } else {
    dat <- as.matrix(assayData(input)$counts)
    mysamples = colnames(assayData(input)$counts)
  }
  
  if (is.null(factor)) {  # per sample  
    print("Biotypes detection is to be computed for:")
    print(colnames(dat))
    biotablas = vector("list", length = NCOL(dat))
    names(biotablas) = colnames(dat)
    
  } else {  # per condition
    mifactor = pData(input)[,factor]
    niveles = levels(mifactor)
    print("Biotypes detection is to be computed for:")
    print(niveles)
    biotablas = vector("list", length = length(niveles))
    names(biotablas) = niveles
    dat = sapply(niveles, function (k) rowSums(dat[,grep(k, mifactor)]))
  }
  
    
  infobio <- as.character(featureData(input)@data$Biotype)
  
  genome <- 100*table(infobio)/sum(table(infobio))
  ordre <- order(genome, decreasing = TRUE)
  
  for (i in 1:length(biotablas)) {
    
    detect <- dat[,i] > k
    
    perdet1 <- genome*table(infobio, detect)[names(genome),2]/
      table(infobio)[names(genome)]
    
    perdet2 <- 100*table(infobio, detect)[names(genome),2] /
      sum(table(infobio, detect)[,2])    
    
    biotablas[[i]] <- as.matrix(rbind(perdet1[ordre], perdet2[ordre]))
    rownames(biotablas[[i]]) <- c("detectionVSgenome", "detectionVSsample")    
  }
  
  mybiotable = list("genome" = genome[ordre], "biotables" = biotablas)  
  mybiotable
}







#############################################################################################
#############################################################################################
#############################################################################################



biodetection.plot <- function(dat, samples = c(1,2), ...) {  

  if (length(samples) > 2) {
    stop("ERROR: This function cannot generate plots for more than 2 samples.\n 
         Please, use it as many times as needed to generate the plots for all your samples.\n")
  }
    
  biotable1 <- rbind(dat$genome, dat$biotables[[samples[1]]], rep(0, length(dat$genome)))
  
    
  
  # Computing ylim for left and right axis
     
  if (ncol(biotable1) >= 3) {
    ymaxL <- ceiling(max(biotable1[,1:3], na.rm = TRUE))
    ymaxR <- max(biotable1[,-c(1:3)], na.rm = TRUE)        
  } else {
    ymaxL <- ceiling(max(biotable1, na.rm = TRUE))
    ymaxR = 0 
  } 
  
  
  if (length(samples) == 2) {
    
    biotable2 <- rbind(dat$genome, dat$biotables[[samples[2]]], rep(0, length(dat$genome)))
    
    if (ncol(biotable2) >= 3) {
      ymax2 <- ceiling(max(biotable2[,1:3], na.rm = TRUE))
      ymax2sin <- max(biotable2[,-c(1:3)], na.rm = TRUE)
      ymaxR <- ceiling(max(ymaxR, ymax2sin))
    } else {
      ymax2 <- ceiling(max(biotable2, na.rm = TRUE))      
    }    
    
    ymaxL = max(ymaxL, ymax2)   
  }
  
  
  
  # Rescaling biotables (datos2)
  if (length(samples) == 2) {        
    if (ncol(biotable2) >= 3) biotable2[,-c(1:3)] <- biotable2[,-c(1:3)]*ymaxL/ymaxR    
  }  
  
  # Rescaling biotables (datos1)  
  if (ncol(biotable1) >= 3) biotable1[,-c(1:3)] <- biotable1[,-c(1:3)]*ymaxL/ymaxR
  
  
  ## PLOTS
  
  if (length(samples) == 1) {   # Plot (1 sample) - 2 scales    
    
    par(mar = c(11, 4, 2, 2))
    
    barplot(biotable1[c(1,3),], main = names(dat$biotables)[samples[1]],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))
    
    barplot(biotable1[c(2,4),], main = names(dat$biotables)[samples[1]],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(2, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 2, add = TRUE)    
    
    if (ymaxR > 0) {  # if number of biotypes >= 3 so we have left and right axis
      axis(side=4, at = pretty(c(0,ymaxL), n = 5), 
           labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))    
      abline(v = 9.5, col = 3, lwd = 2, lty = 2)        
    }    
    
    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 2, 2), density = c(NA,30,NA),
           border = c("grey", 2, 2),
           legend = c("% in genome", "detected", "% in sample"))
    
    
  } else {   # Plot (2 samples)    
    
    par(mar = c(11, 4, 2, 2))    
    
    # Datos1
    barplot(biotable1[c(1,3),], main = names(dat$biotables)[samples[1]],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))
    
    barplot(biotable1[c(2,4),], main = names(dat$biotables)[samples[1]],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(2, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 2, add = TRUE)
    
    if (ymaxR > 0) {   # if number of biotypes >= 3 so we have left and right axis
      axis(side=4, at = pretty(c(0,ymaxL), n = 5), 
           labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))      
      abline(v = 9.5, col = 3, lwd = 2, lty = 2)
    }  
    
    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 2, 2), density = c(NA,30,NA), border = c("grey", 2, 2),
           legend = c("% in genome", "detected", "% in sample"))
    
    
    # Datos2        
    barplot(biotable2[c(1,3),], main = names(dat$biotables)[samples[2]],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 4), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 4))
    
    barplot(biotable2[c(2,4),], main = names(dat$biotables)[samples[2]],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(4, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 4, add = TRUE)
    
    if (ymaxR > 0) {   # if number of biotypes >= 3 so we have left and right axis
      axis(side=4, at = pretty(c(0,ymaxL), n = 5), 
           labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))
      abline(v = 9.5, col = 3, lwd = 2, lty = 2)
    }      
    
    legend(x = "topright", bty = "n", horiz = FALSE, 
           fill = c("grey", 4, 4), density = c(NA,30,NA),
           border = c("grey", 4, 4),
           legend = c("% in genome", "detected", "% in sample"))
  }
  
  # Reset with the default values
  par(mar = c(5, 4, 4, 4) + 0.1)
  
  
}
