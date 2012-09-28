biodetection.dat <- function(input, data_selection=c(1), k = 0) {

  if (inherits(input,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")

  if (any(!is.na(featureData(input)@data$Biotype)) == FALSE)
    stop ("No biological classification provided.\nPlease run addData() function to add 
          this information\n")

  if (!is.null(assayData(input)$exprs)) {
    dat1 <- rowSums(as.matrix(assayData(input)$exprs[,data_selection[1]]))
    mysamples = colnames(assayData(input)$exprs)[data_selection]
  } else {
    dat1 <- rowSums(as.matrix(assayData(input)$counts[,data_selection[1]]))
    mysamples = colnames(assayData(input)$counts)[data_selection]
  }
    
  
  infobio <- featureData(input)@data$Biotype
  
  detect1 <- dat1 > k
  
  genome <- 100*table(infobio)/sum(table(infobio))
  
  ordre <- order(genome, decreasing = TRUE)
  
  perdet1 <- genome*table(infobio, detect1)[names(genome),2]/
    table(infobio)[names(genome)]
  
  perdet2 <- 100*table(infobio, detect1)[names(genome),2] /
    sum(table(infobio, detect1)[,2])
  
  ceros <- rep(0, length(genome))
  
  biotable <- as.matrix(rbind(genome[ordre], perdet1[ordre],
                              perdet2[ordre], ceros))
  rownames(biotable) <- c("genome", "detectionVSgenome", "detectionVSsample",
                          "ceros")
  
  #ymax1 <- round(max(biotable[,1:3])*1.2,0)
  ymax1 <- ceiling(max(biotable[,1:3], na.rm = TRUE))
  ymax1sin <- max(biotable[,-c(1:3)], na.rm = TRUE)
    
  if (length(data_selection) == 2) {

     if (!is.null(assayData(input)$exprs))
       dat2 <- rowSums(as.matrix(assayData(input)$exprs[,data_selection[2]]))
     else
       dat2 <- rowSums(as.matrix(assayData(input)$counts[,data_selection[2]]))
      
    detect2 <- dat2 > k
    
    perdet3 <- genome*table(infobio, detect2)[names(genome),2]/
      table(infobio)[names(genome)]
    
    perdet4 <- 100*table(infobio, detect2)[names(genome),2] /
      sum(table(infobio, detect2)[,2])
    
    biotable2 <- as.matrix(rbind(genome[ordre], perdet3[ordre],
                                 perdet4[ordre], ceros))
    rownames(biotable2) <- c("genome", "detectionVSgenome", "detectionVSsample",
                             "ceros")
    
    
    #ymax3 <- round(max(biotable2[,1:3])*1.2,0)
    ymax2 <- ceiling(max(biotable2[,1:3], na.rm = TRUE))
    ymax2sin <- max(biotable2[,-c(1:3)], na.rm = TRUE)
    
    ymaxL <- max(ymax1, ymax2)
    ymaxR <- ceiling(max(ymax1sin, ymax2sin))
    
    # scaling data on the right (datos2)
    biotable2b <- biotable2
    biotable2b[,-c(1:3)] <- biotable2b[,-c(1:3)]*ymaxL/ymaxR
    
  } else {
    
    ymaxL <- ymax1
    ymaxR <- ceiling(ymax1sin)
    
  }
  
  # scaling data on the right (datos1)
  biotable1 <- biotable
  biotable1[,-c(1:3)] <- biotable1[,-c(1:3)]*ymaxL/ymaxR
  
  if (length(data_selection) == 2) {
    
    mybiotable <- list("table"=biotable1, "table2"=biotable2b, 
                       "params" = c(ymaxR, ymaxL),
                       "samples"=mysamples)
    
  } else {
    
    mybiotable <- list("table"=biotable1, "params" = c(ymaxR, ymaxL),
                       "samples"=mysamples)
        
  }
  
  mybiotable
}




#############################################################################################
#############################################################################################
#############################################################################################




biodetection.plot <- function(dat) {  

  if (length (dat$samples) > 2) {
    stop("ERROR! This method cannot plot more than 2 different samples\n")
  }

  biotable1 <- dat$table
  ymaxR <- dat$params[1]
  ymaxL <- dat$params[2]

  if (length(dat$samples) == 1) {
    
    # Plot (1 sample) - 2 scales
    par(mar = c(10, 4, 2, 2))
    
    barplot(biotable1[c(1,3),], main = dat$samples,
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))
    
    barplot(biotable1[c(2,4),], main = dat$samples,
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(2, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 2, add = TRUE)
    
    axis(side=4, at = pretty(c(0,ymaxL), n = 5), 
         labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))
    
    abline(v = 9.5, col = 3, lwd = 2, lty = 2)
    
    text(x = 9, y = ymaxL*1.1, "Left axis", col = 3, font = 3, adj = 1)
    text(x = 10, y = ymaxL*1.1, "Right axis", col = 3, font = 3, adj = 0)
    
    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 2, 2), density = c(NA,30,NA),
           border = c("grey", 2, 2),
           legend = c("% in genome", "detected", "% in sample"))
    
    #dev.off()
    
  } else {
    
    # Plot (2 samples)
    
    par(mar = c(10, 4, 2, 2), mfrow = c(1,2))
    
    
    # Datos1
    barplot(biotable1[c(1,3),], main = dat$samples[1],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))
    
    barplot(biotable1[c(2,4),], main = dat$samples[1],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(2, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 2, add = TRUE)
    
    axis(side=4, at = pretty(c(0,ymaxL), n = 5), 
         labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))
    
    abline(v = 9.5, col = 3, lwd = 2, lty = 2)
    
    text(x = 9, y = ymaxL*1.1, "Left axis", col = 3, font = 3, adj = 1)
    text(x = 10, y = ymaxL*1.1, "Right axis", col = 3, font = 3, adj = 0)
    
    legend(x = "topright", bty = "n", horiz = FALSE,
           fill = c("grey", 2, 2), density = c(NA,30,NA), border = c("grey", 2, 2),
           legend = c("% in genome", "detected", "% in sample"))
    
    
    # Datos2
    biotable2b = dat$table2
      
    barplot(biotable2b[c(1,3),], main = dat$samples[2],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 4), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 4))
    
    barplot(biotable2b[c(2,4),], main = dat$samples[2],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c(4, 1), las = 2, density = 30,
            ylim = c(0, ymaxL), border = 4, add = TRUE)
    
    axis(side=4, at = pretty(c(0,ymaxL), n = 5), 
         labels = round(pretty(c(0,ymaxL), n = 5)*ymaxR/ymaxL, 1))
    
    text(x = 9, y = ymaxL*1.1, "Left axis", col = 3, font = 3, adj = 1)
    text(x = 10, y = ymaxL*1.1, "Right axis", col = 3, font = 3, adj = 0)
    
    abline(v = 9.5, col = 3, lwd = 2, lty = 2)
    
    legend(x = "topright", bty = "n", horiz = FALSE, 
           fill = c("grey", 4, 4), density = c(NA,30,NA),
           border = c("grey", 4, 4),
           legend = c("% in genome", "detected", "% in sample"))
  }
  
  # Reset with the default values
  par(mar = c(5, 4, 4, 4) + 0.1)
  layout(1)
  
}
