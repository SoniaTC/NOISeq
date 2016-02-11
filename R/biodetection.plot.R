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
  
  numgenes = nrow(dat)
  
  if (is.null(factor)) {  # per sample  
    cat("Biotypes detection is to be computed for:\n")
    print(colnames(dat))
    biotablas = vector("list", length = NCOL(dat))
    names(biotablas) = colnames(dat)
    
  } else {  # per condition
    mifactor = pData(input)[,factor]
    niveles = levels(mifactor)
    cat("Biotypes detection is to be computed for:\n")
    print(niveles)
    biotablas = vector("list", length = length(niveles))
    names(biotablas) = niveles
    dat = sapply(niveles, function (k) rowSums(as.matrix(dat[,grep(k, mifactor)])))
  }
  
    
  infobio <- as.character(featureData(input)@data$Biotype)
  
  genome <- 100*table(infobio)/sum(table(infobio))
  ordre <- order(genome, decreasing = TRUE)
  
  for (i in 1:length(biotablas)) {
    
    detect <- dat[,i] > k
    
    perdet1 <- genome*table(infobio, detect)[names(genome),"TRUE"]/
      table(infobio)[names(genome)]
    
    perdet2 <- 100*table(infobio, detect)[names(genome),"TRUE"] /
      sum(table(infobio, detect)[,"TRUE"])    
    
    biotablas[[i]] <- as.matrix(rbind(perdet1[ordre], perdet2[ordre]))
    rownames(biotablas[[i]]) <- c("detectionVSgenome", "detectionVSsample")    
  }
  
  mybiotable = list("genome" = genome[ordre], "biotables" = biotablas, "genomesize" = numgenes)  
  mybiotable
}







#############################################################################################
#############################################################################################
#############################################################################################



biodetection.plot <- function(dat, samples = c(1,2), plottype = c("persample", "comparison"), 
                              toplot = "protein_coding", toreport = FALSE,...) {  
  
  mypar = par(no.readonly = TRUE)

  if (length(samples) > 2) {
    stop("ERROR: This function cannot generate plots for more than 2 samples.\n 
         Please, use it as many times as needed to generate the plots for all your samples.\n")
  }
  
  if (is.numeric(samples)) samples = names(dat$biotables)[samples]
    
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
    
    barplot(biotable1[c(1,3),], main = samples[1],
            xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
            beside = TRUE, col = c("grey", 2), las = 2,
            ylim = c(0, ymaxL), border = c("grey", 2))
    
    barplot(biotable1[c(2,4),], main = samples[1],
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
    
    if (plottype == "persample") {    ### A plot for each sample separately
      
      # Datos1
      barplot(biotable1[c(1,3),], main = samples[1],
              xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
              beside = TRUE, col = c("grey", 2), las = 2,
              ylim = c(0, ymaxL), border = c("grey", 2))
      
      barplot(biotable1[c(2,4),], main = samples[1],
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
      barplot(biotable2[c(1,3),], main = samples[2],
              xlab = NULL, ylab = "%features", axis.lty = 1, legend = FALSE,
              beside = TRUE, col = c("grey", 4), las = 2,
              ylim = c(0, ymaxL), border = c("grey", 4))
      
      barplot(biotable2[c(2,4),], main = samples[2],
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

    
    if (plottype == "comparison") {   ## A plot comparing two samples with regard to genome and for % in sample
      
      lefttable = rbind(100*dat$biotables[[samples[1]]][1,]/dat$genome, 100*dat$biotables[[samples[2]]][1,]/dat$genome)     
      
      righttable = rbind(dat$biotables[[samples[1]]][2,], dat$biotables[[samples[2]]][2,])
      
      if (length(toplot) > 1) {
        toplot = toplot[1]
        print("WARNING: More than one biotype was provided, the proportion test will only by applied to the first biotype.")
      }
      
      if ((toplot != 1) && (toplot != "global")) {         
        numgenes = dat$genomesize
        myx = round(righttable[,toplot]*numgenes/100, 0)
        mytest = prop.test(x = myx, n = rep(numgenes, 2), alternative = "two.sided")
        if (is.numeric(toplot)) toplot = colnames(righttable)[toplot]
      }
      
      asumar = colSums(righttable)
      asumar = which(asumar < 0.25)
      if (length(asumar) > 1) {
        righttable = cbind(righttable[,-asumar], rowSums(righttable[,asumar]))
        colnames(righttable)[ncol(righttable)] = "Others"
      }            
      
      # Detection in the genome
      bbb = barplot(lefttable, main = "Biotype detection over genome total",
                    xlab = NULL, ylab = "% detected features", axis.lty = 1, legend = FALSE, cex.names = 0.8,
                    beside = TRUE, col = c(2,4), las = 2, density = 80, border = c(2,4), ylim = c(0,100))    
      bbb = colSums(bbb)/2
      lines(bbb, dat$genome, pch = 20, type = "o", lwd = 2)
      
      
      # %detection in the sample        
      barplot(righttable, main = "Relative biotype abundance in sample",
              xlab = NULL, ylab = "Relative % biotypes", axis.lty = 1, legend = FALSE,
              beside = TRUE, col = c(2, 4), las = 2, border = c(2,4))
                  
      legend(x = "topright", bty = "n", horiz = FALSE, pch = c(15,15,20), lwd = c(NA,NA,1), legend = c(samples, "% in genome"),
             col = c(2,4,1))
      
      if ((toplot != 1) && (toplot != "global")) {        
        print(paste("Percentage of", toplot, "biotype in each sample:"))
        names(mytest$estimate) = samples
        print(round(mytest$estimate*100, 4))
        print(paste("Confidence interval at 95% for the difference of percentages:", samples[1], "-", samples[2]))
        print(round(mytest$conf.int[1:2]*100, 4))
        if (mytest$p.value < 0.05) {
          print(paste("The percentage of this biotype is significantly DIFFERENT for these two samples (p-value =",
                      signif(mytest$p.value, 4), ")."))
        } else {
          print(paste("The percentage of this biotype is NOT significantly different for these two samples (p-value =",
                      signif(mytest$p.value, 4), ")."))
        }       
      }      
    }
    
    
  }
  
  # Reset with the default values
  if (!toreport) par(mypar)
  
  
}
