setClass("Biodetection", representation(dat="list"))
setClass("CD", representation(dat="list"))
setClass("CountsBio", representation(dat="list"))
setClass("GCbias", representation(dat="list"))
setClass("lengthbias", representation(dat="list"))
setClass("Saturation", representation(dat="list"))
setClass("PCA", representation(dat="list"))

setGeneric("explo.plot", function(object, ...) standardGeneric("explo.plot"))

setMethod("explo.plot", "Biodetection", function(object, samples = c(1,2), plottype = c("persample", "comparison"), 
                                                 toplot = "protein_coding", ...) 
  biodetection.plot(object@dat, samples = samples, plottype = plottype, toplot = toplot, ...))

setMethod("explo.plot", "CD", function(object, samples = NULL, ...) cd.plot(object@dat, samples = samples, ...))


setMethod("explo.plot", "CountsBio", function(object, samples = c(1,2), toplot = "global", plottype = c("barplot", "boxplot"),...)
          countsbio.plot(object@dat, samples, toplot, plottype, ...))

setMethod("explo.plot", "GCbias", function(object, samples = NULL, toplot = "global", ...)
          GC.plot(object@dat, samples = samples, toplot = toplot, ...))

setMethod("explo.plot", "lengthbias", function(object, samples = NULL, toplot = "global", ...)
          length.plot(object@dat, samples = samples, toplot = toplot, ...))

setMethod("explo.plot", "Saturation",
          function(object, samples = NULL, toplot = 1, yleftlim = NULL, yrightlim = NULL, ...)
          saturation.plot(object@dat, samples = samples, toplot = toplot, yleftlim = yleftlim, yrightlim = yrightlim, ...))

setMethod("explo.plot", "PCA", function(object, samples = 1:2, plottype = "scores", factor = NULL)
  PCA.plot(object@dat, samples = samples, plottype = plottype, factor = factor))



# Show methods for exploration objects
setMethod("show", "Biodetection",
    function(object) {
      cat("\n Reference genome: \n==========\n")
      names(dimnames(object@dat$genome)) = NULL
      print(object@dat$genome)
      for (i in c(1:length(object@dat$biotables))) {
        cat("\n",names(object@dat$biotables)[i],"\n==========\n")
        print(object@dat$biotables[[i]])
      }
    })
           

setMethod("show", "CD",
    function(object) {      
      cat("\n Confidence intervals for median of M to compare each sample to reference:\n=======\n")
      print(object@dat$DiagnosticTest)
      cat("\n Reference sample is:\n=======\n")
      print(object@dat$refColumn)
    })

setMethod("show", "CountsBio",
          function(object) {            
              cat("\n Summary: \n============\n")              
              print(object@dat$summary[[1]])            
          })


setMethod("show","GCbias",
          function(object) {
            x <- object@dat$RegressionModels
            for (i in 1:length(x)) {
              print(names(x)[i])
              print(summary(x[[i]]))
            }            
          })


setMethod("show","lengthbias",
          function(object) {
            x <- object@dat$RegressionModels
            for (i in 1:length(x)) {
              print(names(x)[i])
              print(summary(x[[i]]))
            }            
          })

setMethod("show","Saturation",
          function(object) {
            x <- dat2save(object)
            cat("\n Number of detected features at each sequencing depth: \n============\n") 
            for (i in 1:length(x)) {
              print(names(x)[i])
              print(x[[i]])
            }            
          })


setMethod("show","PCA",
          function(object) {
            x <- object$result$var.exp
            x = round(x*100,4)
            colnames(x) = c("%Var", "Cum %Var")
            rownames(x) = paste("PC", 1:nrow(x))
            cat("\n Percentage of total variance explained by each component: \n============\n") 
            print(x)                        
          })



# Coercion methods for exploration objects

setGeneric("dat2save", function(object)  standardGeneric("dat2save"))

setMethod("dat2save","Biodetection", function(object) object@dat)

setMethod("dat2save","CD", function(object) object@dat)

setMethod("dat2save","CountsBio", function(object)  object@dat$summary)

setMethod("dat2save","GCbias", function(object) object@dat$data2plot)

setMethod("dat2save","lengthbias", function(object) object@dat$data2plot)

setMethod("dat2save","Saturation", function(object) {
  
  muestras = vector("list", length = length(object@dat$depth))
  names(muestras) = names(object@dat$depth)
  
  for (i in 1:length(muestras)) {
    
    muestras[[i]] = object@dat$depth[[i]]

    for (j in 1:length(object@dat$saturation)) {
      muestras[[i]] = cbind(muestras[[i]], object@dat$saturation[[j]][[i]])
    }

    colnames(muestras[[i]]) <- c("depth", names(object@dat$saturation))
  }

  muestras  
})


setMethod("dat2save","PCA", function(object) object@dat$result)


############################################################################
############################# OUTPUT OBJECT ################################
############################################################################

setClass("myInfo",representation(method="character", k="numeric", lc="numeric", factor="vector",
                  v="numeric",nss="numeric",pnr="numeric",comparison="vector",replicates="character"))
setClass("Output",representation(results="list"),contains="myInfo")

setValidity("Output",
    function(object) {
        if (!(is.character(object@method))) {
          return(paste("Method must be a string"))
        } else if (!(is.numeric(object@k))) {
          return(paste("k must be numeric"))
        } else if (!(is.numeric(object@lc))) {
          return(paste("lc must be numeric"))
        } else if (!(is.vector(object@factor))) {
          return(paste("Factor must be a vector of strings"))
        } else if (!(is.numeric(object@v))) {
          return(paste("v must be numeric"))
        } else if (!(is.numeric(object@nss))) {
          return(paste("nss must be numeric"))
        } else if (!(is.numeric(object@pnr))) {
          return(paste("pnr must be numeric"))
        } else if (!(is.vector(object@comparison))) {
          return(paste("Comparison must be a vector of strings"))
        } else if (!(is.list(object@results))) {
          return(paste("Results must be a list of data.frames"))
        } else {
          return(TRUE)
        }
    })

Output <-
function (data, method, k, lc, factor, v, nss, pnr, comparison, replicates) {
  
  new("Output",results=data, method = method, k = k, lc = lc, factor = factor, v = v, nss = nss, 
      pnr = pnr, comparison = comparison, replicates = replicates)
}

setMethod("show", "Output",
    function(object) {
      
      if (object@method == "n")
        object@method = "none"
      
      for (i in 1:length(object@results)) {
        cat("\nSummary",i,"\n=========\n")
        cat("\nYou are comparing",object@comparison[i],"from", object@factor[i], "\n\n")
        print(head(object@results[[i]][order(object@results[[i]][,5], decreasing = TRUE),]))
      }
      cat("\nNormalization\n")
      cat("\tmethod:", object@method, "\n")
      cat("\tk:", object@k, "\n")
      cat("\tlc:", object@lc, "\n")
      
      # Simulated samples
      if (object@replicates == "no") {
        cat("\nYou are working with simulated replicates:\n")
        cat("\tpnr:",object@pnr,"\n")
        cat("\tnss:",object@nss,"\n")
        cat("\tv:",object@v,"\n")
      } 
      # With biological or technical replicates
      else {
        cat("\nYou are working with",object@replicates, "replicates\n")
      }
           
})
