require(methods)

setClass("Biodetection", representation(dat="list"))
setClass("CD", representation(dat="data.frame"))
setClass("CountsBio", representation(dat="list"))
setClass("DLBio", representation(dat="list"))
setClass("Saturation", representation(dat="list"))

setGeneric("explo.plot", function(object, ...) standardGeneric("explo.plot"))


setMethod("explo.plot", "Biodetection", function(object) biodetection.plot(object@dat))
setMethod("explo.plot", "CD", function(object) cd.plot(object@dat))
setMethod("explo.plot", "CountsBio", function(object, toplot = 1, samples = NULL, ylim = NULL)
          countsbio.plot(object@dat, toplot = toplot, samples = samples, ylim = ylim))
setMethod("explo.plot", "DLBio", function(object, samples = NULL, toplot = "protein_coding", ylim = NULL)
          DLbio.plot(object@dat, samples = samples, toplot = toplot,ylim = ylim))
setMethod("explo.plot", "Saturation",
          function(object, toplot = 1, samples = NULL, ylim = NULL, yrightlim = NULL)
          saturation.plot(object@dat, toplot = toplot, samples = samples, ylim = ylim, yrightlim = yrightlim))

# Show methods for exploration objects
setMethod("show", "Biodetection",
    function(object) {
      for (i in c(1:length(object@dat$samples))) {
        cat("\n",object@dat$samples[i],"\n==========\n")
        print(object@dat[i])
      }
    })
           

setMethod("show", "CD",
    function(object) {
      
      cat("\nSummary\n=======\n")
      print(head(object@dat))
      
    })

setMethod("show", "CountsBio",
          function(object) {
            for (i in c(1:length(object@dat$quart))) {
              cat("\n",names(object@dat$quart)[i],"\n============\n")
              aux <- cbind(object@dat$quart[[i]],object@dat$bionum[rownames(object@dat$quart[[i]])])
              colnames(aux)[4] <- "bionum"
              print(head(aux))
            }
            
          })

setMethod("show","DLBio",
          function(object) {
            x <- dat2save(object)
            print(head(x[[1]]))
          })

setMethod("show","Saturation",
          function(object) {
            x <- dat2save(object)
            print(head(x[[1]]))
          })

# Coercion methods for exploration objects

setGeneric("dat2save", function(object)  standardGeneric("dat2save"))

setMethod("dat2save","Biodetection", function(object) {
  aux <- list()

  if (length(object@dat$samples) == 2)
    aux <- list(object@dat$table, object@dat$table2)
  else
    aux <- list(object@dat$table)
  
  names(aux) <- object@dat$samples

  aux  
})

setMethod("dat2save","CD", function(object) object@dat)

setMethod("dat2save","CountsBio", function(object)  {

  x <- list()
  for (i in c(1:length(object@dat$quart))) {
    aux <- cbind(object@dat$quart[[i]],object@dat$bionum[rownames(object@dat$quart[[i]])])
    colnames(aux)[4] <- "bionum"
    x[[i]] <- aux
  }

  names(x) <- names(object@dat$quart)

  x

  })

setMethod("dat2save","DLBio", function(object) {

  biotipos <- names(object@dat$result)
  muestras <- names(object@dat$depth)

  mat <- matrix(0, nrow = length(biotipos), ncol = 7, dimnames=list(biotipos, c(1:7)))
  lista <- list()

  for (i in 1:length(muestras)) {
    mat.aux <- mat

    for (j in 1:length(biotipos)) {
      mat.aux[j,3:7] <- object@dat$result[[j]][[i]]      
    }

    mat.aux[,1] <- object@dat$bionum
    mat.aux[,2] <- object@dat$biolength

    colnames(mat.aux) <- c("bionum","biolength",
                           paste("depth_",sprintf("%0.3f", object@dat$depth[[i]]/10^6), sep=""))
    
    lista[[i]] <- mat.aux
  }

  names(lista) <- muestras
  lista
  
})

setMethod("dat2save","Saturation", function(object) {

  biotipos <- names(object@dat$saturation)
  muestras <- names(object@dat$depth)

  mat <- matrix(0, nrow = length(biotipos), ncol = 11, dimnames=list(biotipos, c(1:11)))
  lista <- list()

  for ( i in 1:length(muestras)) {
    mat.aux <- mat

    for (j in 1:length(biotipos)) {
      mat.aux[j,2:6]  <- object@dat$saturation[[j]][[i]]
      mat.aux[j,7:11] <- object@dat$newdet[[j]][[i]]
    }

    mat.aux[,1] <- object@dat$bionum

    colnames(mat.aux) <- c("bionum",
                           paste("num", sprintf("%0.3f", object@dat$depth[[i]]/10^6),sep="_"),
                           paste("new", sprintf("%0.3f", object@dat$depth[[i]]/10^6),sep="_"))
    
    lista[[i]] <- mat.aux
  }

  names(lista) <- muestras
  lista  
  
})

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
