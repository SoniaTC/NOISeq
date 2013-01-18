readData <-
function (data = NULL, factors = NULL, length = NULL, biotype = NULL, chromosome = NULL) {

  if (is.null(data))
    stop("Expression information must be provided to the readData function")

  if (is.null(factors))
    stop("Condition information must be provided to the readData funcion")

  if (is.null(length) == FALSE && is.vector(length) == FALSE)
    stop( "The length info should be a vector.")

  if (is.null(chromosome) == FALSE && ncol(chromosome) != 3)
    stop( "The chromosome object should be a matrix or data.frame with 3 columns: chromosome, start position and end position.")

  countData <- as.matrix( data )

  rowNames <- rownames(countData)
  rownames(factors) <- colnames(countData)
     
  pheno <- AnnotatedDataFrame(data=as.data.frame(factors))

  input <- ExpressionSet(
               assayData = countData,
               phenoData = pheno)

  if (!is.null(length))
    input <- addData(data = input, length = length)

  if (!is.null(biotype))
    input <- addData(data = input, biotype = biotype)

  if (!is.null(chromosome))
    input <- addData(data = input, chromosome = chromosome)

  input    

}


######################################################
######################################################
######################################################





addData <- function(data, length = NULL, biotype = NULL, chromosome = NULL, factors = NULL) {

  if (inherits(data,"eSet") == FALSE)
    stop("Error. You must give an eSet object.")
  
  if (is.null(length) == FALSE && is.vector(length) == FALSE)
    stop( "The length info should be a vector.")

  if (is.null(chromosome) == FALSE && ncol(chromosome) != 3)
    stop( "The chromosome object should be a matrix or data.frame with 3 columns: chromosome, start position and end position.")

  if (!is.null(assayData(data)$exprs))
    rowNames <- rownames(assayData(data)$exprs)
  else
    rowNames <- rownames(assayData(data)$counts)

  # If exists length
  if (!is.null(length)) {
    Length <- rep(NA,length(rowNames))
    names(Length) <- rowNames
    Length[rowNames] <- as.numeric(length[rowNames])
    
    featureData(data)@data <- cbind(featureData(data)@data, Length)
  }

  # If exists biotype
  if (!is.null(biotype)) {
    Biotype <- rep(NA,length(rowNames))
     names(Biotype) <- rowNames
    Biotype[rowNames] <- as.character(biotype[rowNames])
    
    featureData(data)@data <- cbind(featureData(data)@data, Biotype)
  }

  # If exists biotype
  if (!is.null(chromosome)) {

    Chromosome <- GeneStart <- GeneEnd <- rep(NA,length(rowNames))
    names(Chromosome) <- names(GeneStart) <- names(GeneEnd) <- rowNames
    
    Chromosome[rowNames] <- chromosome[rowNames,1]
    GeneStart[rowNames] <- as.numeric(chromosome[rowNames,2])
    GeneEnd[rowNames] <- as.numeric(chromosome[rowNames,3])
        
    featureData(data)@data <- cbind(featureData(data)@data, Chromosome, GeneStart, GeneEnd)
  }

  # If exists new factors
  if (!is.null(factors))
    phenoData(data)@data <- cbind(phenoData(data)@data, factors)
  
  data
  
}
