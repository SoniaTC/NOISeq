degenes <-
function (object, q = 0.9, M = NULL) {
  # object = noiseq output object
  # M = "up" (up-regulated in condition 1), "down" (down-regulated in condition 1), NULL (all differentially expressed genes)
  # q = probability threshold (between 0 and 1)

  if (class(object) != "Output")
    stop("You must give the object returned by the noiseq function\n")

  x <- object@results[[1]]
  
  noiseqbio = "theta" %in% colnames(x)[1:4]
  
  if (noiseqbio) {
    y <- na.omit(x[c("theta","prob")])
    colnames(y)[1] = "M"
  } else {
    y <- na.omit(x[c("M","D","prob")])
  }
  

  if (is.null(M)) {

    losdeg <- y[y[,"prob"] > q,]
    print(paste(dim(losdeg)[1], "differentially expressed features"))

  } else if (M == "up") {

    estos <- y[y[,"M"] > 0,]
    losdeg <- estos[estos[,"prob"] > q,]
    print(paste(dim(losdeg)[1], "differentially expressed features (up in first condition)"))

  } else if (M == "down") {

    estos <- y[y[,"M"] < 0,]
    losdeg <- estos[estos[,"prob"] > q,]
    print(paste(dim(losdeg)[1], "differentially expressed features (down in first condition)"))

  } else {

    stop("ERROR! Value for parameter M is not valid. Please, choose among NULL, 'up' or 'down'")

  }

  # Restore the object with the same "results" structure
  losdeg = x[rownames(losdeg),]

  losdeg[order(losdeg[,"prob"], decreasing = TRUE),]

}

