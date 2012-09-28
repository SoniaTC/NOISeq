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

