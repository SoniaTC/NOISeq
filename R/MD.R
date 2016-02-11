MD <-
function (dat = dat, selec = c(1:nrow(dat))) {
  pares <- as.matrix(combn(ncol(dat), 2))

  if (NCOL(pares) > 30) {  
    sub30 <- sample(1:NCOL(pares), size = 30, replace = FALSE)
    pares <- pares[,sub30]
  }
  
  mm <- NULL
  dd <- NULL
  for (i in 1:ncol(pares)) {
    a <- dat[selec,pares[1,i]]
    b <- dat[selec,pares[2,i]]
    mm <- cbind(mm, log(a/b, 2))
    dd <- cbind(dd, abs(a-b))
  }
  list("M" = mm, "D" = dd)
}



###########################################################################################
###########################################################################################



MDbio = function (dat = dat, selec = c(1:nrow(dat)), param = NULL, a0per = 0.9) {
  
  pares <- as.matrix(combn(ncol(dat), 2))
  
  mm <- NULL
  dd <- NULL
  for (i in 1:ncol(pares)) {
    a <- dat[selec,pares[1,i]]
    b <- dat[selec,pares[2,i]]
    mm <- cbind(mm, log(a/b, 2))
    dd <- cbind(dd, (a-b))  
  }
  
  
  ## Correcting (M,D) 
  sd.M = sqrt(param$sd[,1]^2 / (dat[,1]^2 * log(2)^2 * param$n[1]) + 
                param$sd[,2]^2 / (dat[,2]^2 * log(2)^2 * param$n[2]))
  
  sd.D = sqrt(param$sd[,1]^2/param$n[1] + param$sd[,2]^2/param$n[2])
  
  if(is.null(a0per)) {
    a0.M = a0.D = 0
  } else {
    
    if (a0per == "B") {
      B = 100
      a0.M <- B*max(sd.M, na.rm = TRUE)
      a0.D <- B*max(sd.D, na.rm = TRUE)
      
    } else {
      a0per = as.numeric(a0per)
      a0.M <- quantile(sd.M, probs = a0per, na.rm = TRUE)
      a0.D <- quantile(sd.D, probs = a0per, na.rm = TRUE)
    }
  }
  
  mm <- mm / (a0.M + sd.M)
  dd <- dd / (a0.D + sd.D)
  
  
  # Results
  list("M" = mm, "D" = dd)
}

