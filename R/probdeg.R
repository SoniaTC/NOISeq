probdeg <-
function (Mg, Dg, Mn, Dn, prec = 2) {

  # Mg, Dg -> signal
  # Mn, Dn -> noise
  # prec = precission (number of digits to round M and D)

  tot <- length(Mn)   # number of points in noise distribution

  gens <- names(Mg)

  Mruido <- abs(round(Mn, prec))
  Druido <- round(Dn, prec)
  Mgen <- abs(round(Mg, prec))
  Dgen <- round(Dg, prec)

  MDgen <- na.omit(cbind(Mgen, Dgen))

  MDunic <- unique(MDgen)  

  Nres <- apply(MDunic, 1, n.menor, S1 = Mruido, S2 = Druido)
  
  lugares <- apply(MDgen, 1, busca, S = MDunic)

  Nconj <- Nres[lugares]

  names(Nconj) <- names(lugares)

  laprob <- Nconj / tot

  laprob <- laprob[gens]
  names(laprob) <- gens

  Nconj <- Nconj[gens]
  names(Nconj) <- gens

  laprob <- list("prob" = laprob, "numDE" = Nconj, "numNOISE" = tot)

  laprob

}

