MD.plot <-
function (Ms = Ms, Ds = Ds, Mn = Mn, Dn = Dn, 
          xlim = range(Ms,na.rm=TRUE), ylim = c(0,quantile(Ds,0.8,na.rm=TRUE)),
          tit = "") {
  plot(Mn, Dn, pch = ".", main = tit, xlab = "M", ylab = "D",
       xlim = xlim, ylim = ylim)
  points(Ms, Ds, col = 2, pch = ".")

  legend("topright", c("noise", "signal"), col = 1:2, pch = 15,
         bg = "lightgrey")

}

