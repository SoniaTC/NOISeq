sinceros <-
function (datos, k) {
  datos0 <- datos

  if (is.null(k)) {

    mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])

    kc <- mini0/2

    datos0[datos0 == 0] <- kc

  } else {

    datos0[datos0 == 0] <- k

  }
  
  datos0
}

