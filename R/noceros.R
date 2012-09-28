noceros <-
function (x, num = TRUE, k = 0) {
  
  nn <- length(which(x > k))
  
  if (num) {
    nn
    
  } else {
    if(nn > 0) { which(x > k) } else { NULL }
  }
}

