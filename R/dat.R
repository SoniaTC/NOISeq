##### Function to generate data for exploratory plots #####

# By Sonia & Pedro
# Modified: 10-sep-13

dat = function (input, type = c("biodetection","cd","countsbio","GCbias","lengthbias","saturation"), 
                k = 0, ndepth = 6, factor = NULL, norm = FALSE, refColumn = 1) {
  
  type <- match.arg(type)
  
  if (type == "biodetection") {
    
    output = new("Biodetection", dat = biodetection.dat(input, factor = factor, k = k))    
  }
  
  if (type == "cd") {
    
    output = new("CD", dat = cd.dat(input, norm = norm, refColumn = refColumn))
  }
  
  if (type == "countsbio") {
    
    output =  new("CountsBio", dat = countsbio.dat(input, factor = factor, norm = norm))
  }


  if (type == "GCbias") {
    
    output =  new("GCbias", dat = GC.dat(input, factor = factor, norm = norm))
  }

  
  if (type == "lengthbias") {
    
    output =  new("lengthbias", dat = length.dat(input, factor = factor, norm = norm))
  }
  
  if (type == "saturation") {
    
    output = new("Saturation", dat = saturation.dat(input, k = k, ndepth = ndepth))
  }
  
  
  output  
  
}
