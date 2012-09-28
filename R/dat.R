##### Function to generate data for exploratory plots #####

# By Sonia & Pedro
# Created: 20-sep-12

dat = function (input, type = c("biodetection","cd","countsbio","DLbio","saturation"), 
                selection=c(1,2), k = 0, ndepth = 5, newdetections = TRUE) {
  
  type <- match.arg(type)
  
  if (type == "biodetection") {
    
    output = new("Biodetection", dat = biodetection.dat(input, data_selection=selection, k = k))    
  }
  
  if (type == "cd") {
    
    output = new("CD", dat = cd.dat(input, columns = selection))
  }
  
  if (type == "countsbio") {
    
    output =  new("CountsBio", dat = countsbio.dat(input, cols = selection, k = k,
                                 ndepth = ndepth, quartiles = TRUE))
  }
  
  if (type == "DLbio") {
    
    output =  new("DLBio", dat = DLbio.dat(input, k = k, ndepth = ndepth))
  }
  
  if (type == "saturation") {
    
    output = new("Saturation", dat = saturation.dat(input, k = k, ndepth = ndepth, newdetections = newdetections))
  }
  
  
  output  
  
}
