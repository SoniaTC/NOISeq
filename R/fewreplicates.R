####################################################################################################
#########  Algorithm to share information across genes when few replicates are available  ##########
####################################################################################################


## By Sonia Tarazona
## Created: 11-mar-2013




## Function to compute Z for noise when few replicates are available

share.info = function (mydata, n1, n2, r, nclust)   {
  
  # clustering data by k-means algorithm
  # 1.a) Normalized data
  cl = kmeans(mydata, nclust, nstart = 25, iter.max = nclust + 20)
  # 1.b) Log-scale
#   cl = kmeans(log2(mydata+1), nclust, nstart = 25, iter.max = nclust + 20)
  
  print("k-means clustering done")
  print(paste("Size of", nclust, "clusters:"))
  print(cl$size)
  
#   miscolores <- c(colors()[c(554,89,111,512,17,586,132,428,601,568,86,390)], 1:8)
#   pdf("~/Dropbox/noiseq/paperNOISeqBIO/draft/img/clusters.pdf")
#   plot(log2(rowMeans(mydata[,1:n1])+1), log2(rowMeans(mydata[,(n1+1):(n1+n2)])+1), col = miscolores[cl$cluster],
#        pch = cl$cluster)
#   dev.off()
            
  # Creating pseudo-data
  cluster.data = lapply(1:nclust, function (k) { mydata[cl$cluster == k,] })
    
  # Resampling
  npermu = cl$size * r
  
  print("Resampling cluster...")
  
  myres = vector("list", length = nclust)
  
  for (i in 1:nclust) {
    
    print(i)
    
#       if (cl$size[i] > 1) {   # OPTION 2.A
    if (cl$size[i] > 1 && cl$size[i] < 1000) {   # OPTION 2.C: small clusters
      
      myres[[i]] = t(sapply(1:npermu[i], function (j) {
        
        permu = sample(cluster.data[[i]])
        
        nn1 = n1*cl$size[i]
        nn2 = n2*cl$size[i]
        
        mean1 = mean(permu[1:nn1])
        mean2 = mean(permu[(nn1+1):(nn1+nn2)])
        
        sd1 = sd(as.numeric(permu[1:nn1]))
        sd2 = sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
        
        data.frame("M" = log2(mean1/mean2), "D" = mean1-mean2, 
                   "M.sd" = sqrt(sd1^2 / (mean1^2 * log(2)^2 * nn1) + 
                     sd2^2 / (mean2^2 * log(2)^2 * nn2)),
                   "D.sd" = sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))                
        
      }))      
    }
    
    
    if (cl$size[i] >= 1000) {  # OPTION 2.C & 2.D: big clusters
      
#       Option 2.D: clustering big clusters
      cl2 = kmeans(cluster.data[[i]], nclust, nstart = 25, iter.max = nclust + 20)
      
      print(paste("Size of", nclust, "clusters of subcluster:", i))
      print(cl2$size)
      subcluster.data = lapply(1:nclust, function (k) { cluster.data[[i]][cl2$cluster == k,] })
      
      npermu2 = cl2$size * r
      myres2 = vector("list", length = nclust)
      
      for (h in 1:nclust) {
        if (cl2$size[h] > 1) {   
          
          myres2[[h]] = t(sapply(1:npermu2[h], function (j) {
            
            permu = sample(subcluster.data[[h]])
            
            nn1 = n1*cl2$size[h]
            nn2 = n2*cl2$size[h]       
            
            mean1 = mean(permu[1:nn1])
            mean2 = mean(permu[(nn1+1):(nn1+nn2)])
            
            sd1 = sd(as.numeric(permu[1:nn1]))
            sd2 = sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
            
            data.frame("M" = log2(mean1/mean2), "D" = mean1-mean2, 
                       "M.sd" = sqrt(sd1^2 / (mean1^2 * log(2)^2 * nn1) + 
                                       sd2^2 / (mean2^2 * log(2)^2 * nn2)),
                       "D.sd" = sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))              
          }))      
        }
      }
      myres[[i]] = do.call("rbind", myres2)
      
      
      # Option 2.C: dividing big clusters in subclusters of size 500      
#       newsize = 500
#       k1 = floor(cl$size[i]/newsize)
#       selgenes = sample(cl$size[i], newsize*k1)
#       
#       tmp = NULL
#       
#       for (kk in 1:k1) {
#         
#         subsel = selgenes[(1+(kk-1)*newsize):(kk*newsize)]
#         subcluster = cluster.data[[i]][subsel,]
#         
#         subres = t(sapply(1:(newsize*r), function (j) {
#           
#           permu = sample(subcluster)
#           
#           nn1 = n1*newsize
#           nn2 = n2*newsize
#           
#           mean1 = mean(permu[1:nn1])
#           mean2 = mean(permu[(nn1+1):(nn1+nn2)])
#           
#           sd1 = sd(as.numeric(permu[1:nn1]))
#           sd2 = sd(as.numeric(permu[(nn1+1):(nn1+nn2)]))
#           
#           data.frame("M" = log2(mean1/mean2), "D" = mean1-mean2, 
#                      "M.sd" = sqrt(sd1^2 / (mean1^2 * log(2)^2 * nn1) + 
#                        sd2^2 / (mean2^2 * log(2)^2 * nn2)),
#                      "D.sd" = sqrt(sd1^2/sqrt(nn1) + sd2^2/sqrt(nn2)))                           
#         }))
#         
#         tmp = rbind(tmp, subres)                      
#       }
#       myres[[i]] = tmp
      
#     }    
      
  }
    
  }                  
                  
#   
#   print("myres")
#   print(sapply(myres, NROW))
#   
  
  # Computing Zr for noise 
  print("Computing Z for noise...")
  
  # 4.A) a0: Global for all R*G permutations
  myres = do.call("rbind", myres)
#   print(dim(myres))
  
  
  a0.M <- quantile(as.numeric(myres[,"M.sd"]), probs = 0.9, na.rm = TRUE)
  a0.D <- quantile(as.numeric(myres[,"D.sd"]), probs = 0.9, na.rm = TRUE)
  
  M <- as.numeric(myres[,"M"]) / (a0.M + as.numeric(myres[,"M.sd"]))
  D <- as.numeric(myres[,"D"]) / (a0.D + as.numeric(myres[,"D.sd"]))
  (M + D) / 2
  
  
  # 4.B) a0: Per cluster
  
#   myres = sapply(myres, function (x) { 
#     
#               a0.M <- quantile(as.numeric(x[,"M.sd"]), probs = 0.9, na.rm = TRUE)
#               a0.D <- quantile(as.numeric(x[,"D.sd"]), probs = 0.9, na.rm = TRUE)
#               
#               M <- as.numeric(x[,"M"]) / (a0.M + as.numeric(x[,"M.sd"]))
#               D <- as.numeric(x[,"D"]) / (a0.D + as.numeric(x[,"D.sd"]))
#               (M + D) / 2
#           })  
#   unlist(myres)    
  


}



## Toy data

# nr = 2*2
# gk = 1000
# 
# datos = rbind(matrix(rnorm(nr*gk, mean = 0, sd = 1), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 10, sd = 2), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 20, sd = 4), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 50, sd = 5), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 80, sd = 10), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 100, sd = 10), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 150, sd = 12), nrow = gk),
#               matrix(rnorm(nr*gk, mean = 200, sd = 15), nrow = gk))
# 
# datos[which(datos == 0)] = 0.0001
# datos = abs(datos)


## A. fumigatus data

# load("~/Dropbox/transpat_henriette_RNAseq/henriette.RData")
# 
# datos = miscounts[,6:9]
# datos = datos + 1
# 
# mynclust = 12
# 
# miresultado = share.info(mydata = datos, n1 = 2, n2 = 2, r = 5, nclust = mynclust)

