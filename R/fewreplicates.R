####################################################################################################
#########  Algorithm to share information across genes when few replicates are available  ##########
####################################################################################################


## By Sonia Tarazona
## Created: 11-mar-2013




## Function to compute Z for noise when few replicates are available

share.info = function (mydata, n1, n2, r, nclust)   {
  
  # clustering data by k-means algorithm
  # 1.a) Normalized data
  gc()
  cl = suppressWarnings(kmeans(mydata, nclust, nstart = 25, iter.max = nclust + 30))

  cat("...k-means clustering done\n")
  cat(paste("Size of", nclust, "clusters:\n"))
  print(cl$size)
            
  # Creating pseudo-data
  cluster.data = lapply(1:nclust, function (k) { mydata[cl$cluster == k,] })
    
  # Resampling
  npermu = cl$size * r
  npermu = sapply(npermu, function (x) min(x, 1000))  ## modified to reduce the number of permutations to be done
  
  cat("Resampling cluster...")
  
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
      
      cat(paste("Size of", nclust, "subclusters of cluster:", i)); cat("\n")
      print(cl2$size)
      subcluster.data = lapply(1:nclust, function (k) { cluster.data[[i]][cl2$cluster == k,] })
      
      npermu2 = cl2$size * r
      npermu2 = sapply(npermu2, function (x) min(x, 1000))  ## modified to reduce the number of permutations to be done
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
    }
  }                    

  
  # Computing Zr for noise 
  cat("Computing Z for noise...\n")
  
  # 4.A) a0: Global for all R*G permutations
  myres = do.call("rbind", myres)
  
  
  a0.M <- quantile(as.numeric(myres[,"M.sd"]), probs = 0.9, na.rm = TRUE)
  a0.D <- quantile(as.numeric(myres[,"D.sd"]), probs = 0.9, na.rm = TRUE)
  
  M <- as.numeric(myres[,"M"]) / (a0.M + as.numeric(myres[,"M.sd"]))
  D <- as.numeric(myres[,"D"]) / (a0.D + as.numeric(myres[,"D.sd"]))
  (M + D) / 2


}



