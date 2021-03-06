ASCAfun2<-function (X,Desa,Desb,Fac) {

n <- ncol(X) # number of genes
I <- ncol(Desa) # number of levels in the factor TIME
J <- ncol(Desb) # number of levels in the other factor

XK1<-matrix(NA,nrow=I,ncol=n)

for (i in 1:I) {
     sub<-X[Desa[,i]==1,]
 if(is.null(nrow(sub))){ #when there isn't replicates
    XK1[i,] <- sub 
}else{
     XK1[i,]<-apply(sub,2,mean) }
}

XK2<-matrix(NA,nrow=J,ncol=n)

for (j in 1:J) {
     sub<-X[Desb[,j]==1,]
 if(is.null(nrow(sub))){ #when there isn't replicates
    XK2[j,] <- sub 
}else{
     XK2[j,]<-apply(sub,2,mean) }
}

NK<-matrix(NA,nrow=I,ncol=J)
XK<-matrix(NA,nrow=I*J,ncol=n)

k=1
for (j in 1:J){
  for (i in 1:I){
    sub<-X[(Desa[,i]+Desb[,j])==2,]
if(is.null(nrow(sub))){  #when there isn't replicates
    NK[i,j] <- 1
    XK[k,] <- sub-XK1[i,]-XK2[j,]
}else{
    NK[i,j]<-sqrt(nrow(sub))
    XK[k,]<-apply(sub,2,mean)-XK1[i,]-XK2[j,] }
    k=k+1
  }
}

XKw<-XK*(as.numeric(NK))

PCA<-PCA.GENES(XKw)
      scw<-PCA$scores[,1:Fac]
      ld<-PCA$loadings[,1:Fac]
      ssq<-PCA$var.exp
      if(Fac==1) {
      scw<-as.matrix(scw)
      ld<-as.matrix(ld) 
      }      
 if(Fac==0) {
      scw<-as.matrix(rep(0,I*J))
      ld<-as.matrix(rep(0,n))
	}           
# Re-weigth the scores
sc<-scw/(as.numeric(NK))
    
XKrec<-sc%*%t(ld)

Xab<-NULL
TPab<-NULL
for (i in 1:nrow(X)){
     position1<-which(Desa[i,]==1)
     position2<-which(Desb[i,]==1)
     Xab<-rbind(Xab,XK[I*(position2-1)+position1,])
     TPab<-rbind(TPab,XKrec[I*(position2-1)+position1,])
}
Eab<-Xab-TPab

 leverage<-apply(ld^2,1,sum)
 SPE<-apply(Eab^2,2,sum)

output<-list(XK,sc,ld,ssq,Xab,TPab,Eab,leverage,SPE)
names(output)<-c("data","scores","loadings","var.exp","X","TP","E","leverage","SPE")
output


}
