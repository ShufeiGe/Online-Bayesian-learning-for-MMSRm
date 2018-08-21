#step1. Calculate the dissimilarity matrix
dissm=daisy(t(Y[,1:n.min]), metric ="euclidean",
              #c("euclidean", "manhattan", "gower"),
      stand = FALSE, type = list())
#step2. Get the tree based on the hierichical clustering
temp.cluster=hclust(dissm, method = "complete", members = NULL)
#step3. Get the cluster labels
label.temp=cutree(temp.cluster,k=K)

#step4. Get the initial values.
Beta=array(1,c(K,d,1))
pi.K=rep(NA,K)
for( k in 1:K){
  Beta[k,,]=solve(t(Si)%*%Si)%*%rowMeans(t(Si)%*%Y[,which(label.temp==k)])
  pi.K[k]=length(which(label.temp==k))/dim(Y[,1:n.min])[2]
}
