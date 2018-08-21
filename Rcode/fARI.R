# A function to return estiamted labels and ARI(true.label,est.label)
# Y: data to be checked (Y is a mi \times n matrix)
# Y.test.label : true label of Y
# n: number of obs in Y
# K: number of estimated clusters
ARI_f=function(Y,Y.test.label,n,K,Beta.mean,b.mean,sigma.mean,Si,pi.K.mean,mi){
Y.est.label=rep(0,n) 
for (i in 1:n) {
  tau1=tau2=rep(0,K)
  for (k in 1:K){
    #----old mu.star /Sig.star
    #mu.star=Si%*%Beta.mean[[k]]
    #Sig.star=(Si%*%t(Si))*xi.mean[k]+sigma.mean[k]*diag(mi)
    #Sig.star.det.log=as.numeric(determinant(Sig.star,logarithm=TRUE)[[1]])
    #Sig.star.inv=solve(Sig.star)
    #----new mu.star /Sig.star
    mu.star=Si%*%Beta.mean[[k]]+Si%*%b.mean[k,,i]
    Sig.star=sigma.mean[k]*diag(mi)
    Sig.star.det.log=as.numeric(determinant(Sig.star,logarithm=TRUE)[[1]])
    Sig.star.inv=solve(Sig.star)    
    #-----------%%%%%%%%%-------------------------------
    Y.err=Y[,i]-mu.star
    exp.term=t(Y.err)%*%Sig.star.inv%*%Y.err
    tau1[k]=log(pi.K.mean[k])-Sig.star.det.log/2-exp.term/2  #log(pi.k*f(yi,zi=k|...))
  }
  #for (k in 1:K){
  #  tau2[k]=1/sum(exp(tau1-tau1[k]))
  # }
  Y.est.label[i]=which.max(tau1) 
}  
ARI=adjustedRandIndex(Y.test.label,Y.est.label)
return(list(ARI=ARI,est.label=Y.est.label)) 
}
#testARI=ARI_f(Y,Y.test.label,n,K,Beta.mean,b.mean,sigma.mean,Si,pi.K.mean,mi)
ARI_f2=function(Y,Y.test.label,n,K,Beta.mean,xi.mean,sigma.mean,Si,pi.K.mean,mi){
  Y.est.label=rep(0,n) 
  for (i in 1:n) {
    tau1=tau2=rep(0,K)
    for (k in 1:K){
      #----old mu.star /Sig.star
      mu.star=Si%*%Beta.mean[[k]]
      Sig.star=(Si%*%t(Si))*xi.mean[k]+sigma.mean[k]*diag(mi)
      Sig.star.det.log=as.numeric(determinant(Sig.star,logarithm=TRUE)[[1]])
      Sig.star.inv=solve(Sig.star)
      #----new mu.star /Sig.star
      #mu.star=Si%*%Beta.mean[[k]]+Si%*%b.mean[k,,i]
      #Sig.star=sigma.mean[k]*diag(mi)
      #Sig.star.det.log=as.numeric(determinant(Sig.star,logarithm=TRUE)[[1]])
      #Sig.star.inv=solve(Sig.star)    
      #-----------%%%%%%%%%-------------------------------
      Y.err=Y[,i]-mu.star
      exp.term=t(Y.err)%*%Sig.star.inv%*%Y.err
      tau1[k]=log(pi.K.mean[k])-Sig.star.det.log/2-exp.term/2  #log(pi.k*f(yi,zi=k|...))
    }
    #for (k in 1:K){
    #  tau2[k]=1/sum(exp(tau1-tau1[k]))
    # }
    Y.est.label[i]=which.max(tau1) 
  }  
  ARI=adjustedRandIndex(Y.test.label,Y.est.label)
  return(list(ARI=ARI,est.label=Y.est.label)) 
}
#ARI_f2(Y,Y.test.label,n,K,Beta.mean,xi.mean,sigma.mean,Si,pi.K.mean,mi)
  