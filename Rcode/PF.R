#library("matrixcalc")
###Gibbs: pi ###
piGibbs <- function(nn, alpha){
  pi.out=as.vector(rdirichlet(1, nn + alpha))
  return(pi.out)
}

###Gibbs: beta ###
betaGibbs <- function(NP,K,nk,Syb,Sigma0,u0,sigma){
beta_out=list()
for (k in 1:K){
  beta_out[[k]]=matrix(0,nr=NP,ncol=d)
}
for (k in 1:K){
 for ( i in 1:NP){
  sig.betak=solve((solve(Sigma0)+nk[i,k]*Si2/sigma[k])) 
  mu.betak=sig.betak%*%(solve(Sigma0)%*%u0+c(Syb[[k]][i,])/sigma[k])
  beta_out[[k]][i,]=rmvnorm(1,mean=mu.betak,sigma=sig.betak)
 }
}  
 return(beta_out) 
}
betaGibbs2 <- function(NP,K,nk,Syb,Sigma0,u0,sigma){
  beta_out=list()
  for (k in 1:K){
    beta_out[[k]]=matrix(0,nr=NP,ncol=d)
  }
  for (k in 1:K){
    for ( i in 1:NP){
      sig.betak=solve((solve(Sigma0)+nk[i,k]*Si2/sigma[i,k])) 
      mu.betak=sig.betak%*%(solve(Sigma0)%*%u0+c(Syb[[k]][i,])/sigma[i,k])
      beta_out[[k]][i,]=rmvnorm(1,mean=mu.betak,sigma=sig.betak)
    }
  }  
  return(beta_out) 
}
###Gibbs: sigma ###
sigmaGibbs <- function(NP,K,nk,Syb,Syb2,beta,g0,h0){
sigma_out=matrix(0,nr=NP,nc=K)
g=g0+nk*mi/2
for (k in 1:K){
  for (i in 1:NP){
   gi=g[i,k]
   hi=h0+(Syb2[i,k]-2*beta[[k]][i,]%*%Syb[[k]][i,]+nk[i,k]*t(beta[[k]][i,])%*%Si2%*%beta[[k]][i,])/2
   sigma_out[i,k]= rinvgamma(1,shape=gi,scale=hi)
  }
}
  return(sigma_out)
}

###Gibbs: xi ###
xiGibbs <- function(NP,K,t,b2,a0,b0){
  xi_out=matrix(0,nr=NP,nc=K)
  a=a0+t*d/2
  b=b0+b2/2
  for (k in 1:K){
    for (i in 1:NP){
    xi_out[i,k]=rinvgamma(1,shape=a,scale=b[i,k])
    }
  }
  return(xi_out)
}

# W : weight of All particles at time {t-1}
# compute weights of all particles at time {t-1} conditional on theta_{t-1}
# i.e P(y_{t}|theta_{t-1})
# computWeight(NP,K,t,Si,Sit2,pi,beta,xi,sigma)
computWeight=function(NP,K,t,Si,Sit2,pi,beta,xi,sigma){
  w=matrix(0,nr=NP,nc=K)
  W=rep(0,NP)
  for (i in 1:NP){
    for (k in 1:K){
      muk.yi=Si%*%beta[[k]][i,]
      sigk.yi=Sit2*xi[i,k]+sigma[i,k]*diag(mi)  
      w[i,k]=dmvnorm(y[,t],mean=muk.yi,sigma=sigk.yi,log=TRUE)  
    }
  }
  W=rowSums(exp(w-max(w))*pi)/sum(exp(w-max(w))*pi) 
  return(W) 
}
##--------------------------------------------------------------------><<<<<<<<<<<<#
computWeight_v2=function(NP,K,t,Si,Sit2,pi,beta,xi,sigma){
    w=matrix(0,nr=NP,nc=K)
    W.normlized=rep(0,NP)
    for (i in 1:NP){
        for (k in 1:K){
            muk.yi=Si%*%beta[[k]][i,]
            sigk.yi=Sit2*xi[i,k]+sigma[i,k]*diag(mi)
            w[i,k]=dmvnorm(y[,t],mean=muk.yi,sigma=sigk.yi,log=TRUE)
        }
    }
    W.normlized=rowSums(exp(w-max(w))*pi)/sum(exp(w-max(w))*pi)
    #output list
    py=rowSums(exp(w)*pi)
    out_put=list(W=W.normlized,py=py)
    return(out_put)
}

##--------------------------------------------------------------------><<<<<<<<<<<<#

##--------------------------------------------------------------------><<<<<<<<<<<<#
##(1) input (y, pi, beta, sigma, xi, n.min, NP,K)
##    NP : # of particles
##     K : # of clusters
## n.min : 
##    y  : nrow(y)=mi, ncol(y)=n
##  pi   : (pi1,...,piK)
##  xi   : (xi1,...,xiK)
##  beta : nrow(beta)=d, nrow(beta)=K
##(2) output : list(pi=pivec,beta=beta_list,sigma=sigma_out,xi=xi_out,
##           b2=b2,nk=nk,Syb=Syb,Syb2=Syb2)
##      pi: nr=NP, nc=K   ;   beta : K  X [NP X d]
##   sigma: NP X K        ;   xi   : NP X K
##      b2: NP X K        ;   nk   : NP X K
##     Syb: k  X [NP X d] ;   Syb2 : NP X K
pf1 <- function(y, pi, beta, sigma,xi, n.min, NP, K){
  
	nk <- matrix(0, nr = NP, nc = K)
	#yk2 <- matrix(0, nr = NP, nc = K)
	#yk <- matrix(0, nr = NP, nc = K)
	b<- list()
	b2<-matrix(0,nr=NP,nc=K)
	Syb<- list()
	beta_out<-list()
	Syb2<- matrix(0,nr=NP,nc=K)
	for (k in 1:K){
	Syb[[k]]<-matrix(0,nr=NP,nc=d)  
	beta_out[[k]]<-matrix(0,nr=NP,nc=d)  
	}
#	LhMatrix <- matrix(NA, nr = K, nc = n.min)
#	PosteriorMatrix <- matrix(NA, nr = K, nc = n.min)
	logLhMatrix_K <- rep(NA, K)
	logZ.py <- rep(NA, n.min)	
	LhMatrix <- matrix(NA, nr = n.min, nc = K)
	PosteriorMatrix <-matrix(NA, nr = n.min, nc = K)
	pb <- txtProgressBar(min = 0, max =n.min, style = 3)
	for(j in 1:n.min){
	  setTxtProgressBar(pb, j)
	  #logLhMatrix_K = -log(sigma)*mi/2-sum((y[,j]-mu)^2)/sigma 
	  for (k in 1:K){
	    muk.y=Si%*%beta[,k]
	    sigk.y=Sit2*xi[k]+sigma[k]*diag(mi)
	    logLhMatrix_K[k]=dmvnorm(y[,j],mean=muk.y,sigma=sigk.y,log=TRUE)  
	  }
	  
	  if(j==1){
	    logZ.py[j]=log(sum(pi*exp(logLhMatrix_K)))
	  }else{
	    logZ.py[j]=log(sum(pi*exp(logLhMatrix_K)))+ logZ.py[(j-1)]
	  }

	  LhMatrix[j,]=exp(logLhMatrix_K-max(logLhMatrix_K))
	  PosteriorMatrix[j,] <- pi*LhMatrix[j,]/sum(pi*LhMatrix[j,])
	  PosteriorMatrix[j,]
	  Zvec <- sample.int(K, NP, replace = TRUE, prob = PosteriorMatrix[j,])

	  for(i in 1:NP){
	    Zt=Zvec[i]
	    for (k in 1:K){
	      sig_bt=solve(1/xi[k]*diag(d)+(Zt==k)*Si2/sigma[k])
	      mu_bt=sig_bt%*%(t(Si)%*%(y[,j]-Si%*%beta[,k]))*(Zt==k)/sigma[k]
	      bt=rmvnorm(1,mean=mu_bt,sigma=sig_bt)
	      #sig_bk.not=1/xi[k]*diag(d)
	      #mu_bk.not=0	  
	      b2[i,k]=sum(bt^2)+b2[i,k]
	      nk[i,k]=nk[i,k]+1*(Zt==k)
	      
	      #   Syb[[k]][i,]=Syb[[k]][i,]+t(Si)%*%(yt-Si%*%t(bt))   
	      Syb[[k]][i,]=Syb[[k]][i,]+t(Si)%*%(y[,j]-Si%*%t(bt))*(Zt==k)
	      #   Syb2[i,k]=Syb2[i,k]+t(yt-Si%*%t(bt))%*%(yt-Si%*%t(bt))      
	      Syb2[i,k]=Syb2[i,k]+t(y[,j]-Si%*%t(bt))%*%(y[,j]-Si%*%t(bt))*(Zt==k) 
	    }
	  }
	  ## Z_unique=unique(Zvec)
	  ## k_id=list()	  
	  ##   for (k in Z_unique){
	  ##   k_id[[k]]=which(Zvec==k) 
	  ##   sig_bk=solve(1/xi[k]*diag(d)+Si2/sigma[k])
	  ##   mu_bk=sig_bk%*%(t(Si)%*%(y[,j]-Si%*%beta[,k]))/sigma[k]
	  ##   #sig_bk.not=1/xi[k]*diag(d)
	  ##   #mu_bk.not=0	 
	  ##   
	  ##   b[[k]]=rmvnorm(length(k_id[[k]]),mu_bk,sig_bk)
	  ##   b2[k_id[[k]],k]=b2[k_id[[k]],k]+rowSums(b[[k]]^2)
	  ##   #b2[k_id[[k]],k]=rowSums(b[[k]]^2)
	  ##   
	  ##   ## each row of Syb[[k]] : delta(zi=k)%*%t(Si)(yi-Sibik)
	  ##   ## Syb[[k]] : N row * d col, each row corresponds one particle.
	  ##   Syb[[k]][k_id[[k]],]=Syb[[k]][k_id[[k]],]+t(t(Si)%*%(y[,j]-Si%*%t(b[[k]])))
	  ##   ## each element of Syb2 : delta(zi=k)%*%t((yi-Sibik))%*%(yi-Sibik)
	  ##   ## Syb[[k]] : N row * d col, each row corresponds one particle.
	  ##   Syb2[k_id[[k]],k]=Syb2[k_id[[k]],k]+colSums((y[,j]-Si%*%t(b[[k]]))^2) 
	  ##   }
	  ##   nk[cbind(1:NP,Zvec)] <- 1 + nk[cbind(1:NP,Zvec)]
	}
	#mt <- cbind(b2,nk,Syb2); Syb
  pivec <- t(apply(nk, 1, function(x)(piGibbs(x, alpha0))))
  beta_list <-betaGibbs(NP,K,nk,Syb,Sigma0,u0,sigma)
  sigma_out <- sigmaGibbs(NP,K,nk,Syb,Syb2,beta_list,g0,h0)
  xi_out<-xiGibbs(NP,K,n.min,b2,a0,b0)
  
  #output list
  out_put=list(pi=pivec,beta=beta_list,sigma=sigma_out,xi=xi_out,
               b2=b2,nk=nk,Syb=Syb,Syb2=Syb2,logZ.py=logZ.py)
	return(out_put)
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##(1) input : list(yt,t,NP, K,pi,beta,sigma,xi,b2,nk,Syb,Syb2,Si2,Sit2)
##      yt: y[,t]         ;;;     t  : t
##      NP:               ;;;     K  :
##      pi: nr=NP, nc=K   ;;;   beta : K  X [NP X d]
##   sigma: NP X K        ;;;   xi   : NP X K
##      b2: NP X K        ;;;   nk   : NP X K
##     Syb: k  X [NP X d] ;;;   Syb2 : NP X K
##(2) output : list(yt,t,NP, K,pi,beta,sigma,xi,b2,nk,Syb,Syb2,Si2,Sit2)
##      yt: y[,t]         ;;;     t  : t
##      NP:               ;;;     K  :
##      pi: nr=NP, nc=K   ;;;   beta : K  X [NP X d]
##   sigma: NP X K        ;;;   xi   : NP X K
##      b2: NP X K        ;;;   nk   : NP X K
##     Syb: k  X [NP X d] ;;;   Syb2 : NP X K
pf2 <- function(yt,t,NP,K,pi,beta,sigma,xi,b2,nk,Syb,Syb2,Si2,Sit2){
  #compute weights of all particles 1 to NP, at time t-1
  #i.e. W_{t-1} : P(y_{t}|theta_{t-1})
  WResult=computWeight_v2(NP,K,t,Si,Sit2,pi,beta,xi,sigma)
  W=WResult$W
  py=WResult$py
  #sample index
  samp.ind=sample(x=c(1:NP),size=NP,prob=W,replace=TRUE)
# ----------update corresponding sufficient statistics--------
  pi=as.matrix(pi[samp.ind,])
  sigma=as.matrix(sigma[samp.ind,])
  xi=as.matrix(xi[samp.ind,])
  b2=as.matrix(b2[samp.ind,])
  nk=as.matrix(nk[samp.ind,])
  Syb2=as.matrix(Syb2[samp.ind,])
  for (k in 1:K){
  beta[[k]]=beta[[k]][samp.ind,]
  Syb[[k]]=Syb[[k]][samp.ind,]  
  }
# -----------------------------------------------------------
 for(i in 1:NP){
    logLhMatrix_K     <- rep(NA, K)
    LhMatrix        <- rep(NA, K)
    PosteriorMatrix <- rep(NA, K)
   # id=samp.ind[i]
    for (k in 1:K){
      muk.y=Si%*%beta[[k]][i,]
      sigk.y=Sit2*xi[i,k]+sigma[i,k]*diag(mi)
      logLhMatrix_K[k]=dmvnorm(yt,mean=muk.y,sigma=sigk.y,log=TRUE)  
    }
    LhMatrix=exp(logLhMatrix_K-max(logLhMatrix_K))
    PosteriorMatrix<- pi[i,]*LhMatrix/sum(pi[i,]*LhMatrix)
    Zt <- sample.int(K, 1, replace = TRUE, prob = PosteriorMatrix)
    for (k in 1:K){
      sig_bt=solve(1/xi[i,k]*diag(d)+(Zt==k)*Si2/sigma[i,k])
      mu_bt=sig_bt%*%(t(Si)%*%(yt-Si%*%beta[[k]][i,]))*(Zt==k)/sigma[i,k]
      bt=rmvnorm(1,mean=mu_bt,sigma=sig_bt)
      #sig_bk.not=1/xi[k]*diag(d)
      #mu_bk.not=0	  
      b2[i,k]=sum(bt^2)+b2[i,k]
      nk[i,k]=nk[i,k]+1*(Zt==k)      
#   Syb[[k]][i,]=Syb[[k]][i,]+t(Si)%*%(yt-Si%*%t(bt))   
      Syb[[k]][i,]=Syb[[k]][i,]+t(Si)%*%(yt-Si%*%t(bt))*(Zt==k)
#   Syb2[i,k]=Syb2[i,k]+t(yt-Si%*%t(bt))%*%(yt-Si%*%t(bt))      
      Syb2[i,k]=Syb2[i,k]+t(yt-Si%*%t(bt))%*%(yt-Si%*%t(bt))*(Zt==k) 
    }
 }
  
  pi    <- t(apply(nk, 1, function(x)(piGibbs(x, alpha0))))
  beta  <- betaGibbs2(NP,K,nk,Syb,Sigma0,u0,sigma)
  sigma <- sigmaGibbs(NP,K,nk,Syb,Syb2,beta,g0,h0)
  xi    <- xiGibbs(NP,K,t,b2,a0,b0)
  
  out_put=list(pi=pi,beta=beta,sigma=sigma,xi=xi,
               b2=b2,nk=nk,Syb=Syb,Syb2=Syb2,py=py)
  return(out_put)  # added on 4th Apr, 2018
  }  
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>











