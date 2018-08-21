# simulation study example :
#          ----surface clustering via smc(Gibbs sampleing) algorithm----
#          ----see paper section '5.2 Comparison of online SMC with MCMC'
# setp 1:  simulate surfaces from function fk()
rm(list=ls())
path="/Users/gesophie/Desktop/Github/R_package/"
setwd(path)
#please note that the invgamma() here is the one in package "LaplacesDemon" ;
#so don't call library "invgamma",  otherwise it will lead mistakes
library(MASS)
library(LaplacesDemon)
library(mclust)
library(matrixcalc)
library(mvtnorm)
library(cluster)
library(plotly)
seed=859859
set.seed(seed)
n = 10000  # total no. of observations
n.max=10000
n.min=20   # no. of obs gonna be used in this simulation.
#n.max=10000
K =3  #clustering exampler when K=3 
#f_k=1
pi.true=rep((1/K),K)
Z <- sample(1:K, size = n, replace = TRUE, prob = pi.true)
mi=12*12
source("NBFs.R")
source("simulator.R")

# setp 2:  surface clustering via smc algorithm
# parameters initialization
out_mcmc=paste("int_mcmc_clustering_K_",K,".RData",sep="") # output .R file name
source("fARI.R")
source("Algorithm1_k.R")
Beta.mean=list()
b.mean=result_b
for (k in 1:K){
  Beta.mean[[k]]=apply(result_Beta[k,(burn_in+1):it_max,,],2,mean)
}


pi.K.mean=apply(as.matrix(result_piK[(burn_in+1):it_max,]),2,mean)
sigma.mean=apply(as.matrix(result_sigma[(burn_in+1):it_max,]),2,mean)
xi.mean=apply(as.matrix(result_xi[(burn_in+1):it_max,]),2,mean)

## Assign labels and calculated ARI
Y.est.label=rep(0,n)
Y.test=Y
Y.test.label=Y.label[1:n.min]
ARI.all=ARI_f2(Y,Y.test.label,n,K,Beta.mean,xi.mean,sigma.mean,Si,pi.K.mean,mi)
#ARI.all=ARI_f(Y,Y.test.label,n,K,Beta.mean,b.mean,sigma.mean,Si,pi.K.mean,mi)
Y.est.label=ARI.all$est.label
ARI.init=ARI.all$ARI

Beta_init=matrix(1,nr = d, nc  = K)
for (k in 1:K){
  Beta_init[,k]=Beta.mean[[k]]  
}

b_init=b.mean
xi_init=xi.mean
sigma_init=sigma.mean
pi_init=pi.K.mean

#save(Beta_init, b_init, xi_init,sigma_init,pi_init,file="par_init.RData")
#load("par_init.RData")
Si2=t(Si)%*%Si                              #Si_Si=  t(Si)%*%Si
Sit2=Si%*%t(Si)                             #Sit2=  Si%*%t(Si)

NP = 1000
N  = dim(Y.all)[2]
source("PF.R")
y=Y.all
y.true.label=Y.label
pi=pi_init
sigma=sigma_init
NP 
K
beta=Beta_init
xi=xi_init
b=b_init
###-----  update the first n.min -----####
out <- pf1(y, pi, beta, sigma,xi, n.min, NP, K)
if(K==1){
  pi=as.matrix(t(out$pi)) 
}else{
  pi=out$pi 
}

beta=out$beta
sigma=out$sigma
xi=out$xi
b2=out$b2
nk=out$nk
Syb=out$Syb
Syb2=out$Syb2
t1=n.min+1
t2=t1+500  # for test ; 
#t2=n.max
## in order to simulate no. of observations from n.min to n.max,
## we need to  set t2=n.max


pb <- txtProgressBar(min = 0, max = (t2-t1+1), style = 3)
time_cost=rep(0,(t2-t1+1))

# since computing cost of ARI function is high, we only compute it 
# every 20 steps
ARI=rep(NA,(t2-t1+1))
xi.smc =matrix(NA,nr=(t2-t1+1),nc=K)
sig.smc=matrix(NA,nr=(t2-t1+1),nc=K)
pik.smc=matrix(NA,nr=(t2-t1+1),nc=K)

xi.lci.smc =matrix(NA,nr=(t2-t1+1),nc=K)
sig.lci.smc =matrix(NA,nr=(t2-t1+1),nc=K)
pik.lci.smc =matrix(NA,nr=(t2-t1+1),nc=K)

xi.uci.smc =matrix(NA,nr=(t2-t1+1),nc=K)
sig.uci.smc =matrix(NA,nr=(t2-t1+1),nc=K)
pik.uci.smc =matrix(NA,nr=(t2-t1+1),nc=K)

beta.rd.smc.uci=matrix(NA,nr=(t2-t1+1),nc=K)
beta.rd.smc.lci=matrix(NA,nr=(t2-t1+1),nc=K) 
beta.rd.smc=matrix(NA,nr=(t2-t1+1),nc=K) 
beta.smc.ind=sample(1:d,K)

b2.smc.m =matrix(NA,nr=(t2-t1+1),nc=K) 
#here we only save beta1[ind1],beta2[ind2],...,betak[indk]
for ( t in t1 : t2){
  rl.t=(t-t1+1)  
  setTxtProgressBar(pb, rl.t)
  t0=Sys.time()
  yt=y[,t]
  out=pf2(yt,t,NP,K,pi,beta,sigma,xi,b2,nk,Syb,Syb2,Si2,Sit2)
  time_cost[rl.t]=Sys.time()-t0
  #pi=out$pi
  #pi=t(as.matrix(out$pi))
  if(K==1){
    pi=as.matrix(t(out$pi)) 
  }else{
    pi=out$pi 
  }
  beta=out$beta
  sigma=out$sigma
  xi=out$xi
  b2=out$b2
  nk=out$nk
  Syb=out$Syb
  Syb2=out$Syb2
  #-----------save some important pars and calculate ARI for each step-------
  beta.m=list()
  for (k in 1:K){
    beta.m[[k]]=apply(beta[[k]],2,mean) 
    beta.rd.smc[rl.t,k]=beta.m[[k]][beta.smc.ind[k]] #save certain elements in beta
    beta.rd.smc.lci[rl.t,k]=quantile(beta[[k]][,beta.smc.ind[k]] ,0.025)
    beta.rd.smc.uci[rl.t,k]=quantile(beta[[k]][,beta.smc.ind[k]] ,0.975)
    xi.lci.smc [rl.t,k]=quantile(xi[,k]   ,0.025)
    sig.lci.smc[rl.t,k]=quantile(sigma[,k],0.025)
    pik.lci.smc[rl.t,k]=quantile(pi[,k] ,0.025)
    xi.uci.smc [rl.t,k]=quantile(xi[,k] ,0.975)
    sig.uci.smc[rl.t,k]=quantile(sigma[,k],0.975)
    pik.uci.smc[rl.t,k]=quantile(pi[,k] ,0.975)
  }
  pi.K.m=apply(pi,2,mean)
  sig.m=apply(sigma,2,mean)
  xi.m=apply(xi,2,mean)
  # ---save important pars---
  xi.smc[rl.t,]=xi.m
  sig.smc[rl.t,]=sig.m
  pik.smc[rl.t,]=pi.K.m
  b2.smc.m[rl.t,]=apply(b2,2,mean)
  # ---get ARI-----
  get.ARI=(t%%100==0)+(t==t1)+(t==t2)
  if(get.ARI){
    y.tmp=y[,1:t]
    y.tmp.l=y.true.label[1:t]
    ARI.all=ARI_f2(y.tmp,y.tmp.l,t,K,beta.m,xi.m,sig.m,Si,pi.K.m,mi)
    ARI[rl.t]=ARI.all$ARI  
  }
  if(t%%200==0){
    save.image(file="SMC_Sim_K3.RData") 
  }
}

save.image(file="SMC_Sim_K3.RData") 



par(mfrow=c(K,2))
for(k in 1:K){
  plot(x=t1:t2,sig.smc[,k],type="l",xlab="t",ylab=paste("Avg.sig^2.k: k=",k,sep="")) 
  hist(sigma[,k],freq=FALSE,main=paste("sig.k: k=",k,"at time =",t2,sep="")) 
  abline(v=0.1,lty=1,col="red")
  abline(v=mean(c(sigma[,k])),lty=1,col="blue")
  abline(v=quantile(c(sigma[,k]),probs=c(0.025,0.975)),lty=2,col="blue")
}

for(k in 1:K){
  plot(x=t1:t2,pik.smc[,k],type="l",xlab="t",ylab=paste("Avg.pi.k: k=",k,sep="")) 
  hist(pi[,k],freq=FALSE,main=paste("pi.k: k=",k,"at time =",t2,sep=""))
  abline(v=1/K,lty=1,col="red")
  abline(v=mean(c(pi[,k])),lty=1,col="blue")
  abline(v=quantile(c(pi[,k]),probs=c(0.025,0.975)),lty=2,col="blue")
}

for(k in 1:K){
  plot(x=t1:t2,xi.smc[,k],type="l",xlab="t",ylab=paste("Avg.xi^2.k: k=",k,sep="")) 
  hist(xi[,k],freq=FALSE,main=paste("pi.k: k=",k,"at time =",t2,sep=""))
  abline(v=.3,lty=1,col="red")
  abline(v=mean(c(xi[,k])),lty=1,col="blue")
  abline(v=quantile(c(xi[,k]),probs=c(0.025,0.975)),lty=2,col="blue")
}

