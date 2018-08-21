# simulation study example :
#          ----surface clustering via smc(Gibbs sampleing) algorithm----
#          ----see paper section '5.2 Comparison of online SMC with MCMC'
# setp 1:  simulate surfaces from function fk()
rm(list=ls())
path="/Users/gesophie/Desktop/Github/R_package/"
setwd(path)
#please note that the invgamma() here is the one in package "LaplacesDemon" ;
#so don't call library "invgamma",  otherwise it will lead mistakes
rm(list=ls())         
library(MASS)         
library(LaplacesDemon)
library(mclust)       
library(matrixcalc)   
library(mvtnorm)      
library(cluster)      
n = 5000              
n.max=n               
n.min=100             

seed=58068
K=8
d1=8
d2=8
d=d1*d2

out_mcmc="Initial_mcmc_out_example.RData"
out_smc="smc_out_example.RData"
source("read_in_real_app1.R")
source("fARI.R")
# initialize parmaters
source("Algorithm1_k.R")
xi=result_xi
pi=result_piK
Beta=result_Beta
pi.K=result_piK
sigma=result_sigma
Beta.mean=list()
b.mean=b
for (k in 1:K){
  Beta.mean[[k]]=apply(Beta[k,(burn_in+1):it_max,,],2,mean)  
}

pi.K.mean=apply(as.matrix(pi.K[(burn_in+1):it_max,]),2,mean)
sigma.mean=apply(as.matrix(sigma[(burn_in+1):it_max,]),2,mean)
xi.mean=apply(as.matrix(result_xi[(burn_in+1):it_max,]),2,mean)
Beta_init=matrix(1,nr = d, nc  = K)
for (k in 1:K){
  Beta_init[,k]=Beta.mean[[k]]  
}
b_init=b.mean
xi_init=xi.mean
sigma_init=sigma.mean
pi_init=pi.K.mean

Si2=t(Si)%*%Si                                 #Si_Si=  t(Si)%*%Si
Sit2=Si%*%t(Si)                             #Sit2=  Si%*%t(Si)

NP = 1000   # in real data application, we set NP=1000
N=dim(Y.all)[2]

source("PF.R")
y=Y.all
y.true.label=Y.all.label
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
logZ.py=out$logZ.py

beta=out$beta
sigma=out$sigma
xi=out$xi
b2=out$b2
nk=out$nk
Syb=out$Syb
Syb2=out$Syb2

t1=n.min+1
t2=n.min+200
## --------change t2=n.min+3 to t2=n.max for whole smc run 
## ---------(i.e t changes from n.min to T one by one)
#t2=n.max  
pb <- txtProgressBar(min = 0, max = (t2-t1+1), style = 3)
time_cost=rep(0,(t2-t1+1))

# here we remove ARI computation for real data application
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
logZ =c(logZ.py,rep(NA,(t2-t1+1)))
#here we only save beta1[ind1],beta2[ind2],...,betak[indk]
for ( t in t1 : t2){
  rl.t=(t-t1+1)  
  setTxtProgressBar(pb, rl.t)
  t0=Sys.time()
  yt=y[,t]
  out=pf2(yt,t,NP,K,pi,beta,sigma,xi,b2,nk,Syb,Syb2,Si2,Sit2)
  ##if(t==t1){
  ##  logZ[t]=log(mean(out$py)) 
  ##}else{
  ##  logZ[t]=log(mean(out$py))+ logZ[(t-1)]
  ##}
  logZ[t]=log(mean(out$py))+ logZ[(t-1)]
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
    xi.uci.smc [rl.t,k]=quantile(xi[,k]   ,0.975)
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
  if(t%%200==0){
    save.image(file=out_smc) 
  }
}
save.image(file=out_smc) 


#-------- image plot
rotate <- function(x) t(apply(x, 2, rev)) # matrix rotation
par(mfrow=c(3,4),oma=c(0.2,1.5,0.2,1.5),mar=c(0.2,0.2,0.2,0.2),
    cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
for(k in 1:K){
  y.k1=Si%*%beta.m[[k]]
  len=16
  image(rotate(matrix(y.k1,len,len,byrow = TRUE)),col = topo.colors(len*len),yaxt="n",xaxt='n', ann=FALSE)
}



##----calcualte WAIC by using the estiamtes of paramters at Time T
log.py.mean=rep(NA,t2)
log.py.var =rep(NA,t2)
lppd=rep(NA,t2)
pwaic2=rep(NA,t2)
for (t in 1:t2){
  WResult=computWeight_v2(NP,K,t,Si,Sit2,pi,beta,xi,sigma)
  py=WResult$py
  if(t==1){
    log.py.mean[t]=mean(log(py))
    log.py.var[t] = var(log(py))
    lppd[t]=log(mean(py))
  }else{
    log.py.mean[t]=mean(log(py)) + log.py.mean[(t-1)]
    log.py.var[t] = var(log(py)) + log.py.var[(t-1)]
    lppd[t]       = log(mean(py))+ lppd[(t-1)]
  }
}
waic1=2*lppd-4*log.py.mean
waic2=2*log.py.var-2*lppd

save.image(file=out_smc) 


