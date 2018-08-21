# simulation study example :
#          ----surface clustering via mcmc(Gibbs sampleing) algorithm----
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
n.min=200   # no. of obs gonna be used in this simulation.
#n.max=10000
K =3  #clustering exampler when K=3 
#f_k=1
pi.true=rep((1/K),K)
Z <- sample(1:K, size = n, replace = TRUE, prob = pi.true)
mi=12*12
source("NBFs.R")
source("simulator.R")

# setp 2:  surface clustering via mcmc algorithm
out_mcmc=paste("mcmc_clustering_K_",K,".RData",sep="") # output .R file name
source("Algorithm1_k.R")
sigma.mean=apply(result_sigma,2,mean)
xi.mean=apply(result_xi,2,mean)
pi.K.mean=apply(result_piK,2,mean)
#--------------------------------------------------------
Beta.mean=list()
b.mean=b
Y.fitted=matrix(0,nrow=mi,ncol=K)
for (k in 1:K){
  Beta.mean[[k]]=apply(result_Beta[k,(burn_in+1):it_max,,],2,mean)
  Y.fitted[,k]=Si%*%Beta.mean[[k]]
}

#plot histogram
par(mfrow=c(2,K))
for(k in 1:K){
  # histogram of xi^2
  hist(result_xi[-c(1:1000),k], breaks=400,xlab = "", ylab=" " ,main=" ")
  abline(v=0.3,lty=1,col="red")
  abline(v=mean(c(result_xi[-c(1:1000),k])),lty=1,col="blue")
  abline(v=quantile(c(result_xi[-c(1:1000),k]),probs=c(0.025,0.975)),lty=2,col="blue")
  # histogram of sigma^2
  hist(result_sigma[-c(1:1000),k], breaks=400,xlab = "", ylab=" " ,main=" ")
  abline(v=0.1,lty=1,col="red")
  abline(v=mean(c(result_sigma[-c(1:1000),k])),lty=1,col="blue")
  abline(v=quantile(c(result_sigma[-c(1:1000),k]),probs=c(0.025,0.975)),lty=2,col="blue")
}
#output x1,x2, y.fitted, y.mean, then we can use the data to make plots in Matlab.
for(f_k in 1:K){
  file1=paste("mcmc_clustering_Y_true_k",f_k,".csv",sep="")
  file2=paste("mcmc_clustering_Y_fitted_k",f_k,".csv",sep="")
  write.csv(matrix(c(Y.fitted[,k]),12,12),file = file1)
  write.csv(matrix(c(Y.mean[,k]),12,12),file = file2)  
}

file1=paste("X1",".csv",sep="")
file2=paste("X2",".csv",sep="")
write.csv(matrix(c(X1),12,12),file = file1)
write.csv(matrix(c(X2),12,12),file = file2)


