# simulation study example :
#          ----surface fitting via mcmc(Gibbs sampleing) algorithm----
#          ----see paper section '5.1 surface fitting'
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
n = 100
n.min=100
n.max=100
K =1

f_k=1   # f_k=1,...6, indicates which function is used to do simulation
pi.true=rep((1/K),K)
Z=rep(f_k,n)
#Z <- sample(1:K, size = n, replace = TRUE, prob = pi.true)
mi=12*12
source("NBFs.R")
source("simulator.R")

# setp 2:  surface fitting via mcmc algorithm
out_mcmc=paste("mcmc_fitting_k_",f_k,".RData",sep="") # output .R file name
source("Algorithm1_k.R")
sigma.mean=apply(result_sigma,2,mean)
xi.mean=apply(result_xi,2,mean)
pi.K.mean=apply(result_piK,2,mean)
#--------------------------------------------------------
Beta.mean=list()
b.mean=b
for (k in 1:K){
  Beta.mean[[k]]=apply(result_Beta[k,(burn_in+1):it_max,,],2,mean)
}

Y.fitted=Si%*%Beta.mean[[1]]
plot_ly(z =matrix(c(Y.fitted),12,12)) %>% add_surface()
plot_ly(z =matrix(c(Y.mean),12,12)) %>% add_surface()

#plot histogram
par(mfrow=c(2,K))
# histogram of xi^2
hist(result_xi[-c(1:1000)], breaks=400,xlab = "", ylab=" " ,main=" ")
abline(v=0.3,lty=1,col="red")
abline(v=mean(c(result_xi[-c(1:1000)])),lty=1,col="blue")
abline(v=quantile(c(result_xi[-c(1:1000)]),probs=c(0.025,0.975)),lty=2,col="blue")
# histogram of sigma^2
hist(result_sigma[-c(1:1000)], breaks=400,xlab = "", ylab=" " ,main=" ")
abline(v=0.1,lty=1,col="red")
abline(v=mean(c(result_sigma[-c(1:1000)])),lty=1,col="blue")
abline(v=quantile(c(result_sigma[-c(1:1000)]),probs=c(0.025,0.975)),lty=2,col="blue")


#output x1,x2, y.fitted, y.mean, then we can use the data to make plots in Matlab.
file1=paste("Y_true_k",f_k,".csv",sep="")
file2=paste("Y_fitted_k",f_k,".csv",sep="")
write.csv(matrix(c(Y.fitted),12,12),file = file1)
write.csv(matrix(c(Y.mean),12,12),file = file2)
file1=paste("X1",".csv",sep="")
file2=paste("X2",".csv",sep="")
write.csv(matrix(c(X1),12,12),file = file1)
write.csv(matrix(c(X2),12,12),file = file2)


