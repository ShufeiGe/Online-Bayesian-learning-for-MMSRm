## Simulate Y
## Simulation 1: Get Y(y1,..,yn) ;yi=(yi1,...,yimi)'
## n: number of surfaces (observations)
## mi: dimension of yi ;here m1=m2...=mn=m
x1=x2=seq(-1,1,length=sqrt(mi))
x1_x2=expand.grid(x1,x2)
X1=x1_x2$Var1
X2=x1_x2$Var2


## d: No. of NBFs basis ; 
## Create Si (nrow(Si)=mi, ncol(Si)=d))
## d=6 X 6 NBFs
## X1 X X2 =[-1,1] X [-1,1]
## c=(c1,c2,...,cd)
## ci=(ci1,ci2)'
## Create d (center) nodes

d1=6
d2=6
d=d1*d2

##  c :centers of each grid
delta_1=(max(x1)-min(x1))/(d1-1)
delta_2=(max(x2)-min(x2))/(d2-1)
c1=seq(min(x1),max(x1),by=delta_1)
c2=seq(min(x2),max(x2),by=delta_2)
c=as.matrix(t(expand.grid(c1,c2)))
row.names(c)=c("x1.center","x2.center")
rm(c1,c2)
# S1=S2=...=Sn; nrow(Si)=mi; ncol(Si)=d
# Si=matrix(0,nrow=mi,ncol=d)

Si=NBFs(X1,X2,c,delta_1,delta_2)  #return Si= mi x d covariates matrix
Si_Si=crossprod(Si,Si)            #Si_Si=  t(Si)%*%Si  


ux.1=function(x1,x2){
  z=(x1^3+x2^3+3)/sqrt(1+x1^2+x2^2) 
  return(z)}

ux.2=function(x1,x2){ 
  z=(x1^2+1+x2^2)/sqrt(x1^2+x2/4+4) 
  return(z)
}

ux.3=function(x1,x2){ 
  z=1-sin(x1^2+1)+cos(1+x2^2)/2 
  return(z)
}

ux.4=function(x1,x2){ 
  z=sin(x1*x2)
  return(z)
}

ux.5=function(x1,x2){ 
  z=cos(x1+x2)+sin(x1^2)+cos(x2^2)
  return(z)
}

ux.6=function(x1,x2){
  z=x1+x2
  return(z)
  }


b=matrix(0,nrow=d,ncol=n)
err=matrix(0,nrow=mi,ncol=n)
for (i in 1:n)
{
  b[,i]=rnorm(d,0,sqrt(0.3))
  err[,i]=rnorm(mi,0,sqrt(0.1))
}

Y=matrix(0,nrow=mi,ncol=n)

for (i in 1:n){
  Y[,i]=ux.1(X1,X2)*(Z[i]==1)+ux.2(X1,X2)*(Z[i]==2)+ux.3(X1,X2)*(Z[i]==3)+ux.4(X1,X2)*(Z[i]==4)+ux.5(X1,X2)*(Z[i]==5)+ux.6(X1,X2)*(Z[i]==6)+Si%*%b[,i]+err[,i]
}
Y.mean=matrix(0,nrow=mi,ncol=K)
if(K==1){
  Y.mean=ux.1(X1,X2)*(f_k==1)+ux.2(X1,X2)*(f_k==2)+ux.3(X1,X2)*(f_k==3)+ux.4(X1,X2)*(f_k==4)+ux.5(X1,X2)*(f_k==5)+ux.6(X1,X2)*(f_k==6)
}else{
  for (f_k in 1:K){
    Y.mean[,f_k]=ux.1(X1,X2)*(f_k==1)+ux.2(X1,X2)*(f_k==2)+ux.3(X1,X2)*(f_k==3)+ux.4(X1,X2)*(f_k==4)+ux.5(X1,X2)*(f_k==5)+ux.6(X1,X2)*(f_k==6)
  } 
}


#call algorithm1
it_max=5000
burn_in=1000
#seed=4636890
Y.all=Y
Y.all.label=Z
Y.label=Z[1:n.min]
Y=Y.all[,1:n.min]
n=n.min
## hyper parameter

alpha0=rep(1,K)
d0=u0=rep(0,d)
Sigma0=diag(d)
g0=10
h0=1
a0=10
b0=1
hyper=list(alpha0=alpha0,g0=g0,h0=h0,a0=a0,b0=b0)

xi=rep(0.01,K)    # here xi    is a vector with length k    
sigma=rep(0.01,K) # here sigma is a vector with length k

b=array(0,c(K,d,n))   
alpha=rep(1/K,K)
#get initial values for Beta and pi.K
source("par_initializationbyhrclu.R")
mod_par=list(xi= xi,sigma=sigma,pi.K=pi.K,Beta=Beta,b=b,K=K,alpha=alpha)


