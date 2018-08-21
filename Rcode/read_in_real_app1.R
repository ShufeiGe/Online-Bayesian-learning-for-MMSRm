source("NBFs.R")
Y1=read.table("zip.test")
Y2=read.table("zip.train")
Y.all.label=c(Y1[,1],Y2[,1])
set.seed(123)
Y.all=t(rbind(Y1[,2:257],Y2[,2:257]))
Y.all=scale(Y.all,center = TRUE,scale=TRUE)
##  #randomly shuffle the observations
ind.y=sample(1:9000,replace = FALSE)
Y.all=Y.all[,ind.y]
Y.all.label=Y.all.label[ind.y]
Y.all=Y.all[,1:n.max]
Y.all.label=Y.all.label[1:n.max]

 

#------------------------------
mi=16*16
x1=x2=seq(-1,1,length=sqrt(mi))
x1_x2=expand.grid(x1,x2)
X1=x1_x2$Var1
X2=x1_x2$Var2

## parameters ...
## No. of NBFs basis ; Create Si (nrow(Si)=mi, ncol(Si)=d))
## Simulation (1.1)  : d=5 X 5 NBFs
## X1 X X2 =[-10,10] X [-10,10]
## c=(c1,c2,...,cd)
## ci=(ci1,ci2)'
## Create d (center) nodes


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

#system.time({Si_Si=crossprod(Si,Si)})
#system.time({Si_Si=  t(Si)%*%Si } )

n=n.min
it_max=2000
burn_in=1000

Y=Y.all[,1:n.min]
Y.label=Y.all.label[1:n.min]
## hyper parameter

alpha0=rep(1,K)
d0=u0=rep(0,d)
Sigma0=diag(d)
g0=10
h0=1
a0=10
b0=1
hyper=list(alpha0=alpha0,g0=g0,h0=h0,a0=a0,b0=b0)

#set.seed(seed)
#pi=rdirichlet(1,alpha0)
xi=rep(0.01,K)
#rinvgamma(K,a0/10,b0*1000)           # here xi    is a vector with length k    
sigma=rep(0.01,K)
#rinvgamma(K,g0/10,h0*1000)      # here sigma is a vector with length k

#Beta=array(1,c(K,d,1))
#load('sim_v1_betaint.RData')
#Beta[1,,]=beta1
#Beta[2,,]=beta2
#Beta[3,,]=beta3
b=array(0,c(K,d,n))   
alpha=rep(1/K,K)
#pi.K=alpha
#get initial values for Beta and pi.K
source("par_initializationbyhrclu.R")
mod_par=list(xi= xi,sigma=sigma,pi.K=pi.K,Beta=Beta,b=b,K=K,alpha=alpha)

