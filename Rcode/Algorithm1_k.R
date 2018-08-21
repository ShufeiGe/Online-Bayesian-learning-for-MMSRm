##Algorithm1 
## input : 
##       it_max     - max # of iterations
##       seed       - sample seed
##       n          - number of yi (i=1,2,...,n)
##       mi         - length of yi (length(yi)=mi)
##       d          - d=d1*d2 
##                  - d1:number of rows of nodes
##                  - d2:number of cols of nodes
##       hyper parameters : d0,u0,Sigma0,g0,h0,a0,b0
##  d0=u0=rep(0,25)
##  Sigma0=diag(25)
##  g0=50
##  h0=0.01
##  a0=50
##  b0=0.01
##  hyper=list(g0= g0,h0=h0,a0=a0,b0=b0)
##  model parameters
##   set.seed(seed)
##   xi0=rinvgamma(1,a0,b0)
##   sigma0=rinvgamma(1,g0,h0)
##   Beta0=mvrnorm(1,u0,Sigma0)
##   Sigma_bi0=diag(25)*xi0
##   bi0=mvrnorm(1,d0,Sigma_bi0)
##   b0=matrix(rep(bi0,times=n),nrow=d,ncol=n)
##   mod_par=list( xi0= xi0,sigma0=sigma0,Beta0=Beta0,b=b0)


 # Algorithm1.K=function(hyper,mod_par,Y,Si,Si_Si,n,mi,d,it_max,seed){

          t0 <- proc.time()
      
          g0=hyper$g0
          h0=hyper$h0
          a0=hyper$a0
          b0=hyper$b0
     
          xi    = mod_par$xi
          sigma = mod_par$sigma
          Beta  = mod_par$Beta
          pi.K  = mod_par$pi.K
          b     = mod_par$b 
          alpha = mod_par$alpha
          
          result_xi    =  matrix(NA,nrow=it_max,ncol=K)
          result_sigma =  matrix(NA,nrow=it_max,ncol=K)
          result_Beta  =  array(NA,c(K,it_max,d,1))
          result_piK   =  matrix(NA,nrow=it_max,ncol=K)
          labZ         =  rep(NA,n)
          result_Z     =  array(NA,c(it_max,n))
          result_b     =  array(0,c(K,d,n))
          
   set.seed(seed)
   total <- it_max
   # create progress bar
   pb <- txtProgressBar(min = 0, max = total, style = 3)
   Z=matrix(0,ncol=K,nrow=n)
   tau1=tau2=rep(0,K)
 
   
for (t in 1:it_max){
  # Sys.sleep(0.001)
  # update progress bar
  setTxtProgressBar(pb, t)

  ## STEP0. SAMPLE Zi and mixing proportions pi.k
  for (i in 1:n){
     for (k in 1:K){
       mu.star=Si%*%(Beta[k,,]+b[k,,i])
       Y.err=Y[,i]-mu.star
       exp.term=t(Y.err)%*%Y.err/sigma[k]/2
       tau1[k]=log(pi.K[k])-mi/2*log(sigma[k])-exp.term
       rm(list=c("mu.star","Y.err","exp.term"))
      }
     for (k in 1:K){
     tau2[k]=1/sum(exp(tau1-tau1[k]))
     }
    tau2
     Z[i,]=rmultinom(1, size = 1, prob =tau2)
     labZ[i]=which(Z[i,]==1)
  }  
  n.K=colSums(Z)
  pi.K=rdirichlet(1, (alpha+n.K))
  result_piK[t,]=pi.K
  result_Z[t,]=labZ 
    
  for (k in 1:K)
  {
    # Step1.  Sample the random-effects variance: xi(t)~IG(shape, rate)
    # invgamma in pacakge "LaplacesDemon"  IG(shape, scale)
    shape=a0+(n*d)/2
    rate=b0+sum(b[k,,]*b[k,,])/2         #b(t-1)=(b1,...,bn)
    xi[k]=LaplacesDemon::rinvgamma(1, shape=shape, scale =rate) 
   # Step2.  Sample the noise variance  :sigma(t)~IG(shape, rate)
   shape=g0+(n.K[k]*mi)/2
   #Y.err=(Y-matrix(Si%*%Beta[k,,],nrow=mi,ncol=n)-Si%*%b[k,,])*matrix(Z[,k],ncol=n,nrow=mi,byrow=TRUE)
   #rate=h0+sum(Y.err*Y.err)/2
   Y.err=(Y-matrix(Si%*%Beta[k,,],nrow=mi,ncol=n,byrow=FALSE)-Si%*%b[k,,])
   Y.err.2=colSums(Y.err^2)
   rate=h0+sum(Y.err.2*Z[,k])/2
   sigma[k]=LaplacesDemon::rinvgamma(1, shape=shape, scale =rate)   
   #rm(list=c("Y.err"))
   rm(list=c("Y.err","Y.err.2"))
   
  # Step3.  Sample the fixed-effects coefficients vector   
  #                                     Beta(t)~N(v0(t),V0(t))
  V0=solve(solve(Sigma0)+n.K[k]*Si_Si/sigma[k])   # sigma(t)   
  v0=V0%*%((t(Si)%*%rowSums(Y*matrix(Z[,k],ncol=n,nrow=mi,byrow=TRUE))-Si_Si%*%rowSums(b[k,,]*matrix(Z[,k],ncol=n,nrow=d,byrow=TRUE)))/sigma[k]+solve(Sigma0)%*%u0)   #b=b(t-1);
  Beta[k,,]=mvrnorm(1,v0,V0) 
  result_Beta[k,t,,]=Beta[k,,]
  # Step4. Sample the random-effects coefficients vector bi(t) ~ N(v1(t),V1(t))
  # b(t)=(b1(t),...,bn(t)) ;shorthand b=(b1,..bt)
  for (i in 1:n){
    V1=solve(Si_Si/sigma[k]*Z[i,k]+diag(d)/xi[k])
    v1=Z[i,k]*V1%*%(t(Si)%*%Y[,i]-Si_Si%*%Beta[k,,])/sigma[k]
   b[k,,i]=mvrnorm(1,v1,V1)   
  }
  den =it_max-burn_in  # denominator
#  den=5
  if (t>burn_in) {
    result_b[k,,]=result_b[k,,]+b[k,,]/den 
   }
  #result_b[k,t,,]=b[k,,]
  }
 result_xi[t,]=xi
 result_sigma[t,]=sigma
 updateRdata=t%%1000
 if (updateRdata==0){
   save.image(out_mcmc)  
 }
   } 
close(pb) 
t0 <- proc.time()-t0 
save.image(out_mcmc)
#}
