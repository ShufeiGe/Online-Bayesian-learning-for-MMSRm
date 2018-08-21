## Nodal basis functions 
piecewise=function(x1,x2,c1,c2,delta_1,delta_2){
  c1.0=c1
  c1.1=c1-delta_1
  c1.2=c1+delta_1
  
  c2.0=c2
  c2.1=c2-delta_2
  c2.2=c2+delta_2
  c2.3=delta_2*(x1-c1)/delta_1+c2
  c2.4=delta_2*(x1-c1)/delta_1+c2-delta_2
  c2.5=delta_2*(x1-c1)/delta_1+c2+delta_2
  
  term1=(c1.0<x1 && x1<=c1.2)*(c2.3<=x2 && x2<c2.2)
  term2=(c1.0<x1 && x1<=c1.2)*(c2.0<=x2 && x2<c2.3)
  term3=(c1.0<x1 && x1<=c1.2)*(c2.4<=x2 && x2<c2.0)
  
  term4=(c1.1<x1 && x1<=c1.0)*(c2.1<=x2 && x2<=c2.3)
  term5=(c1.1<x1 && x1<=c1.0)*(c2.3<x2 && x2<=c2.0)
  term6=(c1.1<x1 && x1<=c1.0)*(c2.0<x2 && x2<=c2.5)
  
  y=term1*(1+(c2-x2)/delta_2)+
    term2*(1+(c1-x1)/delta_1)+
    term3*(1+(c1-x1)/delta_1+(x2-c2)/delta_2)+
    term4*(1+(x2-c2)/delta_2)+
    term5*(1+(x1-c1)/delta_1)+
    term6*(1+(x1-c1)/delta_1+(c2-x2)/delta_2)
 return(y) 
}

## Input: vectors x1,x2
## length(x1)=length(x2)=m 
## d :number of regularly spaced points
## c=(c1,...,cd)
## ci=(ci1,ci2)'
NBFs=function(x1,x2,c,delta_1,delta_2){
  
  m=length(x1)
  d=ncol(c)
  S=matrix(0,nrow=m,ncol=d)
  for ( i in 1:m){          
    for (j in 1:d){         #c[,j]
      c1=c[1,j]
      c2=c[2,j]
  S[i,j]=piecewise(x1[i],x2[i],c1,c2,delta_1,delta_2)
    }
  }
  return(S)
}

