library(pracma)
clear()
tic()
# Definition of parameters
#Target ARL0=250. 
#In-control
u0=0
sd0=1
#out-of-control 
u1=0
sd1=1
lambda1=0.1
lambda2=0.1
n=5
#Limits obtained by experimentation to reach target ARL0
LCL1=-0.25957
UCL1=0.25957
LCL2=-0.7456
UCL2=0.20497
c=log(sd0^2)-1/(n-1)-1/(3*(n-1)^2)+2/(15*(n-1)^4) #starting value
runs=10000
Result<-matrix(0, runs, 1)

for(i in 1:runs){
  print(i)
  ua=0 #starting value
  sa=c #starting value
  s1=0
  s2=0
  s1a=0
  s2a=0
  D<-c()
  while(s1<2 & s2<2 & s1a<2 & s2a<2){ #klein rule
    RR<-rnorm(n,u1,sd1)
    R=mean(RR)
    V=var(RR)
    T=lambda1*R+(1-lambda1)*ua
    W=lambda2*log(V)+(1-lambda2)*sa
    
    if (T>UCL1 & W>UCL2){
      D<-rbind(D,1)
      s1=s1+1
      s2=0
      s1a=s1a+1
      s2a=0
    }
    
    if (T>UCL1 & W<LCL2){
      D<-rbind(D,1)
      s1=s1+1
      s2=0
      s1a=0
      s2a=s2a+1
    }
    if (T<LCL1 & W>UCL2){
      D<-rbind(D,1)
      s1=0
      s2=s2+1
      s1a=s1a+1
      s2a=0
    }
    
    if (T<LCL1 & W<LCL2){
      D<-rbind(D,1)
      s1=0
      s2=s2+1
      s1a=0
      s2a=s2a+1
    }
    
    if (T>=LCL1 & T<= UCL1 & W>=LCL2 & W<= UCL2 ){
      D<-rbind(D,1)
      s1=0
      s2=0
      s1a=0
      s2a=0
    }
    
    
    if (T>=LCL1 & T<= UCL1 & W<=LCL2){
      D<-rbind(D,1)
      s1=0
      s2=0
      s1a=0
      s2a=s2a+1
    }
    
    if (T>=LCL1 & T<= UCL1 & W>=UCL2){
      D<-rbind(D,1)
      s1=0
      s2=0
      s1a=s1a+1
      s2a=0
    }
    
    if (W>=LCL2 & W<= UCL2 & T>=UCL1){
      D<-rbind(D,1)
      s1=s1+1
      s2=0
      s1a=0
      s2a=0
    }
    
    if (W>=LCL2 & W<= UCL2 & T<=LCL1){
      D<-rbind(D,1)
      s1=0
      s2=s2+1
      s1a=0
      s2a=0
    }
    
    ua=T
    sa=W
  }
  Result[i,1]=sum(D)
  
  
}

ARL<-mean(Result[,1])
cat("ARL1=",ARL,"\n")
toc()