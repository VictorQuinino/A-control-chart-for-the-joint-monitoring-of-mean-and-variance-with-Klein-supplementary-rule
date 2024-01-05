library(pracma)
clear()
tic()
u0=0
sd0=1
lambda1=0.1
lambda2=0.1
k1=2.53
k2=2.58
n=5
LCL1=u0-k1*((lambda1/(2-lambda1))^0.5)*(sd0/(n^0.5))
UCL1=u0+k1*((lambda1/(2-lambda1))^0.5)*(sd0/(n^0.5))
c=log(sd0^2)-1/(n-1)-1/(3*(n-1)^2)+2/(15*(n-1)^4)
d=2/(n-1)+2/((n-1)^2)+4/(3*(n-1)^3)-16/(15*(n-1)^5)
LCL2=c-k2*(d*lambda2/(2-lambda2))^0.5
UCL2=c+k2*(d*lambda2/(2-lambda2))^0.5

corridas=500000

U<- c(0.00,0.05,0.1,0.25,0.50,1,2) #Vector containing different out-of-control averages
Un<- size(U)
Un=Un[2]
S<- c(0.25,0.50,0.95,1,1.05,1.50,2)#Vector containing different out-of-control standard deviations
Sn<- size(S)
Sn=Sn[2]


  for(j1 in 1:Un){
    for(j2 in 1:Sn){
u1=U[j1] #out-of-control average (defined in the previous vector)
sd1=S[j2]# out-of-control standard deviation (defined in the previous vector)


GG<-matrix(0, corridas, 1)

for(i in 1:corridas){
  #print(i)
  ua=0
  sa=c
  s1=0
  s2=0
  s1a=0
  s2a=0
  
  D<-c()
  
  while(s1<2 & s2<2 & s1a<2 & s2a<2){
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
  GG[i,1]=sum(D)
  
  
}

ARL<-mean(GG[,1])
options(digits=5)
cat(n,'\t',u1,'\t',sd1,'\t',ARL,'\n')
    

      }
  
    }
  
 



toc()

