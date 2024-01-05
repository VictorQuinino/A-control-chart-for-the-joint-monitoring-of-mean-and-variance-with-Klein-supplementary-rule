#A control chart for the joint monitoring of mean and variance 
rm(list = ls())
##Markov chain

u0=0 #In-control average
u1=2 #out-of-control average (defined in the previous vector)
s0=1 #in-control standard deviation 
s1=0.25 # out-of-control standard deviation (defined in the previous vector)
n=5  #Sample size
ARL0=250 #Target ARL0


OtiUCL <- function(alfa){
  alfag=1-(1-alfa[1])*(1-alfa[2])
  ARL=1/alfag
  ARLphi=(ARL-ARL0)^2
  return(ARLphi)
}
#upper limit to be used in the optimize function. Change as needed.
LSa=1 
LSb=1
#initial value. Change as needed.
ivxbar=0.04
ivs2=0.04
par_optim <- nlminb(c(ivxbar,ivs2),OtiUCL,lower=c(0,0),upper =c(LSa,LSb))
alfaa=par_optim[[1]][1] #Prob Xbar
alfab=par_optim[[1]][2] #Prob S2

#ARL1
LSCxb=qnorm((1-alfaa/2),u0,s0/(n^0.5))
LICxb=qnorm(alfaa/2,u0,s0/(n^0.5))
LSCqui=qchisq((1-alfab/2),(n-1))
LICqui=qchisq(alfab/2,(n-1))

pc=pnorm(LICxb,u1,s1/(n^0.5))
pa=1-pnorm(LSCxb,u1,s1/(n^0.5))
pb=1-pa-pc

pas=1-pchisq(LSCqui*(s0^2/s1^2),(n-1))
pcs=pchisq(LICqui*(s0^2/s1^2),(n-1))
pbs=1-pas-pcs


ARL1<-1/(1-pb*pbs)

cat('LSCxb=',LSCxb,"\n")
cat('LICxb=',LICxb,"\n")
cat('LSCqui=',LSCqui,"\n")
cat('Probxbar=',alfaa,"\n")
cat('Probs2=',alfab,"\n")
cat('ARL1=',ARL1,"\n")