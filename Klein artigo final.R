#A control chart for the joint monitoring of mean and variance with Klein's supplementary rule
require(pracma)
rm(list = ls())
####################################################################################################
u0=0 #In-control average
u1=1.2 #out-of-control average (defined in the previous vector)
s0=1 #in-control standard deviation 
s1=1.05 # out-of-control standard deviation (defined in the previous vector)
n=5  #Sample size
ARL0=370.4 #Target ARL0

####################################################################################################
OtiUCL <- function(U){#Optimization function used to find the UCL and LCL
  UCLxb=qnorm((1-U[1]/2),u0,s0/(n^0.5))
  LCLxb=qnorm(U[1]/2,u0,s0/(n^0.5))
  LCs2=qchisq((1-U[2]),(n-1))
  #X-bar control chart 
  pxi=pnorm(LCLxb,u0,s0/(n^0.5))
  pxs=1-pnorm(UCLxb,u0,s0/(n^0.5))
  pxc=1-pxs-pxi
  
  #S2 control chart 
  psf=1-pchisq(LCs2,(n-1))
  psc=1-psf
  #Markov chain
  size<- 15 #Size of the markov chain
  MarkovChain<- matrix(0,nrow=size,ncol=size,byrow=TRUE)
  MarkovChain[1,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[2,]<- c(pxc*psc,0,pxs*psc,pxi*psc,0,pxc*psf,0,pxs*psf,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[3,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[4,]<- c(pxc*psc,pxs*psc,0,0,pxi*psc,pxc*psf,pxs*psf,0,0,pxi*psf,0,0,0,0,0) 
  MarkovChain[5,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[6,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,0,0,0,0,0,pxc*psf,pxs*psf,0,pxi*psf,0) 
  MarkovChain[7,]<- c(pxc*psc,0,pxs*psc,pxi*psc,0,0,0,0,0,0,pxc*psf,0,pxs*psf,pxi*psf,0)   
  MarkovChain[8,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0)   
  MarkovChain[9,]<- c(pxc*psc,pxs*psc,0,0,pxi*psc,0,0,0,0,0,pxc*psf,pxs*psf,0,0,pxi*psf)   
  MarkovChain[10,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0)   
  MarkovChain[11,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[12,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[13,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[14,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  MarkovChain[15,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
  
  #The stationary distribution
  A = t(MarkovChain) - eye(15)
  A[15,] = ones(1,15)
  B = zeros(15, 1)
  B[15,1] = 1
  Solved_markov_chain = solve(A)%*%B
  
  ARL<- 1/sum(Solved_markov_chain[3,],Solved_markov_chain[5,],
              Solved_markov_chain[8,],Solved_markov_chain[10,],
              Solved_markov_chain[11,],Solved_markov_chain[12,],
              Solved_markov_chain[13,],Solved_markov_chain[14,],
              Solved_markov_chain[15,])
  
  ARLphi=(ARL-ARL0)^2
  
  return(ARLphi)
}

####################################################################################################
#Limit used in the optimize function for the x-bar control chart 
LSa=1 
LSb=1 
#initial value. Change as needed.
ivxbar=0.04
ivs2=0.04
par_optim <- nlminb(c(ivxbar,ivs2),OtiUCL,lower=c(0,0),upper =c(LSa,LSb))#
Ua=par_optim[[1]][1] #Prob Xbar
Ub=par_optim[[1]][2] #Prob S2
UCLxb=qnorm((1-Ua/2),u0,s0/(n^0.5)) #upper control limit for the x-bar
LCLxb=qnorm(Ua/2,u0,s0/(n^0.5)) #lower control limit for the x-bar
LCs2=qchisq((1-Ub),(n-1)) #control limit for the s2

#X-bar control chart
pxi=pnorm(LCLxb,u1,s1/(n^0.5))
pxs=1-pnorm(UCLxb,u1,s1/(n^0.5))
pxc=1-pxs-pxi

#S2 control chart
psf=1-pchisq(LCs2*(s0^2/s1^2),(n-1))
psc=1-psf

####################################################################################################
#Now that we have the probabilities for this specific case, we solve it through a Markov Chain again
size<- 15 #Size of the markov chain
MarkovChain<- matrix(0,nrow=size,ncol=size,byrow=TRUE)

MarkovChain[1,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[2,]<- c(pxc*psc,0,pxs*psc,pxi*psc,0,pxc*psf,0,pxs*psf,pxi*psf,0,0,0,0,0,0) 
MarkovChain[3,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[4,]<- c(pxc*psc,pxs*psc,0,0,pxi*psc,pxc*psf,pxs*psf,0,0,pxi*psf,0,0,0,0,0) 
MarkovChain[5,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[6,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,0,0,0,0,0,pxc*psf,pxs*psf,0,pxi*psf,0) 
MarkovChain[7,]<- c(pxc*psc,0,pxs*psc,pxi*psc,0,0,0,0,0,0,pxc*psf,0,pxs*psf,pxi*psf,0)   
MarkovChain[8,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0)   
MarkovChain[9,]<- c(pxc*psc,pxs*psc,0,0,pxi*psc,0,0,0,0,0,pxc*psf,pxs*psf,0,0,pxi*psf)   
MarkovChain[10,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0)   
MarkovChain[11,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[12,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[13,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[14,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 
MarkovChain[15,]<- c(pxc*psc,pxs*psc,0,pxi*psc,0,pxc*psf,pxs*psf,0,pxi*psf,0,0,0,0,0,0) 

#The stationary distribution
A = t(MarkovChain) - eye(15)
A[15,] = ones(1,15)
B = zeros(15, 1)
B[15,1] = 1
Solved_markov_chain = solve(A)%*%B

ARL1<- 1/sum(Solved_markov_chain[3,],Solved_markov_chain[5,],
             Solved_markov_chain[8,],Solved_markov_chain[10,],
             Solved_markov_chain[11,],Solved_markov_chain[12,],
             Solved_markov_chain[13,],Solved_markov_chain[14,],
             Solved_markov_chain[15,])


options(digits=5)
cat('UCLxb=',UCLxb,"\n")
cat('LCLxb=',LCLxb,"\n")
cat('LCs2=',LCs2,"\n")
cat('ARL1=',ARL1,"\n")