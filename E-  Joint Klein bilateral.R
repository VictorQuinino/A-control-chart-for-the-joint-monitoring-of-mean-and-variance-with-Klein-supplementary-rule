#A control chart for the joint monitoring of mean and variance with Klein's supplementary rule

require(pracma)
rm(list = ls())

cat('N','\t', 'u1','\t','s1','\t','ARL1','\n')

m1f<- c(0.00,0.25,0.50,1,2) #Vector containing different out-of-control averages
m1fs<- size(m1f)
s1f<- c(0.25,0.50,1,1.50,2)#Vector containing different out-of-control standard deviations
s1fs<- size(s1f)
nf<-  c(5) #Vector containing different sample sizes
nfs<- size(nf)

####################################################################################################
for(k in 1:nfs[2]){
  for(i in 1:m1fs[2]){
    for(j in 1:s1fs[2]){
      u0=0 #In-control average
      u1=m1f[i] #out-of-control average (defined in the previous vector)
      s0=1 #in-control standard deviation 
      s1=s1f[j]# out-of-control standard deviation (defined in the previous vector)
      n=nf[k] #Sample size
      ARL0=250 #Target ARL0
      
      ####################################################################################################
      OtiUCL <- function(U){#Optimization function used to find the UCL and LCL
        UCLxb=qnorm((1-U[1]/2),u0,s0/(n^0.5))
        LCLxb=qnorm(U[1]/2,u0,s0/(n^0.5))
        LCs2a=qchisq((1-U[2]),(n-1))
        LCs2b=qchisq(U[3],(n-1))
        
        #X-bar control chart 
        pxl=pnorm(LCLxb,u0,s0/(n^0.5))
        pxu=1-pnorm(UCLxb,u0,s0/(n^0.5))
        pxc=1-pxu-pxl
        
        #S2 control chart 
        psu=1-pchisq(LCs2a,(n-1))
        psl=pchisq(LCs2b,(n-1))
        psc=1-psu-psl
        
        #Markov chain
        size<- 25 #Size of the markov chain
        MarkovChain<- matrix(0,nrow=size,ncol=size,byrow=TRUE)
        
        MarkovChain[1,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[2,]<- c(pxc*psc,0,pxc*psl,pxc*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,pxl*psc,0,pxl*psl,pxl*psu,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[3,]<- c(pxc*psc,pxc*psu,0,0,pxc*psl,pxu*psc,pxu*psu,0,0,pxu*psl,pxl*psc,pxl*psu,0,0,pxl*psl,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[4,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[5,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[6,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,0,0,0,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,0,0,0,0,0)
        MarkovChain[7,]<- c(pxc*psc,0,pxc*psl,pxc*psu,0,0,0,0,0,0,pxl*psc,0,pxl*psl,pxl*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,0,0,0,0,0)
        MarkovChain[8,]<- c(pxc*psc,pxc*psu,0,0,pxc*psl,0,0,0,0,0,pxl*psc,pxl*psu,0,0,pxl*psl,pxu*psc,pxu*psu,0,0,pxu*psl,0,0,0,0,0)
        MarkovChain[9,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[10,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[11,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,0,0,0,0,0,0,0,0,0,0,pxl*psc,pxl*psu,pxl*psl,0,0)
        MarkovChain[12,]<- c(pxc*psc,0,pxc*psl,pxc*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,0,0,0,0,0,0,0,0,0,0,pxl*psc,0,pxl*psl,pxl*psu,0)
        MarkovChain[13,]<- c(pxc*psc,pxc*psu,0,0,pxc*psl,pxu*psc,pxu*psu,0,0,pxu*psl,0,0,0,0,0,0,0,0,0,0,pxl*psc,pxl*psu,0,0,pxl*psl)
        MarkovChain[14,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[15,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[16,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[17,]<- c(pxc*psc,0,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[18,]<- c(pxc*psc,pxc*psu,0,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[19,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[20,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[21,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[22,]<- c(pxc*psc,0,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[23,]<- c(pxc*psc,pxc*psu,0,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[24,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
        MarkovChain[25,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0) 
        
        
        
        A = t(MarkovChain) - eye(25)
        A[25,] = ones(1,25)
        B = zeros(25, 1)
        B[25,1] = 1
        Solved_markov_chain = solve(A)%*%B
        
        ARL<- 1/sum(Solved_markov_chain[4,],Solved_markov_chain[5,],
                    Solved_markov_chain[9,],Solved_markov_chain[10,],
                    Solved_markov_chain[14,],Solved_markov_chain[15,],
                    Solved_markov_chain[16,],Solved_markov_chain[17,],
                    Solved_markov_chain[18,],Solved_markov_chain[19,],
                    Solved_markov_chain[20,],Solved_markov_chain[21,],
                    Solved_markov_chain[22,],Solved_markov_chain[23,],
                    Solved_markov_chain[24,],Solved_markov_chain[25,])
        
        ARLphi=(ARL-ARL0)^2
        #ARLphi=ARL
        return(ARLphi)
      }
      
      ####################################################################################################
      LSa=0.06462 #Limit used in the optimize function for the x-bar control chart 
      LSb=0.11456 #Limit used in the optimize function for the s2 control chart
      par_optim <- nlminb(c(0.04,0.04,0.04),OtiUCL,lower=c(0,0,0),upper =c(LSa,LSb,LSb))#
      Ua=par_optim[[1]][1] #Prob Xbar
      Ub=par_optim[[1]][2] #Prob S2
      Uc=par_optim[[1]][3]
      
      UCLxb=qnorm((1-Ua/2),u0,s0/(n^0.5)) #upper control limit for the x-bar
      LCLxb=qnorm(Ua/2,u0,s0/(n^0.5)) #lower control limit for the x-bar
      LCs2a=qchisq((1-Ub),(n-1)) #control limit for the s2
      LCs2b=qchisq(Uc,(n-1)) #control limit for the s2
      
      #X-bar control chart
      pxl<-pnorm(LCLxb,u1,s1/(n^0.5))
      pxu=1-pnorm(UCLxb,u1,s1/(n^0.5))
      pxc=1-pxu-pxl
      
      #S2 control chart
      psu=1-pchisq(LCs2a*(s0^2/s1^2),(n-1))
      psl=pchisq(LCs2b*(s0^2/s1^2),(n-1))
      psc=1-psu-psl
      ####################################################################################################
      #Now that we have the probabilities for this specific case, we solve it through a Markov Chain again
      size<- 25 #Size of the markov chain
      MarkovChain<- matrix(0,nrow=size,ncol=size,byrow=TRUE)
      
      MarkovChain[1,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[2,]<- c(pxc*psc,0,pxc*psl,pxc*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,pxl*psc,0,pxl*psl,pxl*psu,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[3,]<- c(pxc*psc,pxc*psu,0,0,pxc*psl,pxu*psc,pxu*psu,0,0,pxu*psl,pxl*psc,pxl*psu,0,0,pxl*psl,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[4,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[5,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[6,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,0,0,0,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,0,0,0,0,0)
      MarkovChain[7,]<- c(pxc*psc,0,pxc*psl,pxc*psu,0,0,0,0,0,0,pxl*psc,0,pxl*psl,pxl*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,0,0,0,0,0)
      MarkovChain[8,]<- c(pxc*psc,pxc*psu,0,0,pxc*psl,0,0,0,0,0,pxl*psc,pxl*psu,0,0,pxl*psl,pxu*psc,pxu*psu,0,0,pxu*psl,0,0,0,0,0)
      MarkovChain[9,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[10,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[11,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,0,0,0,0,0,0,0,0,0,0,pxl*psc,pxl*psu,pxl*psl,0,0)
      MarkovChain[12,]<- c(pxc*psc,0,pxc*psl,pxc*psu,0,pxu*psc,0,pxu*psl,pxu*psu,0,0,0,0,0,0,0,0,0,0,0,pxl*psc,0,pxl*psl,pxl*psu,0)
      MarkovChain[13,]<- c(pxc*psc,pxc*psu,0,0,pxc*psl,pxu*psc,pxu*psu,0,0,pxu*psl,0,0,0,0,0,0,0,0,0,0,pxl*psc,pxl*psu,0,0,pxl*psl)
      MarkovChain[14,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[15,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[16,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[17,]<- c(pxc*psc,0,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[18,]<- c(pxc*psc,pxc*psu,0,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[19,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[20,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[21,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[22,]<- c(pxc*psc,0,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[23,]<- c(pxc*psc,pxc*psu,0,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[24,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      MarkovChain[25,]<- c(pxc*psc,pxc*psu,pxc*psl,0,0,pxu*psc,pxu*psu,pxu*psl,0,0,pxl*psc,pxl*psu,pxl*psl,0,0,0,0,0,0,0,0,0,0,0,0)
      
      #Solved_markov_chain<-MarkovChain%^%100000
      
      a<- replicate(25,0)
      for(l in 1:size){
        
        a[l]<- sum(MarkovChain[i,])
      }
      
      
      A = t(MarkovChain) - eye(25)
      A[25,] = ones(1,25)
      B = zeros(25, 1)
      B[25,1] = 1
      Solved_markov_chain = solve(A)%*%B
      
      ARL1<- 1/sum(Solved_markov_chain[4,],Solved_markov_chain[5,],
                   Solved_markov_chain[9,],Solved_markov_chain[10,],
                   Solved_markov_chain[14,],Solved_markov_chain[15,],
                   Solved_markov_chain[16,],Solved_markov_chain[17,],
                   Solved_markov_chain[18,],Solved_markov_chain[19,],
                   Solved_markov_chain[20,],Solved_markov_chain[21,],
                   Solved_markov_chain[22,],Solved_markov_chain[23,],
                   Solved_markov_chain[24,],Solved_markov_chain[25,])
      
      
      options(digits=5)
      cat(n,'\t',u1,'\t',s1,'\t',ARL1,'\n')
      
    }}}