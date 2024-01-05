library(pracma)
clear()
tic()

# Definition of parameters
# Target ARL0=370.4
# In-control
u0 = 0
sd0 = 1

# Out-of-control - Table 2
u1 = 0.5
sd1 = 1.5
n = 5

# Limits obtained by experimentation to reach target ARL0
# LCL1=-0.98188; UCL1=0.98188; UCL2=8.3347 #n=4
LCL1 = -0.87822; UCL1 = 0.87822; UCL2 = 10.051 #n=5
# LCL1=-0.8017; UCL1=0.8017; UCL2=11.671 #n=6
# LCL1=-0.7422; UCL1=0.7422; UCL2=13.227 #n=7

runs = 1000000 # number of simulations for Monte Carlo
Result <- matrix(0, runs, 1)

for (i in 1:runs) {
  s1 = 0
  s2 = 0
  s1a = 0
  D <- c()
  
  while (s1 < 2 & s2 < 2 & s1a < 2) { # klein rule
    RR <- rnorm(n, u1, sd1)
    R = mean(RR)
    V = var(RR)
    T = R
    W = V * (n - 1) / sd0
    
    if (T > UCL1 & W > UCL2) {
      D <- rbind(D, 1)
      s1 = s1 + 1
      s2 = 0
      s1a = s1a + 1
    }
    
    if (T < LCL1 & W > UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = s2 + 1
      s1a = s1a + 1
    }
    
    if (T >= LCL1 & T <= UCL1 & W > UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = 0
      s1a = s1a + 1
    }
    
    if (T > UCL1 & W < UCL2) {
      D <- rbind(D, 1)
      s1 = s1 + 1
      s2 = 0
      s1a = 0
    }
    
    if (T < LCL1 & W < UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = s2 + 1
      s1a = 0
    }
    
    if (T >= LCL1 & T <= UCL1 & W < UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = 0
      s1a = 0
    }
    
    ua = T
    sa = W
  }
  Result[i, 1] = sum(D)
}

ARL <- mean(Result[, 1])
cat("ARL1=", ARL, "\n")
toc()
