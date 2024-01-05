library(pracma)
clear()
tic()

# Definition of parameters
# Target ARL0=250. 
# In-control
u0 = 0
sd0 = 1

# Out-of-control - Table 3. 
u1 = 1
sd1 = 1
n = 5

# Limits obtained by experimentation to reach target ARL0
LCL1 = -0.91685; UCL1 = 0.91685; LCL2 = 0.63447; UCL2 = 9.9741  # n = 5
runs = 1000000
Result <- matrix(0, runs, 1)

for (i in 1:runs) {
  
  s1 = 0
  s2 = 0
  s1a = 0
  s2a = 0
  D <- c()
  
  while (s1 < 2 & s2 < 2 & s1a < 2 & s2a < 2) {  # Klein rule
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
      s2a = 0
    }
    
    if (T > UCL1 & W < LCL2) {
      D <- rbind(D, 1)
      s1 = s1 + 1
      s2 = 0
      s1a = 0
      s2a = s2a + 1
    }
    
    if (T < LCL1 & W > UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = s2 + 1
      s1a = s1a + 1
      s2a = 0
    }
    
    if (T < LCL1 & W < LCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = s2 + 1
      s1a = 0
      s2a = s2a + 1
    }
    
    if (T >= LCL1 & T <= UCL1 & W >= LCL2 & W <= UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = 0
      s1a = 0
      s2a = 0
    }
    
    if (T >= LCL1 & T <= UCL1 & W <= LCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = 0
      s1a = 0
      s2a = s2a + 1
    }
    
    if (T >= LCL1 & T <= UCL1 & W >= UCL2) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = 0
      s1a = s1a + 1
      s2a = 0
    }
    
    if (W >= LCL2 & W <= UCL2 & T >= UCL1) {
      D <- rbind(D, 1)
      s1 = s1 + 1
      s2 = 0
      s1a = 0
      s2a = 0
    }
    
    if (W >= LCL2 & W <= UCL2 & T <= LCL1) {
      D <- rbind(D, 1)
      s1 = 0
      s2 = s2 + 1
      s1a = 0
      s2a = 0
    }
    
    ua = T
    sa = W
  }
  Result[i, 1] = sum(D)
  
}

ARL <- mean(Result[, 1])
cat("ARL1=", ARL, "\n")

toc()
