setwd("~/github/domestication_selection_model")
library(dplyr)

T_ <- 5
phi_ <-0.9
N0 <- 100
q1 <- 0.1
E1 <- rep(1,T_)
q2 <- 0.05
E2 <- rep(1,T_)

phi <- rep(0.9,T_) # survival rate


plot(N_t)



p1 <- q1*E1
p2 <- q2*E2

sigma_12 <- 0.05
sigma_21 <- 0.05
A <- function(t){
  cbind(c(1,0,0,0,0),
        
        c(1-phi[t], 
          phi[t]*(1-sigma_12)*(1-p1[t]), 
          phi[t]*sigma_12*(1-p2[t]), 
          phi[t]*(1-sigma_12)*p1[t], 
          phi[t]*sigma_12*p2[t]),
        
        c(1-phi[t], 
          phi[t]*sigma_21*(1-p1[t]), 
          phi[t]*(1-sigma_21)*(1-p2[t]), 
          phi[t]*sigma_21*p1[t], 
          phi[t]*(1-sigma_21)*p2[t]),
        
        c(1-phi[t], 
          phi[t]*(1-sigma_12)*(1-p1[t]), 
          phi[t]*sigma_12*(1-p2[t]), 
          phi[t]*(1-sigma_12)*p1[t], 
          phi[t]*sigma_12*p2[t]),
        
        c(1-phi[t], 
          phi[t]*sigma_21*(1-p1[t]), 
          phi[t]*(1-sigma_21)*(1-p2[t]), 
          phi[t]*sigma_21*p1[t], 
          phi[t]*(1-sigma_21)*p2[t])
        
        )
}




prop1 <- 0.75
for(i in 1:N){
  ind <- 1*(runif(1) < prop1)
  state <- rep(0,5)
  state[4+ind] <- 1
  
}
