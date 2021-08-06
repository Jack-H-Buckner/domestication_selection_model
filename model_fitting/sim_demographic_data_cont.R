library(deSolve)
library(dplyr)

dx <- function(t,x,par){
  # state
  N <- x[1]
  z <- x[2]
  
  # parameters
  r <- par[1]
  K <- par[2]
  s <- par[3]
  z1 <- par[4]
  
  dN <- N*(r*(1-N/K) - s/2*(1+z^2)) + m(t)
  dz <- -s*z - m(t)*(z-z1)/N
  
  return(list(c(dN,dz)))
}


m <- function(t){
  0.2*exp(-0.5*(t-50)^2)+0.2*exp(-0.2*(t-55)^2)
}

parameters <- c(0.25,1,0.025,2)
state <- c(1,0)
times <- seq(0, 90, by = 0.01)
out <- ode(y = state, times = times, func = dx, parms = parameters)

plot(out)


