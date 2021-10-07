library(dplyr)
setwd("~/github/domestication_selection_model/model_fitting")


x_t <- function(x, z_t1, m_t, m_t1, a, b, s, z1, epsilon_Nt, epsilon_zt, n_samps){
  N <- x[1]
  z <- x[2]
  
  mu <- m_t1/N


  N_t <- m_t + N*exp(a - 1/2*s/(1+s)*z^2 - b*N+epsilon_Nt)
  z_t <- (N_t-m_t)*z/(N_t*(1+s)) + m_t/N_t*z1+epsilon_zt
 
  p1 <- (1-mu)^2 *exp(- 1/2 *s/(1+s)*z_t1^2)  
  p2 <- 2*mu*(1-mu)*exp( - 1/2 *s/(1+s)*((z_t1+z1)/2)^2)
  p3 <- mu^2*exp(- 1/2 *s/(1+s)*z1^2 )
  

  p1_prime <- p1/(p1+p2+p3)
  p2_prime <- p2/(p1+p2+p3)
  p3_prime <- p3/(p1+p2+p3)
  
  
  samp <- rmultinom(n=1,size = n_samps,prob=c(p1_prime,p2_prime,p3_prime))
  
  return(list(x = c(N_t,z_t), samp = samp))
}






data <- function(x0, x0_lag, m_t, n_t,a, b, s, z1, sigma_Nt, sigma_zt, sigma_obs){
  N_t <- c()
  z_t <- c()
  y_t <- c()
  n1_ls <- c()
  n2_ls <- c()
  n3_ls <- c()
  x <- x0
  x_lag <- x0_lag
  for(i in 1:length(m_t)){
 
    N_t1 <- x[1]

    if(i > 1){
      x_new <- x_t(x, x_lag[2], m_t[i], m_t[i-1], a, b, s, z1, rnorm(1,0,sigma_Nt), rnorm(1,0,sigma_zt), n_t[i])
    }else{
      x_new <- x_t(x, x_lag[2], m_t[i], m_t[1], a, b, s, z1, rnorm(1,0,sigma_Nt), rnorm(1,0,sigma_zt), n_t[i])
    }

    x_lag <- x
    x <- x_new$x
    samp <- x_new$samp
    y <- log((x[1] - m_t[i])/N_t1)
    y_t <- append(y_t, y)
    y_mod2 <- append(y_t, y)
    N_t <- append(N_t, x[1])
    z_t <- append(z_t,x[2])
    n1_ls <- append(n1_ls, samp[1])
    n2_ls <- append(n2_ls, samp[2])
    n3_ls <- append(n3_ls, samp[3])
  }

  dat <- data.frame(y = y_t, N = exp(log(N_t) + rnorm(length(N_t),0,sigma_obs)), 
                    m = m_t, z = z_t,
                    n1 = n1_ls, n2 = n2_ls, 
                    n3 = n3_ls)
  
  return(list(dat = dat, x = x, x_lag = x_lag))
}

generate_data <- function(path, x0, m_t, n_t,a, b, s, z1, sigma_Nt, sigma_zt, sigma_obs){
  d_ls <- data(x0, x0, m_t, n_t,a, b, s, z1, sigma_Nt, sigma_zt, sigma_obs)
  dat <- d_ls$dat
  write.csv(dat, path)
  return(list(x = d_ls$x, x_lag = d_ls$x_lag))
}





add_rows <- function(dat, x0, x0_lag, 
                     m_t, n_t,a, b, s, z1, 
                     sigma_Nt, sigma_zt, sigma_obs){
              
  d_ls <- data(x0, x0_lag, m_t, n_t,a, b, s, z1, sigma_Nt, sigma_zt, sigma_obs)
  dat <- rbind(dat,d_ls$dat)
  
  return(list(dat = dat, x = d_ls$x, x_lag = d_ls$x_lag))
}

