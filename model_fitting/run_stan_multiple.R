options(mc.cores = parallel::detectCores())

setwd("~/github/domestication_selection_model/model_fitting")
library(dplyr)
library(rstan)
library(shinystan)




f <- function(i){
  dat <- read.csv("data/set6_n100.csv")
  dat <- dat[1:i,]
  stan_data <- list(
    N = nrow(dat),
    
    y = dat$y,
    N_t = dat$N,
    m_t = dat$m,
    #m_t_lag = c(0,dat$m[1:(length(dat$m)-1)]),
    n = cbind(dat$n1,dat$n2,dat$n3),
    
    mu_a = 0.25,
    sigma_a = 0.5,
    alpha_b = 2,
    beta_b = 2,
    alpha_s = 2,
    beta_s = 10,
    lambda_z = 1,
    lambda_sigma = 1,
    a_rho = 1,
    b_rho = 1,
    sigma_z = 0.05,
    alpha_sigma_nu = 1,
    beta_sigma_nu = 1
  )
  
  fit = stan(file = "fit_discrete_tufto_2001.stan",
             data = stan_data, 
             chains = 1,
             iter = 10000,
             control=list(max_treedepth = 13, adapt_delta=0.995),
             include = TRUE,
             pars = c("s","z1","r","K"))
  #launch_shinystan(fit)
  dat <- as.data.frame(fit)
  
  write.csv(dat,paste("samples/fit_data3_n100_",i,".csv"))
}


parallel::mclapply(seq(12,35,3), FUN = f, mc.cores = 8)


