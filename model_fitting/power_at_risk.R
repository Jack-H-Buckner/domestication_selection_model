library(rstan)
library(dplyr)
library(ggplot2)
library(randomForest)
setwd("~/github/domestication_selection_model")
pred <- get(load("two_populations_conservaiton_hatchery/outputs/predict_levels.RData"))
rf_v1 <- get(load("two_populations_conservaiton_hatchery/outputs/rf_v1.RData"))
rf_v2 <- get(load("two_populations_conservaiton_hatchery/outputs/rf_v2.RData"))
rf_v3 <- get(load("two_populations_conservaiton_hatchery/outputs/rf_v3.RData"))
source("model_fitting/sim_demographic_data.R")
### 
# sample from priors
# simulated data up to start time
#
# fit model and evaluate probability of causing decline 
# if probability above threshold then simulate for another 
# T years w/o releases.
#
# If proability below threshold stop fitting and simulate 
# for T years with releases save min value
#
# Save true paramter values 
#


evalutate_probs <- function(dat, m){
  
  stan_data <- list(
    N = nrow(dat),
    
    y = dat$y,
    N_t = dat$N,
    m_t = dat$m,
    #m_t_lag = c(0,dat$m[1:(length(dat$m)-1)]),
    n = cbind(dat$n1,dat$n2,dat$n3),
    
    alpha_a = 2,
    beta_a = 10,
    alpha_b = 2,
    beta_b = 2,
    alpha_s = 2,
    beta_s = 10,
    alpha_z = 5,
    beta_z = 2,
    lambda_sigma = 1,
    a_rho = 1,
    b_rho = 1,
    sigma_z = 0.05,
    alpha_sigma_nu = 1,
    beta_sigma_nu = 1
  )
  
  samples <- stan(file = "~/github/domestication_selection_model/model_fitting/fit_discrete_tufto_2001.stan",
              data = stan_data, 
              chains = 1,
              iter = 1000,
              control=list(max_treedepth = 11, adapt_delta=0.95),
              include = TRUE,
              pars = c("s","z1","r","K"))

  x <- as.data.frame(samples) %>% select(r,K,s,z1)%>%mutate(mu = m/K)
  x$m <- m
  preds <- pred(x)
  
  p <- sum(preds < 2)/length(preds)
  
  return(p)
  
}

evaluate_min <- function(x, x_lag, m_t, a, b, s, z1, 
                         sigma_Nt, sigma_zt){
  N_t <- c(x[1])
  z_t <- c(x[2])
  # x <- x[1,]
  # x_lag <- x_lag[1,]
  for(i in 1:length(m_t)){
    
    if(i > 1){
      x_new <- x_t(x, x_lag[2], m_t[i], m_t[i-1], a, b, s, z1, rnorm(1,0,sigma_Nt), rnorm(1,0,sigma_zt), 0)
    }else{
      x_new <- x_t(x, x_lag[2], m_t[i], m_t[1], a, b, s, z1, rnorm(1,0,sigma_Nt), rnorm(1,0,sigma_zt), 0)
    }
    x_lag <- x
    x <- x_new$x

    N_t <- append(N_t, x[1])
    z_t <- append(z_t, x[2])
  }

  return(c(min(N_t), max(z_t),mean(z_t)))
}


#' alpha - confidence level to stop
#' beta - confidence to stop monitoring 
#' priors - funciton that samples from prior distribution
#' t - number of years to simulate after sampling
#' N_samples - number for parentage analysis 
#' presc_N - precision of abundnace estimates
samples <- function(alpha, beta, m, n_prior, n_after, interval, 
                    priors, N_samples, presc_N, max_fits){
  pars <- priors()
  m_t <- c(rep(0,n_prior),rep(m, interval*max_fits + n_after))
  n_t <- rep(N_samples,n_prior+interval*max_fits + n_after)
  sigma_obs <- 1/presc_N
  x0 <- c(1.0,0)

  
  dat_ls <- data(x0, x0, m_t, n_t, pars$r, pars$r/pars$K, 
                       pars$s, pars$z1, pars$sigma, 0.01, sigma_obs)


  dat <- dat_ls$dat
  
  true <- pred(data.frame(r=pars$r, K=pars$K, s=pars$s, z1=pars$z1, m=m, mu=m/pars$K))

  x0 <- dat_ls$x
  x0_lag <- dat_ls$x_lag
  p <- evalutate_probs(dat[1:n_prior,], m)


  # reset m_t and n_t
  m_t <- rep(m,interval)
  n_t <- rep(N_samples,interval)
  print(p)
  i <- 1
  while(((p > alpha) & ((1-p) > beta) & i < max_fits)|i==1){

    write.csv(dat,"~/github/domestication_selection_model/model_fitting/data/test_dat.csv")
    
    x0 <- dat[(n_prior+i*interval),c("N","z")]
    x0 <- c(x0$N,x0$z)
    x0_lag <- dat[(n_prior+i*interval - 1),c("N","z")]
    x0_lag <- c(x0_lag$N,x0_lag$z)
    p <- evalutate_probs(dat[1:(n_prior+i*interval),], m)
    print(p)
    i <- i + 1
  }

  if(p <= alpha){
    m_t <- rep(0,n_after)
    n_t <- rep(0,n_after)
    vals <- evaluate_min(x0, x0_lag, m_t, pars$r, pars$r/pars$K,
                         pars$s, pars$z1, pars$sigma, 0.01)
    
    return(c(min_N = vals[1]/pars$K, max_z = vals[2], mean_z = vals[3], exit = 0, iter = i, 
             true = true, alpha = alpha, beta = beta))
  }
  else{
    m_t <- rep(m,n_after)
    n_t <- rep(0,n_after)
    vals <- evaluate_min(x0, x0_lag, m_t, pars$r, pars$r/pars$K,
                          pars$s, pars$z1, pars$sigma, 0.01)

    
    return(c(min_N = vals[1]/pars$K, max_z = vals[2], mean_z = vals[3],exit = 1, iter = i,
             true = true, alpha = alpha, beta = beta))
  }
}

alpha_a = 2
beta_a = 10
alpha_b = 2
beta_b = 2
alpha_s = 2
beta_s = 10
alpha_z = 5
beta_z = 2
lambda_sigma = 0.025


priors <- function(){
  r <- rgamma(1,alpha_a,beta_a)
  if(r < 0.025){
    r <- 0.025
  }else if(r > 2.0){
    r <- 2.0
  }
  K <- 1#rgamma(1,alpha_b,beta_b)
  s <- rgamma(1,alpha_s,beta_s)
  if(s > 0.75){
    s <- 0.75
  }
  z <- rgamma(1,alpha_z,beta_z)
  if(z > 5){
    z <- 5
  }
  sigma <- 0.1#rexp(1,1/lambda_sigma)
  return(list(r=r,K=K,s=s,z1=z,sigma=sigma))
}


priors_fixed <- function(){
  r <- 0.15
  K <- 1
  s <- 0.15
  z <- 3.5
  sigma <- 0.1
  return(list(r=r,K=K,s=s,z1=z,sigma=sigma))
}


alpha <- 0.10
beta <- 0.05
m <- 0.1  
n_prior <- 10
n_after <- 40
interval <- 5 
N_samples <- 100 
presc_N <- 2
max_fits <- 15

model <- stan_model("~/github/domestication_selection_model/model_fitting/fit_discrete_tufto_2001.stan")
s <- function(pars){
  return(samples(pars[1], pars[2], m, n_prior, n_after, interval, 
                 priors, N_samples, presc_N, max_fits))
}


N <- 800
alpha_ls <- runif(N)*0.6+0.2
beta_ls <- runif(N) * (1-alpha_ls)*0.25

sums <- c()
for(i in 1:N){
  sums <- append(sums, alpha_ls[i]+ beta_ls[i])
}
plot(alpha_ls,beta_ls)
abline(a = 1, b = -1)
hist(sums)

pars_ls <- list(c(alpha_ls[1],beta_ls[1]))
for(i in 2:N){
  pars_ls <- append(pars_ls,list(c(alpha_ls[i],beta_ls[i])))
}

ls <- parallel::mclapply(pars_ls,
                         s,mc.cores = 9)
dat <- t(as.data.frame(ls))

write.csv(dat, "~/github/domestication_selection_model/model_fitting/samples/test_dat_ab_100_samples_250_pres.csv")



