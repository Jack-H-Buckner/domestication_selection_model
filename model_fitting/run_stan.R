#options(mc.cores = parallel::detectCores())
setwd("~/github/domestication_selection_model/model_fitting")
library(dplyr)
library(rstan)
library(shinystan)





dat <- read.csv("~/github/domestication_selection_model/model_fitting/data/set3_n500.csv")


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
  lambda_sigma = 2,
  a_rho = 3,
  b_rho = 1,
  sigma_z = 0.001,
  alpha_sigma_nu = 1,
  beta_sigma_nu = 1
)

fit1 = stan(file = "~/github/domestication_selection_model/model_fitting/fit_discrete_tufto_2001.stan",
           data = stan_data, 
           chains = 1,
           iter = 8000,
           control=list(max_treedepth = 13, adapt_delta=0.995),
           include = TRUE,
           pars = c("s","z1","r","K", "z"))

fit2 = stan(file = "~/github/domestication_selection_model/model_fitting/fit_discrete_tufto_2001_N.stan",
           data = stan_data, 
           chains = 1,
           iter = 8000,
           control=list(max_treedepth = 13, adapt_delta=0.995),
           include = TRUE,
           pars = c("s","z1","r","K"))


fit3 = stan(file = "~/github/domestication_selection_model/model_fitting/fit_discrete_tufto_2001_samples.stan",
           data = stan_data, 
           chains = 1,
           iter = 8000,
           control=list(max_treedepth = 13, adapt_delta=0.995),
           include = TRUE,
           pars = c("s","z1","r","K"))





dat <- as.data.frame(fit1)
write.csv(dat,paste("samples/fit_data3_n500_all.csv"))

dat <- as.data.frame(fit2)
write.csv(dat,paste("samples/fit_data3_n500_N.csv"))

dat <- as.data.frame(fit3)
write.csv(dat,paste("samples/fit_data3_n500_samples.csv"))
# 
# d <- as.data.frame(fit1)
# 
# ggplot(d, aes(x = s, y = z1))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()
# 
# 
# ggplot(d, aes(x = log(z1), y = log(`z[21]`)))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()
# 
# 
# 
# ggplot(d, aes(x = log(s), y = log(z1)))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()+
#   geom_smooth(method = "lm", color = "red")
# 
# 
# orthoganal_projection <- function(x,y,u){
#   return(sum(u*c(x,y))/sum(u*u))
# }
# 
# 
# 
# 
# ### extra stuff
# 
# p1 <- ggplot(as.data.frame(fit1),
#              aes(x = s, y = z1))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()+
#   ylim(0,8)+
#   xlim(0,1.4)+
#   theme(legend.position =  "None")
# 
# p2 <- ggplot(as.data.frame(fit2),
#              aes(x = s, y = z1))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()+
#   ylim(0,8)+
#   xlim(0,1.4)+
#   theme(legend.position =  "None")
# 
# p3 <- ggplot(as.data.frame(fit3),
#              aes(x = s, y = z1))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()+
#   ylim(0,8)+
#   xlim(0,1.4)+
#   theme(legend.position =  "None")
# 
# ggsave(file = "~/github/domestication_selection_model/figs/fit1_corr_dat3.png",
#        p1,
#        height = 3,
#        width = 3)
# 
# ggsave(file = "~/github/domestication_selection_model/figs/fit2_corr_dat3.png",
#        p2,
#        height = 3,
#        width = 3)
# 
# ggsave(file = "~/github/domestication_selection_model/figs/fit3_corr_dat3.png",
#        p3,
#        height = 3,
#        width = 3)
# 
# 
# pc <- prcomp(log(d[c("s","z1")]))
# 
# d$s_scale <- log(d$s) - mean(log(d$s))
# d$z1_scale <- log(d$z1) - mean(log(d$z1))
# u <- pc$rotation[,1]
# 
# proj <- c()
# for(i in 1:nrow(d)){
#   proj <- append(proj, orthoganal_projection(log(d$s[i]),log(d$z1[i]),u))
# }
# d$proj <- proj
# ggplot(d, aes(x = proj))+
#   geom_density()+
#   theme_classic()
# 
# ggplot(d, aes(x = proj))+
#   geom_density(linetype = 2)+
#   theme_classic()+
#   geom_function(fun = dnorm,
#               args = list(mean = mean(proj), sd = sd(proj)),
#               linetype = 1, color = "red")
# 
# ggplot(d, aes(x = s_scale, y = z1_scale))+
#   geom_hex()+
#   viridis::scale_fill_viridis()+
#   theme_classic()+
#   geom_abline(aes(intercept = 0, slope = u[2]/u[1]),
#               color = "red")
# 
# 
# 
# pc_comps <- princomp(log(d[c("s","z1")]))
# 
# hist(pc_comps$loadings[,1])
# 
