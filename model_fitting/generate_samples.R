setwd("~/github/domestication_selection_model")
source("model_fitting/sim_demographic_data.R")

path <- "data/set0_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,20),c(0.1,0.1,0.1,0.1))#,c(0.1,0.1,0.1,0.1)
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.000
sigma_obs <- 0.00
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)

path <- "data/set0_n100_prime.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,20))#,c(0.1,0.1,0.1,0.1)
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.00
sigma_obs <- 0.0
dat_ls <- data(x0, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
               sigma_Nt, sigma_zt, sigma_obs)
dat <- dat_ls$dat
dat <- add_rows(dat, dat_ls$x, dat_ls$x_lag, 
                c(0.1,0.1,0.1,0.1), c(100,100,100,100),a, b, s, z1, 
                sigma_Nt, sigma_zt, sigma_obs)
write.csv(dat$dat,"data/set0_n100_prime.csv")



path <- "model_fitting/data/set1.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,15),rep(0.0,15),rep(0.2,15),rep(0.0,10))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.005
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(0, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)





path <- "model_fitting/data/set2.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,15),rep(0.0,10))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(0, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)


path <- "model_fitting/data/set3.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,7))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 2
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(0, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)


path <- "model_fitting/data/set4.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,4),seq(0,0.2,0.025),rep(0.2,10), rep(0.0,4))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 2
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(0, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)










path <- "model_fitting/data/set1_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,15),rep(0.0,15),rep(0.2,15),rep(0.0,10))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.005
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)





path <- "model_fitting/data/set2_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,15),rep(0.0,10))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)


path <- "model_fitting/data/set3_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,7))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 2
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)



path <- "model_fitting/data/set4_n50.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.1,0.025),rep(0.1,20))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 2
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(50, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)


path <- "model_fitting/data/set4_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.1,0.025),rep(0.1,20))
a <- 0.25
b <- 0.25
s <- 0.1
z1 <- 2
sigma_Nt <- 0.05
sigma_zt <- 0.05
sigma_obs <- 0.025
generate_data(path, x0, m_t, rep(50, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)


path <- "model_fitting/data/set5_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.1,0.025),rep(0.1,20))
a <- 0.15
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.1
sigma_zt <- 0.05
sigma_obs <- 0.05
generate_data(path, x0, m_t, rep(50, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)



path <- "model_fitting/data/set6_n100.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.1,0.025),rep(0.1,20))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 2
sigma_Nt <- 0.05
sigma_zt <- 0.05
sigma_obs <- 0.025
generate_data(path, x0, m_t, rep(50, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)









path <- "model_fitting/data/set1_n500.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,15),rep(0.0,15),rep(0.2,15),rep(0.0,10))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.15
sigma_zt <- 0.005
sigma_obs <- 0.1
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)





path <- "model_fitting/data/set2_n500.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,15),rep(0.0,10))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 4
sigma_Nt <- 0.15
sigma_zt <- 0.05
sigma_obs <- 0.1
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)


path <- "model_fitting/data/set3_n500.csv"
x0 <- c(1.0,0.0)
m_t <- c(rep(0,10),seq(0,0.2,0.05),rep(0.2,7))
a <- 0.25
b <- 0.25
s <- 0.2
z1 <- 2
sigma_Nt <- 0.15
sigma_zt <- 0.05
sigma_obs <- 0.1
generate_data(path, x0, m_t, rep(100, length(m_t)),a, b, s, z1, 
              sigma_Nt, sigma_zt, sigma_obs)





