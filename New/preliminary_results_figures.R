library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
setwd("~/github/domestication_selection_model/New")

#######################################################
##### make plots for effect of RRS at equilibrium #####
#######################################################

# check axis values match those used in figures.jl
# fitness for DSI
dat <- read.csv("data/Fitness_RSSByPim_s01.csv")

Rim <- seq(0.0,0.55,0.05)
RRS <- seq(0.1,0.9,0.05)^2
dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$Rim <- Rim
names(dat) <- c(RRS, "Rim")

dat_long <- dat%>%melt(id.var = "Rim")
dat_long$RRS <- as.numeric(as.character(dat_long$variable))

ggplot(dat_long, aes(x = RRS, y = value, color = Rim, group = as.factor(Rim)))+
  geom_line()+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  xlab("Reletive Fitness (Hatchery:Natrual)")+
  ylab("Prop. Survive Selection")

ggsave("figures/RRS_fitness_eq_DSI.png",
       height = 4.0,
       width = 5.0)
# Abundance for DSI
dat <- read.csv("data/Abundance_RSSByPim_s01.csv")

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$Rim <- Rim
names(dat) <- c(RRS, "Rim")

dat_long <- dat%>%melt(id.var = "Rim")
dat_long$RRS <- as.numeric(as.character(dat_long$variable))

ggplot(dat_long, aes(x = RRS, y = value, color = Rim, group = as.factor(Rim)))+
  geom_line()+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  xlab("Reletive Fitness (Hatchery:Natrual)")+
  ylab("Spawning Stock Size")

ggsave("figures/RRS_abundance_eq_DSI.png",
       height = 4.0,
       width = 5.0)

# fitness for DIS
dat <- read.csv("data/Fitness_RSSByPim_s01_DIS.csv")

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$Rim <- Rim
names(dat) <- c(RRS, "Rim")

dat_long <- dat%>%melt(id.var = "Rim")
dat_long$RRS <- as.numeric(as.character(dat_long$variable))

ggplot(dat_long, aes(x = RRS, y = value, color = Rim, group = as.factor(Rim)))+
  geom_line()+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  xlab("Reletive Fitness (Hatchery:Natrual)")+
  ylab("Prop. Survive Selection")

ggsave("figures/RRS_fitness_eq_DIS.png",
       height = 4.0,
       width = 5.0)



# Abundance for DIS
dat <- read.csv("data/Abundance_RSSByPim_s01_DIS.csv")

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$Rim <- Rim
names(dat) <- c(RRS, "Rim")

dat_long <- dat%>%melt(id.var = "Rim")
dat_long$RRS <- as.numeric(as.character(dat_long$variable))

ggplot(dat_long, aes(x = RRS, y = value, color = Rim, group = as.factor(Rim)))+
  geom_line()+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  xlab("Reletive Fitness (Hatchery:Natrual)")+
  ylab("Spawning Stock Size")

ggsave("figures/RRS_abundance_eq_DIS.png",
       height = 4.0,
       width = 5.0)



#########################################
##### make plots for ratchet effect #####
#########################################


dat <- read.csv("data/Fitness_pimByRSS_s01.csv")
Rim = seq(0.1,0.7,0.025)^2
RRS = seq(0.01,0.91,0.3)

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$Rim <- Rim
names(dat) <- c(RRS, "Rim")


dat_long <- dat%>%melt(id.var = "Rim")
dat_long$RRS <- as.numeric(as.character(dat_long$variable))

ggplot(dat_long, aes(x = Rim, y = value, color = RRS, group = as.factor(RRS)))+
  geom_line()+
  geom_point()+
  scale_color_viridis(option = "plasma")+
  theme_classic()+
  xlab("Immigration: Recruitment")+
  ylab("Prop. Survive Selection")



ggsave("figures/Rim_fitness_eq.png",
       height = 4.0,
       width = 5.0)

## make plot for abundance

dat <- read.csv("data/Abundance_pimByRSS_s01.csv")
Rim = seq(0.1,0.7,0.025)^2
RRS = seq(0.01,0.91,0.3)

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$Rim <- Rim
names(dat) <- c(RRS, "Rim")


dat_long <- dat%>%melt(id.var = "Rim")
dat_long$RRS <- as.numeric(as.character(dat_long$variable))
ggplot(dat_long, aes(x = Rim, y = value, color = RRS, group = as.factor(RRS)))+
  geom_line()+
  geom_point()+
  scale_color_viridis(option = "plasma")+
  theme_classic()+
  xlab("Immigration: Recruitment")+
  ylab("Prop. Survive Selection")



ggsave("figures/Rim_abundance_eq.png",
       height = 4.0,
       width = 5.0)






######################################
##### Effect of program duration #####
######################################



dat <- read.csv("data/Fitness_TimByPim_RRS02.csv")
Rim = seq(0.01,0.31,0.05)
T_ = 2^seq(-2,6,0.5)

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$T_ <- T_
names(dat) <- c(Rim, "T_")


dat_long <- dat%>%melt(id.var = "T_")
dat_long$Rim<- as.numeric(as.character(dat_long$variable))

ggplot(dat_long, aes(x = T_, y = value, color = Rim, group = as.factor(Rim)))+
  geom_point()+
  geom_line()+
  scale_color_viridis()+
  theme_classic()+
  xlab("Duration of releases (generations)")+
  ylab("Prop. Survive Selection")


ggsave("figures/T_fitness_eq.png",
       height = 4.0,
       width = 5.0)




ggplot(dat_long %>% filter(T_ < 5), aes(x = T_, y = value, color = Rim, group = as.factor(Rim)))+
  geom_point()+
  geom_line()+
  scale_color_viridis()+
  theme_classic()+
  xlab("Duration of releases (generations)")+
  ylab("Prop. Survive Selection")


ggsave("figures/T_fitness_eq_zoom.png",
       height = 4.0,
       width = 5.0)



###################################################
##### Effect of program duration on threshold #####
###################################################




dat <- read.csv("data/Fitness_TimByPim_s01_RSS01.csv")


Rim <- seq(0.01,0.81,0.02) 

Tgen =  DemographicParameters.Smyth_2016_T1
T_ <- 2^seq(-2,5,1.0)

dat <- rbind(as.numeric(substr(names(dat),2,nchar(names(dat))-7)), dat)
dat$T_ <- T_
names(dat) <- c(Rim, "T_")


dat_long <- dat%>%melt(id.var = "T_")
dat_long$Rim<- as.numeric(as.character(dat_long$variable))




ggplot(dat_long, aes(x = Rim, y = value, color = T_, group = as.factor(T_)))+
  geom_point()+
  geom_line()+
  scale_color_viridis(option = "cividis", trans = "log")+
  theme_classic()+
  xlab("Duration of releases (generations)")+
  ylab("Prop. Survive Selection")


ggsave("figures/Rim_T_fitness_eq_zoom.png",
       height = 4.0,
       width = 5.0)
