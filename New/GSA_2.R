setwd("~/github/domestication_selection_model/New")
library(dplyr)
library(modeest)
library(PNWColors)
library(ggplot2)

### load name and clean data
data <- read.csv("data/GSA_data.csv")

names(data) <- c("Wmin", "Nmin","recovery50", "recovery90", "pim", "RRS", "T", "Policy", "k","r","s","density_dependence")

data <- data %>% filter(Wmin > 0)


### simple plots
pal=pnw_palette("Cascades",2, type = "discrete")
data1 <- data
data1$`T` <- plyr::mapvalues(data$`T`, sort(unique(data$`T`)), c("0.5 generations", "1 generation", "2 generations", "4 generations"))
data1$Policy <- plyr::mapvalues(data$Policy, c(1,0), c("Fixed","Feedback"))
p <- ggplot(data1, aes(y = Wmin, x = as.factor(RRS), fill = as.factor(Policy)))+
  geom_boxplot()+
  facet_wrap(~`T`)+
  theme_classic()+
  scale_fill_manual(values =pal, name = "Policy")+
  scale_color_manual(values =pal)+
  xlab("Reletive Fitness")+
  ylab("Minimum Fitness")+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )


ggsave(file = "figures/Effect of feedback policy.png",
       p,
       height = 7.5,
       width = 7.5)


pal=pnw_palette("Cascades",5, type = "discrete")
p <- ggplot(data1 %>% 
  group_by(`T`, pim, Policy, RRS) %>%
  summarize(Wmin = mean(Wmin)) %>%
  reshape2::dcast( RRS + pim + `T` ~ Policy) %>%
  mutate(difference = Feedback - Fixed),
  aes(x = pim, y = difference, color = as.factor(RRS)))+
  geom_point(size = 3)+
  geom_line(size = 1.5)+
  facet_wrap(~`T`)+
  theme_classic()+
  scale_color_manual(values =pal, name = "RRS")+
  xlab("Immigraiton rate")+
  ylab("differnce in minimum fitness")+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )

ggsave(file = "figures/Differnce between feedback and fixed policy.png",
       p,
       height = 7.5,
       width = 7.5)
  


pal=pnw_palette("Anemone",3, type = "discrete")
data1 <- data
data1$`T` <- plyr::mapvalues(data$`T`, sort(unique(data$`T`)), c("0.5", "1.0", "2.0", "4.0"))
data1$Policy <- plyr::mapvalues(data$Policy, c(1,0), c("Fixed","Feedback"))
p <- ggplot(data1, aes(y = recovery50, x = as.factor(`T`), fill = as.factor(s)))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values =pal, name = "Selection \n strength")+
  scale_color_manual(values =pal)+
  xlab("Duration (generations)")+
  ylab("Tiem to 50% recovery")+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )
p

ggsave(file = "figures/50 recovery.png",
       p,
       height = 7.5,
       width = 7.5)


pal=pnw_palette("Anemone",3, type = "discrete")
data1 <- data
data1$`T` <- plyr::mapvalues(data$`T`, sort(unique(data$`T`)), c("0.5", "1.0", "2.0", "4.0"))
data1$Policy <- plyr::mapvalues(data$Policy, c(1,0), c("Fixed","Feedback"))
p <- ggplot(data1, aes(y = recovery90, x = as.factor(`T`), fill = as.factor(s)))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values =pal, name = "Selection \n strength")+
  scale_color_manual(values =pal)+
  xlab("Duration (generations)")+
  ylab("Tiem to 90% recovery")+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15)
  )
p

ggsave(file = "figures/90 recovery.png",
       p,
       height = 7.5,
       width = 7.5)



KLdivergence <- function(marginal_samples, conditional_samples){
  marginal_range <- range(marginal_samples)
  marginal_pdf <- density(marginal_samples, n = 512, from = marginal_range[1],to = marginal_range[2]) 
  conditional_pdf <- density(conditional_samples, n = 512, from = marginal_range[1],to = marginal_range[2]) 
  step =  marginal_pdf$x[2] - marginal_pdf$x[1]
  conditional_y <- conditional_pdf$y[conditional_pdf$y > 0.0]
  marginal_y <- marginal_pdf$y[conditional_pdf$y > 0.0]
  return(step*sum(conditional_y*log(conditional_y/marginal_y,2)))
} 


EKLdivergence <- function(data, quantity, condition){
  marginal_samples <- data[,quantity]
  acc <- 0
  n <- 0
  for(c in unique(data[,condition])){
    conditional_samples <- data[data[,condition] == c,quantity]
    n <- n + 1
    v <- KLdivergence(marginal_samples, conditional_samples)

    acc <- acc + v
  }
  return(acc/n) 
}








variable_importance <- function(data, outcomes, variables){
  nout <- length(outcomes)
  nvar <- length(variables)
  nrows <-  nout * nvar 
  VarImp <- data.frame(outcome = rep(0, nrows), variable= rep(0,nrows), EKL = rep(0,nrows))
  
  for(i in 1:nrows){
    outcome_ind <- floor((i-1)/nvar)+1
    var_ind <- i -floor((i-1)/nvar) * nvar
    
    VarImp[i,] <- c(outcomes[outcome_ind],variables[var_ind], EKLdivergence(data,outcomes[outcome_ind], variables[var_ind]))
    
    
  }
  return(VarImp)
}




outcomes <- c("Wmin", "Nmin","recovery50", "recovery90")
variables <- c("pim", "RRS", "T", "Policy", "k","r","s","density_dependence")

VarImp <- variable_importance(data, outcomes,variables)

#"Cascades"
#"Starfish"
pal=pnw_palette("Starfish",4, type = "discrete")
VarImp$outcome <- ordered(VarImp$outcome, c("Wmin", "Nmin", "recovery50", "recovery90") )
VarImp$outcome <- plyr::mapvalues(VarImp$outcome, c("Wmin", "Nmin", "recovery50", "recovery90") , c("Min. fitness", "Min. abundance", "50% recovery", "90% recovery"))

VarImp$variable <- ordered(VarImp$variable, rev(c("RRS", "pim", "T", "k", "density_dependence","r","s","Policy") ))
VarImp$variable <- plyr::mapvalues(VarImp$variable, c("RRS", "pim", "T", "k", "density_dependence","r","s","Policy"), c("RRS", "Prop. im", "Duration", "k", "Density","r","s","Policy"))

p <- ggplot(VarImp, aes(x = as.numeric(EKL), y = variable, fill = outcome))+
  geom_bar(stat = "identity")+
  facet_wrap(~outcome)+
  scale_fill_manual(values =pal)+
  theme_classic()+
  xlab("Variable Importance (bits)")+
  ylab("Variable")+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.position = "none"
  )




ggsave(file = "figures/Variable Importance.png",
       p,
       height = 7.5,
       width = 7.5)



outcomes <- c("Wmin", "Nmin","recovery50", "recovery90")
variables <- c("pim", "RRS", "T", "Policy", "k","r","s","density_dependence")

VarImp_fixed <- variable_importance(data %>% filter(Policy == 1, `T` > 200), outcomes,variables)
VarImp_var <- variable_importance(data %>% filter(Policy == 0, `T` > 200), outcomes,variables)
VarImp_fixed$Policy <- "Fixed"
VarImp_var$Policy <- "Adaptive"
VarImp <- rbind(VarImp_fixed,VarImp_var)
#"Cascades"
#"Starfish"
pal=pnw_palette("Lake",2, type = "discrete")
VarImp$outcome <- ordered(VarImp$outcome, c("Wmin", "Nmin", "recovery50", "recovery90") )
VarImp$outcome <- plyr::mapvalues(VarImp$outcome, c("Wmin", "Nmin", "recovery50", "recovery90") , c("Min. fitness", "Min. abundance", "50% recovery", "90% recovery"))

VarImp$variable <- ordered(VarImp$variable, rev(c("RRS", "pim", "T", "k", "density_dependence","r","s","Policy") ))
VarImp$variable <- plyr::mapvalues(VarImp$variable, c("RRS", "pim", "T", "k", "density_dependence","r","s","Policy"), c("RRS", "Prop. im", "Duration", "k", "Density","r","s","Policy"))

p <- ggplot(VarImp, aes(x = as.numeric(EKL), y = variable, fill = Policy))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~outcome)+
  scale_fill_manual(values =pal, "Policy")+
  theme_classic()+
  xlab("Variable Importance (bits)")+
  ylab("Variable")+
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    strip.text = element_text(size = 20)
  )




ggsave(file = "figures/Variable Importance adaptive policy.png",
       p,
       height = 7.5,
       width = 7.5)
