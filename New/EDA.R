setwd("~/github/domestication_selection_model/New")
library(dplyr)
library(ggplot2)
GSA_data <- read.csv("data/Removals_GSA_data.csv")
names(GSA_data) <- c("age_params","order","T_", "RRS", "s", "k", "p_im","F_","W")
GSA_data <-  GSA_data %>% filter(T_ != 0)



ggplot(GSA_data,
       aes(x =F_,y = W, color = as.factor(T_))) +
  geom_point()+geom_smooth()+
  facet_wrap(~order)+
  theme_classic()
nrow(GSA_data)