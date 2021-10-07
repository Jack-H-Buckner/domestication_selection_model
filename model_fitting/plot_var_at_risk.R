library(dplyr)
library(ggplot2)
setwd("~/github/domestication_selection_model")
dat <- read.csv("model_fitting/samples/test_dat.csv")

ggplot(dat %>% 
         dplyr::filter(!is.na(min)) %>%
         dplyr::group_by(alpha, true)%>%
         dplyr::mutate(n = n())%>%
         dplyr::filter(n > 5), aes(x = as.factor(alpha), 
                                               y = min, color = as.factor(true)))+
  geom_boxplot()



ggplot(dat %>% 
         dplyr::filter(!is.na(min)) %>%
         dplyr::group_by(alpha, true)%>%
         dplyr::summarize(p = sum(min < 0.2)/n()), 
       aes(x = as.factor(alpha), y = p, fill = as.factor(true)))+
  geom_bar(stat = "identity", position = "dodge")



ggplot(dat, 
       aes(x = true, y = exit))+
  geom_bin2d()+
  facet_wrap(~alpha)







dat <- read.csv("model_fitting/samples/test_dat_ab_fixed_priors.csv")


ggplot(dat, aes(x = alpha, y = beta, color = max_z))+
  geom_point()+
  viridis::scale_color_viridis()+
  theme_classic()

ggplot(dat, aes(x = alpha, y = beta, color = iter))+
  geom_point()+
  viridis::scale_color_viridis()+
  theme_classic()

ggplot(dat %>%
         dplyr::mutate(alpha = round(alpha,1), beta = round(beta*2, 1)/2)%>%
         dplyr::group_by(alpha,beta)%>%
         dplyr::summarize(max_z = mean(max_z),
                          mean_z = mean(mean_z),
                          min_N = mean(min_N)), 
       aes(x = alpha, y = beta, fill = max_z))+
  geom_tile()+
  viridis::scale_fill_viridis()+
  theme_classic()











library(ggplot2)
dat <- read.csv("model_fitting/samples/test_dat_ab.csv")






ggplot(dat %>% 
         dplyr::select(-X)%>%
         reshape2::melt(id.vars = c("exit","iter","true","alpha","beta"))%>%
         dplyr::mutate(alpha = 1-round(alpha/2,1)*2)%>%
         dplyr::filter(variable %in% c("min_N", "max_z")), 
       aes(x = as.factor(alpha), y = value))+
  geom_boxplot()+
  facet_wrap(~variable)+
  theme_classic()

ggsave(file = "~/github/domestication_selection_model/figs/alpha_N_z.png",
       width = 4,
       height = 3)







dat %>% 
  dplyr::select(-X)%>%
  reshape2::melt(id.vars = c("exit","iter","true","alpha","beta"))%>%
  dplyr::mutate(alpha = round(alpha/2,1)*2)%>%
  dplyr::group_by(variable,alpha)%>%
  dplyr::summarize(mean = round(mean(value),2),
                   sd = round(sd(value),2))%>%
  dplyr::ungroup()%>%
  dplyr::select(-variable)%>%
  kableExtra::kbl()%>%
  kableExtra::kable_classic_2(html_font = "\"Times New Roman\"",
                              full_width = F )%>%
  kableExtra::kable_styling(font_size = 12)%>%
  kableExtra::pack_rows("Min N", 1, 4) %>%
  kableExtra::pack_rows("Max z", 5, 8) %>%
  kableExtra::pack_rows("Mean z", 9, 12)



ggplot(dat %>% 
         dplyr::select(-X)%>%
         reshape2::melt(id.vars = c("exit","iter","true","alpha","beta"))%>%
         dplyr::mutate(alpha = 1-round(alpha/2,1)*2,
                       beta = round(beta*2,1)/2)%>%
         dplyr::filter(variable == "min_N")%>%
         dplyr::group_by(alpha)%>%
         dplyr::summarize(p = sum(value < 0.5)/n()), 
       aes(x = alpha, y = p))+
  geom_point()+
  geom_line()+
  theme_classic()+
  ylab("P(min{ N_t } < 0.5xK)")

ggsave(file = "~/github/domestication_selection_model/figs/alpha_pN.png",
       width = 3,
       height = 3)



ggplot(dat %>% 
         dplyr::select(-X)%>%
         reshape2::melt(id.vars = c("exit","iter","true","alpha","beta"))%>%
         dplyr::mutate(alpha = 1-round(alpha/2,1)*2)%>%
         dplyr::filter(variable == "max_z")%>%
         dplyr::group_by(alpha)%>%
         dplyr::summarize(p = sum(value > 2.0)/n()), 
       aes(x = alpha, y = p))+
  geom_point()+
  geom_line()+
  theme_classic()+
  ylab("P(max{ z_t } > 2.5)")

ggsave(file = "~/github/domestication_selection_model/figs/alpha_pZ.png",
       width = 3,
       height = 3)



library(ggplot2)
dat1 <- read.csv("~/github/domestication_selection_model/model_fitting/samples/test_dat_ab.csv")
dat2 <- read.csv("~/github/domestication_selection_model/model_fitting/samples/test_dat_ab_20_samples.csv")
#dat3 <- read.csv("~/github/domestication_selection_model/model_fitting/samples/test_dat_ab_100_samples_250_pres.csv")
dat1$design <- 1
dat2$design <- 2
#dat3$design <- 3
dat <- rbind(dat1,dat2)#,dat3)
ggplot(dat %>% 
         dplyr::select(-X)%>%
         reshape2::melt(id.vars = c("design","exit","iter","true","alpha","beta"))%>%
         dplyr::mutate(alpha = 1-round(alpha/2,1)*2,
                       beta = round(beta*2,1)/2)%>%
         dplyr::filter(variable == "min_N")%>%
         dplyr::group_by(alpha,design)%>%
         dplyr::summarize(p = sum(value < 0.5)/n()), 
       aes(x = alpha, y = p, 
           linetype = as.factor(design)))+
  geom_point()+
  geom_line()+
  scale_linetype(name = "Design")+
  theme_classic()+
  ylab("P(min{ N_t } < 0.5xK)")


ggsave(file = "~/github/domestication_selection_model/figs/alpha_pN_samples.png",
       width = 5,
       height = 3)





ggplot(dat %>% 
         dplyr::select(-X)%>%
         reshape2::melt(id.vars = c("design","exit","iter","true","alpha","beta"))%>%
         dplyr::mutate(alpha = 1-round(alpha/2,1)*2,
                       beta = round(beta*2,1)/2)%>%
         dplyr::filter(variable == "max_z")%>%
         dplyr::group_by(alpha,design)%>%
         dplyr::summarize(p = sum(value > 2.0)/n()), 
       aes(x = alpha, y = p, 
           linetype = as.factor(design)))+
  geom_point()+
  geom_line()+
  scale_linetype(name = "Design")+
  theme_classic()+
  ylab("P(min{ z_t } > 2.0)")


ggsave(file = "~/github/domestication_selection_model/figs/alpha_pZ_samples.png",
       width = 5,
       height = 3)


