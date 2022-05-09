library(ggplot2)
library(dplyr)
library(reshape2)
library(PNWColors)
setwd("~/github/domestication_selection_model/New")

dat_DSI <- read.csv("data/temporary_long_run_comparison_DSI.csv")
dat_DIS <- read.csv("data/temporary_long_run_comparison_DIS.csv")


names(dat_DSI) <- c("Ratio. im", "RRS", "Duration", "W eq.", "W min.", "SSB eq.", "SSB min.")
names(dat_DIS) <- c("Ratio. im", "RRS", "Duration", "W eq.", "W min.", "SSB eq.", "SSB min.")


### make plor for immigration after selection ####

dat_DSI$Duration <- plyr::mapvalues(as.factor(dat_DSI$Duration),c("0.5","1","2","4"),
                                    c("Duration = 0.5 gen.","Duration = 1.0 gen.","Duration = 2.0 gen.","Duration = 4.0 gen.")) 
ggplot(dat_DSI %>% filter(Duration == "Duration = 1.0 gen."),
       aes(x = `Ratio. im`, y = RRS, fill = `W eq.` ))+
  geom_tile()+
  theme_classic()+
  xlab("Hatchery immigration : Natrual recruitment")+
  ylab("Reletive Fitness ")+
  scale_fill_gradientn(colours = pnw_palette("Lake",n = 100,type="continuous"),
                       name = "Eq. Fitness") +
  theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 16),
        legend.position = "bottom")

ggsave(file = "figures/equilibrium_heatmap_s005.png",
       height = 7,
       width = 7)


dat <- dat_DSI %>% mutate(`W diff.` = `W min.` - `W eq.`)

ggplot(dat,
       aes(x = `Ratio. im`, y = RRS, fill = `W diff.` ))+
  geom_tile()+
  facet_wrap(~Duration, ncol = 2)+
  theme_classic()+
  scale_fill_gradientn(colours = pnw_palette("Lake",n = nrow(dat_DSI),type="continuous")) #



### use linear interpolation to smooth surface ###


# define grid of interpolation points
ratioim_ls <- c()
RRS_ls <- c()
RRS <- seq(0.1,0.9,0.005)
ratioim <- seq(0.025,1.45,0.005)
for(rrs in RRS){
  for(rim in ratioim){
    ratioim_ls <- append(ratioim_ls,rim)
    RRS_ls <- append(RRS_ls,rrs)
  }
}


### define gam interpolation ###
dat <- dat_DSI %>% 
  mutate(Wdiff = `W min.` - `W eq.`,
         Ratioim = `Ratio. im`) %>% 
  filter(Duration == "Duration = 2.0 gen.")

z_dat <- dat %>% 
  dcast(RRS~Ratioim, value.var = "Wdiff") %>%
  as.matrix()
  
pred <- pracma::interp2(unique(dat$Ratioim), 
                        unique(dat$RRS),
                        z_dat[,2:ncol(z_dat)],
                        ratioim_ls, RRS_ls,
                        method = "linear")


dat_grid <- data.frame(Ratioim = ratioim_ls,RRS =RRS_ls)
dat_interp <- dat_grid
dat_interp$pred <- pred 

ggplot(dat_interp,
       aes(x = Ratioim, y = RRS, fill = pred))+
  geom_tile()+
  theme_classic()+
  scale_fill_gradientn(colours = pnw_palette("Lake",n = nrow(dat_interp),type="continuous"))

### make plot for immigration before selection ####

dat <- dat_DIS %>% mutate(`W diff.` = `W min.` - `W eq.`)

ggplot(dat,
       aes(x = `Ratio. im`, y = RRS, fill = `W diff.` ))+
  geom_tile()+
  facet_wrap(~Duration, ncol = 2)+
  theme_classic()+
  scale_fill_gradientn(colours = pnw_palette("Lake",n = nrow(dat_DSI),type="continuous"))





