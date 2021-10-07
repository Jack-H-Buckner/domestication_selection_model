library(randomForest)
library(ggplot2)
library(dplyr)
setwd("~/github/domestication_selection_model/two_populations_conservaiton_hatchery")

dat <- read.csv("outputs/ricker_GSA.csv")

dat <- dat %>% dplyr::mutate(v1 = (eq1-m) < K,
                             v2 = eq1 < K,
                             v3 = eq2 + 0.05 < eq1,
                             mu = m/K) 

dat1 <- dat
dat1$v1 <- as.factor(as.character(data$v1))
inds <- sample(1:nrow(dat1), round(nrow(dat1)/10), replace = F)
y <- as.factor(dat1$v1)
x <- dat1 %>% select(s,z1,m,r,K,mu)

rf_v1 <- randomForest(y = y[-inds],
                   x = x[-inds,],
                   y_test = y[inds],
                   x_test = x[inds,])

save(rf_v1,file = "outputs/rf_v1.RData")

preds <- predict(rf_v1, x[inds,])
sum(as.numeric(preds) == as.numeric(y[inds]))/length(inds)


partialPlot(rf_v1, x[inds,], "s")
partialPlot(rf_v1, x[inds,], "z1")
partialPlot(rf_v1, x[inds,], "m")
partialPlot(rf_v1, x[inds,], "r")
partialPlot(rf_v1, x[inds,], "K")
partialPlot(rf_v1, x[inds,], "mu")

dat2 <- dat %>% filter(v1)
inds <- sample(1:nrow(dat2), round(nrow(dat2)/10), replace = F)
y <- as.factor(as.character(dat2$v2))
x <- dat2 %>% select(s,z1,m,r,K,mu)

rf_v2 <- randomForest(y = y[-inds],
                   x = x[-inds,],
                   y_test = y[inds],
                   x_test = x[inds,])
save(rf_v2,file = "outputs/rf_v2.RData")


preds <- predict(rf_v2, x[inds,])
sum(as.numeric(preds) == as.numeric(y[inds]))/length(inds)


partialPlot(rf, x[inds,], "s")
partialPlot(rf, x[inds,], "z1")
partialPlot(rf, x[inds,], "m")
partialPlot(rf, x[inds,], "r")
partialPlot(rf, x[inds,], "K")
partialPlot(rf, x[inds,], "mu")



dat3 <- dat %>% filter(v2)
inds <- sample(1:nrow(dat3), round(nrow(dat3)/10), replace = F)
y <- as.factor(as.character(dat3$v3))
x <- dat3 %>% select(s,z1,m,r,K,mu)

rf_v3 <- randomForest(y = y[-inds],
                       x = x[-inds,],
                       y_test = y[inds],
                       x_test = x[inds,])

save(rf_v3,file = "outputs/rf_v3.RData")

preds <- predict(rf_v3, x[inds,])
sum(as.numeric(preds) == as.numeric(y[inds]))/length(inds)


partialPlot(rf_v3, x[-inds,], "s")
partialPlot(rf_v3, x[-inds,], "z1")
partialPlot(rf_v3, x[-inds,], "m")
partialPlot(rf_v3, x[-inds,], "r")
partialPlot(rf_v3, x[-inds,], "K")
partialPlot(rf_v3, x[-inds,], "mu")



pred <- function(x){
  preds1 <- predict(rf_v1,x)
  preds2 <- predict(rf_v2, x)
  preds3 <- predict(rf_v3, x)
  level1 <- (preds1 == "TRUE") & !(preds2 == "TRUE"|preds3 == "TRUE")
  level2 <- preds1 == "TRUE" & preds2 == "TRUE" & !(preds3 == "TRUE")
  level3 <- preds1 == "TRUE" & preds2 == "TRUE" & preds3 == "TRUE"
  values <- 1*level1 + 2*level2 + 3*level3
  return(values)
}





save(pred,file = "outputs/predict_levels.RData")








