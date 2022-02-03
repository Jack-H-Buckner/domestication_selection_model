setwd("~/github/domestication_selection_model/White_sturgeon")

df <- read.csv("model_output/DSI_param_sensetivity.csv")

library(ggplot2)
library(randomForest)

names(df) <- c("R_H", "Duration", 
               "s", "RRS", "k_ind", 
               "min_fittness","min_abundance",
               "T_min_abundnace","integrated_fitness")

rf_min_fittness <- randomForest::randomForest(y = df[,c("min_fittness")], 
                                              x = df[,c("R_H", "Duration", 
                                                       "s", "RRS", "k_ind")])

randomForest::varImpPlot(rf_min_fittness)

rf_min_abundance <- randomForest::randomForest(y = df[,c("min_abundance")], 
                                              x = df[,c("R_H", "Duration", 
                                                        "s", "RRS", "k_ind")])

randomForest::varImpPlot(rf_min_abundance)



rf_integrated_fitness <- randomForest::randomForest(y = df[,c("integrated_fitness")], 
                                               x = df[,c("R_H", "Duration", 
                                                         "s", "RRS", "k_ind")])

randomForest::varImpPlot(rf_integrated_fitness)


rf_T_min_abundnace <- randomForest::randomForest(y = df[,c("T_min_abundnace")], 
                                                    x = df[,c("R_H", "Duration", 
                                                              "s", "RRS", "k_ind")])

randomForest::varImpPlot(rf_T_min_abundnace)


ggplot(df, 
       aes(y = min_fittness, 
           fill = as.factor(RRS), 
           x = as.factor(Duration)))+
  geom_boxplot()+
  theme_classic()


ggplot(df, 
       aes(y = min_abundance, 
           fill = as.factor(RRS), 
           x = as.factor(Duration)))+
  geom_boxplot()+
  theme_classic()

randomForest::partialPlot(rf_integrated_fitness, x.var = "RRS",
                          pred.data = df[,c("R_H", "Duration", 
                                            "s", "RRS", "k_ind")])


