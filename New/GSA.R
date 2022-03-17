setwd("~/github/domestication_selection_model/New")
library(dplyr)
library(ggplot2)
GSA_data <- read.csv("data/Removals_GSA_data.csv")


names(GSA_data) <- c("age_params","order","T_", "RRS", "s", "k", "p_im","F_","W")

GSA_data <-  GSA_data %>% filter(T_ != 0)
nrow(GSA_data)
min_W <- 0.85

GSA_data_filtered <- GSA_data %>% 
  dplyr::group_by(age_params,order,T_, RRS,s,k, p_im) %>%
  dplyr::mutate(max_W = max(W)) %>%
  dplyr::ungroup()%>%
  dplyr::filter(W == max_W | W > min_W )%>% 
  dplyr::group_by(age_params,order,T_, RRS,s,k, p_im) %>%
  dplyr::summarize(outcome = min(F_))
  

# computes the expected increase in informaiton 
# about the outcome when you condition on the 
# optimal policy 
expected_KLD <- function(GSA_data_filtered, variable){
  
  # compute  the marginal distribuiton 
  dsn <- GSA_data_filtered %>%
    dplyr::group_by()%>%
    dplyr::mutate(N = n())%>%
    dplyr::group_by_("outcome")%>%
    dplyr::summarize(p = n()/mean(N))
  
  # compute the conditional distributions 
  cnd_dsn <- GSA_data_filtered %>%
    dplyr::group_by_(variable)%>%
    dplyr::mutate(N = n())%>%
    dplyr::ungroup()%>%
    dplyr::group_by_(variable,"outcome")%>%
    dplyr::summarize(p = n()/mean(N))

  # variable and level names 
  outcome_levels <- dsn %>% group_by_("outcome") %>% summarize(n = n())
  variable_levels <- cnd_dsn %>% group_by_(variable) %>% summarize(n = n())
  
  levels <- dsn$outcome
  dsn <- dsn$p
  
  # set accumulator
  EKLD <- 0
 
  i <- 0
 
  # loop over each level of the variable 
  # compute the KL divergence between the marginal 
  # and conditional for each level of the variable 
  for(j in 1:nrow(variable_levels[variable])){

    i <-  i + 1
    v <- as.numeric(variable_levels[j,variable])
  
    cdsn <- cnd_dsn$p[cnd_dsn[variable] == v] 
    clevels <- cnd_dsn$outcome[cnd_dsn[variable] == v] 
 
    inds <- c()
    # loop over outcome levels to find values
    # with non-zero probabilty in the conditional distribution  
    for(i in 1:length(levels)){ 
      if(levels[i] %in% clevels){
          inds <- append(inds,i)
        }
      }
    
    # calcualte KL divergence for current variable level 
    # add it to accumumulator 
    dsni <- dsn[inds]
    EKLD <- EKLD + sum(cdsn * log(cdsn/dsni)) 
  }
  # return average KLD
  return(EKLD/i)
}


# computes the expected increase in informaiton 
# about the outcome when you condition on the 
# optimal policy 
expected_KLD_2d <- function(GSA_data_filtered, variable1, variable2){
  if(variable1 == variable2){
    return(expected_KLD(GSA_data_filtered, variable1))
  }else{
    # compute  the marginal distribuiton 
    dsn <- GSA_data_filtered %>%
      dplyr::group_by()%>%
      dplyr::mutate(N = n())%>%
      dplyr::group_by_("outcome")%>%
      dplyr::summarize(p = n()/mean(N))
    
    # compute the conditional distributions 
    cnd_dsn <- GSA_data_filtered %>%
      dplyr::group_by_(variable1, variable2)%>%
      dplyr::mutate(N = n())%>%
      dplyr::ungroup()%>%
      dplyr::group_by_(variable1, variable2,"outcome")%>%
      dplyr::summarize(p = n()/mean(N))
  
    # variable and level names 
    outcome_levels <- dsn %>% group_by_("outcome") %>% summarize(n = n())
    variable_levels <- cnd_dsn %>% group_by_(variable1, variable2) %>% summarize(n = n())
    
    levels <- dsn$outcome
    dsn <- dsn$p
    
    # set accumulator
    EKLD <- 0
    
    i <- 0
    
    # loop over each level of the variable 
    # compute the KL divergence between the marginal 
    # and conditional for each level of the variable 
    for(j in 1:nrow(variable_levels[c(variable1, variable2)])){
      
      i <-  i + 1
      v <- as.numeric(variable_levels[j,c(variable1, variable2)])
      v1 <-  v[1]
      v2 <- v[2]
      cdsn <- cnd_dsn$p[cnd_dsn[variable1] == v1 & cnd_dsn[variable2] == v2 ] 
      clevels <- cnd_dsn$outcome[cnd_dsn[variable1] == v1 & cnd_dsn[variable2] == v2 ] 
   
      inds <- c()
      # loop over outcome levels to find values
      # with non-zero probabilty in the conditional distribution  
      for(i in 1:length(levels)){ 
        if(levels[i] %in% clevels){
          inds <- append(inds,i)
        }
      }
      
      # calcualte KL divergence for current variable level 
      # add it to accumumulator 
      dsni <- dsn[inds]
      EKLD <- EKLD + sum(cdsn * log(cdsn/dsni)) 
    }
    # return average KLD
    return(EKLD/i)
  }
}

outcome <- "F_opt"
variable <- "k"
vars <- c("age_params","order","T_", "RRS", "s", "k", "p_im")
var_imp <- c()
for(var in vars){
  var_imp <- append(var_imp, expected_KLD(GSA_data_filtered, var))
}

var_imp_dat <- data.frame(var = vars, Imp = var_imp)

ggplot(var_imp_dat,
       aes(x = Imp, y = var))+
  geom_bar(stat = "identity")+
  theme_classic()


# 2d variable importance

var_imp_2d <- c()
vars1 <- c()
vars2 <- c()
i <- 0
for(var1 in vars){
  i <- i + 1
  j <- 0
  for(var2 in vars){
    j <- j + 1
    if(i < j){
      vars1 <- append(vars1, var1)
      vars2 <- append(vars2, var2)
      var_imp_2d <- append(var_imp_2d, expected_KLD_2d(GSA_data_filtered, var1, var2))
    }
  }
}

var_imp_dat_2d <- data.frame(var1 = vars1, var2 = vars2 , Imp = var_imp_2d)
var_imp_dat_2d$var1 <- factor(var_imp_dat_2d$var1, vars)
var_imp_dat_2d$var2 <- factor(var_imp_dat_2d$var2, rev(vars))
ggplot(var_imp_dat_2d,
       aes(x = var2, y = var1, fill = Imp))+
  geom_tile()+
  theme_classic()+
  viridis::scale_fill_viridis()


