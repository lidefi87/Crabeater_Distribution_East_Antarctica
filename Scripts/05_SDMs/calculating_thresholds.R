## Threshold calculation
library(tidyverse)
library(mgcv)
library(SDMtune)
library(randomForest)
library(gbm)
library(prg)
library(pROC)
source("Scripts/05_SDMs/useful_functions.R")



# Confusion matrix function -----------------------------------------------
#Adapted from SDMTune
confMatrix_adap <- function(model,
                       data,
                       th = NULL,
                       test = NULL,
                       type = NULL) {
  
  if (is.null(test)) {
    data <- data
  } else {
    data <- test
  }
  
  n_p <- sum(data$presence == 1)
  n_a <- sum(data$presence == 0)
  pred <- predict(model, data, type = type)
  p_pred <- pred[1:n_p]
  a_pred <- pred[(n_p + 1):(n_p + n_a)]
  
  if (is.null(th)) {
    th <- sort(unique(pred))
    th <- c(0, th, 1)
  }
  
  tp <- fp <- vector(mode = "numeric", length = length(th))
  
  for (i in seq_along(th)) {
    tp[i] <- sum(p_pred >= th[i])
    fp[i] <- sum(a_pred >= th[i])
  }
  
  fn <- n_p - tp
  tn <- n_a - fp
  
  data.frame(th = th,
             tp = tp,
             fp = fp,
             fn = fn,
             tn = tn)
  
}


# Threshold function ------------------------------------------------------
#Adapted from SDMTune
thresholds_adap <- function(model,
                       data,
                       type = NULL,
                       test = NULL) {
  
  n_pres <- sum(data$presence == 1)
  
  cm_train <- confMatrix_adap(model, data, type = type)
  tpr <- cm_train$tp / (cm_train$tp + cm_train$fn)
  tnr <- cm_train$tn / (cm_train$fp + cm_train$tn)
  fpr <- cm_train$fp / (cm_train$fp + cm_train$tn)

  mtp <- min(predict(model, data, type = type))
  ess <- cm_train$th[which.min(abs(tpr - tnr))]
  mss <- cm_train$th[which.max(tpr + tnr)]

  ths <- c(mtp, ess, mss)
  rownames <- c("Minimum training presence",
                "Equal training sensitivity and specificity",
                "Maximum training sensitivity plus specificity")
  colnames <- c("Threshold",
                paste(stringr::str_to_title(type), "value"),
                "Fractional predicted area",
                "Training omission rate")
  
  if (!is.null(test)) {
    cm_test <- confMatrix_adap(model, data, type = type, test = test)
    tpr_test <- cm_test$tp / (cm_test$tp + cm_test$fn)
    tnr_test <- cm_test$tn / (cm_test$fp + cm_test$tn)
    
    ess <- cm_test$th[which.min(abs(tpr_test - tnr_test))]
    mss <- cm_test$th[which.max(tpr_test + tnr_test)]
    
    ths <- c(ths, ess, mss)
    rownames <- c(rownames,
                  "Equal test sensitivity and specificity",
                  "Maximum test sensitivity plus specificity")
    colnames <- c(colnames, "Test omission rate", "P-values")
    n_test <- nrow(data)
    or_test <- vector(mode = "numeric", length = 5)
    p_values <- vector(mode = "numeric", length = 5)
  }
  
  or_train <- vector(mode = "numeric", length = length(ths))
  fpa <- vector(mode = "numeric", length = length(ths))
  
  for (i in seq_along(ths)) {
    index <- which.min(abs(cm_train$th - ths[i]))
    or_train[i] <- cm_train[index, ]$fn / n_pres
    fpa[i] <- fpr[index]
    
    if (!is.null(test)) {
      index <- which.min(abs(cm_test$th - ths[i]))
      or_test[i] <- cm_test[index, ]$fn / n_test
      p_values[i] <- stats::binom.test((round((1 - or_test[i]), 0) * n_test),
                                       n_test, fpa[i],
                                       alternative = "greater")$p.value
    }
  }
  
  output <- data.frame(th = rownames, val = ths, fpa = fpa, or = or_train,
                       stringsAsFactors = FALSE)
  
  if (!is.null(test)) {
    output$or_test <- or_test
    output$pv <- p_values
  }
  
  colnames(output) <- colnames
  
  output
}



# Loading data (ACCESS-OM2-01 matching observations) ----------------------
mod_match_obs <- read_csv("Environmental_Data/mod-match-obs_env_pres_bg_20x_Indian_weaning.csv") %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month))

#Preparing training data and testing data for GAM
mod_match_obs <- prep_data(mod_match_obs, "month", split = T)

#Applying SWD format for all other algorithms
mod_data_sdm <- mod_match_obs$baked_test %>% 
  select(!year) %>% 
  sdm_format() 


#Calculating thresholds
#GAM
gam_mod <- readRDS("SDM_outputs/GAM/Mod_match_obs/best_GAM_mod_match_obs.rds")
gam_thres <- thresholds_adap(gam_mod, mod_match_obs$baked_train, 
                             type = "response", test = mod_match_obs$baked_test)
gam_thres %>% 
  write_csv("SDM_outputs/GAM/Mod_match_obs/threshold_calculation.csv")

#Maxent
maxent_mod <- readRDS("SDM_outputs/Maxent/Mod_match_obs/reduced_Maxent_model/best_red_maxent_model.rds")
max_thres <- thresholds_adap(maxent_mod, mod_match_obs$baked_train, 
                             type = "cloglog", test = mod_match_obs$baked_test)
max_thres %>% 
  write_csv("SDM_outputs/Maxent/Mod_match_obs/threshold_calculation.csv")

#Random Forest
rf_mod <- readRDS("SDM_outputs/RandomForest/Mod_match_obs/reduced_RF_mod_match_obs.rds")
rf_thres <- thresholds_adap(rf_mod, mod_match_obs$baked_train, 
                            type = "response", test = mod_match_obs$baked_test)
rf_thres %>% 
  write_csv("SDM_outputs/RandomForest/Mod_match_obs/threshold_calculation.csv")

#Boosted Regression Trees
brt_mod <- readRDS("SDM_outputs/BoostedRegressionTrees/Mod_match_obs/best_BRT_mod_match_obs.rds")
brt_thres <- thresholds_adap(brt_mod, mod_match_obs$baked_train, 
                            type = "response", test = mod_match_obs$baked_test)
brt_thres %>% 
  write_csv("SDM_outputs/BoostedRegressionTrees/Mod_match_obs/threshold_calculation.csv")



# Loading models and predicting presence (ACCESS-OM2-01 full set) ---------
#Loading data
mod_full <- read_csv("Environmental_Data/model_env_pres_bg_20x_Indian_weaning.csv") %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month))

#Preparing training data and testing data for GAM
mod_full <- prep_data(mod_full, "month", split = T)

#Applying SWD format for all other algorithms
mod_full_data_sdm <- mod_full$baked_test %>% 
  select(!year) %>% 
  sdm_format() 

#Calculating thresholds
#GAM
gam_mod <- readRDS("SDM_outputs/GAM/Mod_full/best_GAM_mod_full.rds")
gam_thres <- thresholds_adap(gam_mod, mod_full$baked_train, 
                             type = "response", test = mod_full$baked_test)
gam_thres %>% 
  write_csv("SDM_outputs/GAM/Mod_full/threshold_calculation.csv")

#Maxent
maxent_mod <- readRDS("SDM_outputs/Maxent/Mod_full/initial_Maxent_model/model.Rds")
max_thres <- thresholds_adap(maxent_mod, mod_full$baked_train, 
                             type = "cloglog", test = mod_full$baked_test)
max_thres %>% 
  write_csv("SDM_outputs/Maxent/Mod_full/threshold_calculation.csv")

#Random Forest
rf_mod <- readRDS("SDM_outputs/RandomForest/Mod_full/model.Rds")
rf_thres <- thresholds_adap(rf_mod, mod_full$baked_train, 
                            type = "response", test = mod_full$baked_test)
rf_thres %>% 
  write_csv("SDM_outputs/RandomForest/Mod_full/threshold_calculation.csv")

#Boosted Regression Trees
brt_mod <- readRDS("SDM_outputs/BoostedRegressionTrees/Mod_full/model.Rds")
brt_thres <- thresholds_adap(brt_mod, mod_full$baked_train, 
                             type = "response", test = mod_full$baked_test)
brt_thres %>% 
  write_csv("SDM_outputs/BoostedRegressionTrees/Mod_full/threshold_calculation.csv")


# Loading models and predicting presence (Observations) -------------------
obs <- read_csv("Environmental_Data/obs_env_pres_bg_20x_Indian_weaning.csv") %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month))

#Preparing training data and testing data for GAM
obs <- prep_data(obs, "month", split = T)

#Applying SWD format for all other algorithms
obs_data_sdm <- obs$baked_test %>% 
  select(!year) %>% 
  sdm_format() 

#Calculating thresholds
#GAM
gam_mod <- readRDS("SDM_outputs/GAM/Obs/best_GAM_obs.rds")
gam_thres <- thresholds_adap(gam_mod, obs$baked_train, type = "response", 
                             test = obs$baked_test)
gam_thres %>% 
  write_csv("SDM_outputs/GAM/Obs/threshold_calculation.csv")


#Maxent
maxent_mod <- readRDS("SDM_outputs/Maxent/Obs/reduced_Maxent_model/best_red_maxent_model.rds")
max_thres <- thresholds_adap(maxent_mod, obs$baked_train, type = "cloglog", 
                             test = obs$baked_test)
max_thres %>% 
  write_csv("SDM_outputs/Maxent/Obs/threshold_calculation.csv")


#Random Forest
rf_mod <- readRDS("SDM_outputs/RandomForest/Obs/model.Rds")
rf_thres <- thresholds_adap(rf_mod, obs$baked_train, type = "response", 
                            test = obs$baked_test)
rf_thres %>% 
  write_csv("SDM_outputs/RandomForest/Obs/threshold_calculation.csv")


#Boosted Regression Trees
brt_mod <- readRDS("SDM_outputs/BoostedRegressionTrees/Obs/model.Rds")
brt_thres <- thresholds_adap(brt_mod, obs$baked_train, type = "response", 
                             test = obs$baked_test)
brt_thres %>% 
  write_csv("SDM_outputs/BoostedRegressionTrees/Obs/threshold_calculation.csv")



