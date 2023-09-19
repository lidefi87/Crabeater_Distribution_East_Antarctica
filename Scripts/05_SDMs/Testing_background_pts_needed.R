
# Loading libraries -------------------------------------------------------
library(SDMtune)
library(tidyverse)
# library(stars)
# library(cmocean)
library(rsample)
library(recipes)


# Loading data ------------------------------------------------------------
base_env_folder <- "Environmental_Data/ACCESS-OM2-01"

#Background 5x
bg5_file <- file.path(base_env_folder, "Obs_BG_5x_Indian_weaning_LowVIF_match-obs.csv")
mod_bg_5 <- read_csv(bg5_file) %>% 
  drop_na() %>% 
  select(c(xt_ocean, yt_ocean, presence, bottom_slope_deg:dist_ice_edge_km))
#Getting column names of variables to kept in other ACCESS-OM2-01 datasets
env_var <- names(mod_bg_5)

#Background 10x
# bg10_file <- file.path(base_env_folder, "unique_crabeater_obs_all_env.csv")
# mod_bg_10 <- read_csv(bg10_file) %>% 
#   #Selecting observations for the Indian sector during the weaning period
#   filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
#   bind_rows(read_csv(file.path(base_env_folder, "unique_background_10x_obs_all_env.csv"))) %>% 
#   select(all_of(env_var)) %>% 
#   drop_na()
# 
# #Background 20x
# bg20_file <- file.path(base_env_folder, "unique_crabeater_obs_all_env.csv")
# mod_bg_20 <- read_csv(bg20_file) %>% 
#   #Selecting observations for the Indian sector during the weaning period
#   filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
#   bind_rows(read_csv(file.path(base_env_folder, "unique_background_20x_obs_all_env.csv"))) %>% 
#   select(all_of(env_var)) %>% 
#   drop_na()

# Data transformation -----------------------------------------------------
#Defining function to prepare training data
prep_train <- function(da){
  set.seed(42)
  #Dividing dataset (training and testing) - Keep training only
  split <- initial_split(da, prop = 0.75, strata = "presence") 
  train <- training(split) %>% 
    arrange(desc(presence))
  test <- testing(split) %>% 
    arrange(desc(presence))
  #Creating recipe
  recipe_bg <- recipe(presence ~ ., train) %>%
    #Coordinates (xt_ocean and yt_ocean) are not pre-processed
    update_role(xt_ocean, new_role = "coords") %>% 
    update_role(yt_ocean, new_role = "coords") %>% 
    #Scaled and center all predictors - exclude coordinates
    step_center(all_predictors(), -has_role("coords")) %>% 
    step_scale(all_predictors(), -has_role("coords"))
  #Applying recipe to training data
  baked_train <- prep(recipe_bg, training = train) %>% 
    bake(new_data = train)
  #Applying recipe to testing data
  baked_test <- prep(recipe_bg, new_data = test) %>% 
    bake(new_data = test)
  out <- list(baked_train = baked_train,
              baked_test = baked_test)
  return(out)
}

#Baked training data
data_bg5 <- prep_train(mod_bg_5) 
# training_bg10 <- prep_train(mod_bg_10) 
# training_bg20 <- prep_train(mod_bg_20) 


# Using SDM tune formatting -----------------------------------------------
sdm_format <- function(baked_list){
  #Create an empty list to store SDM formatted data
  out <- list(train = NULL, test = NULL)
  for(i in 1:length(baked_list)){
    coords <- baked_list[[i]] %>% 
      arrange(desc(presence)) %>% 
      select(xt_ocean, yt_ocean) %>% 
      rename("X" = "xt_ocean", "Y" = "yt_ocean")
    data <- baked_list[[i]] %>% 
      select(!c(presence, xt_ocean, yt_ocean))
    sdmt_ready <- SWD(species = "Crabeater seals",
                      coords = coords,
                      data = data,
                      pa = baked_list[[i]]$presence)
    if(str_detect(names(baked_list)[i], "train")){
      out$train <- sdmt_ready
    }else if(str_detect(names(baked_list)[i], "test")){
      out$test <- sdmt_ready
      }else{print("testing/training not defined")}
  }
  return(out)
}

sdmt_bg5 <- sdm_format(data_bg5)
# sdmt_bg10 <- sdm_format(training_bg10$train_split, training_bg10$baked_train)
# sdmt_bg20 <- sdm_format(training_bg20$train_split, training_bg20$baked_train)


# Modelling ---------------------------------------------------------------
default_model_bg5 <- train(method = "Maxent", data = sdmt_bg5$train, iter = 5000)
cat("Training auc background 5x: ", auc(default_model_bg5))
plotROC(default_model_bg5)
plotROC(default_model_bg5, test = sdmt_bg5$test)

# default_model_bg10 <- train(method = "Maxent", data = sdmt_bg10, iter = 5000)
# cat("Training auc background 10x: ", auc(default_model_bg10))
# plotROC(default_model_bg10)
# 
# default_model_bg20 <- train(method = "Maxent", data = sdmt_bg20, iter = 5000)
# cat("Training auc background 20x: ", auc(default_model_bg20))
# plotROC(default_model_bg20)
