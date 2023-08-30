
# Loading libraries -------------------------------------------------------
library(SDMtune)
library(tidyverse)
library(stars)
library(cmocean)
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
bg10_file <- file.path(base_env_folder, "unique_crabeater_obs_all_env.csv")
mod_bg_10 <- read_csv(bg10_file) %>% 
  #Selecting observations for the Indian sector during the weaning period
  filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
  bind_rows(read_csv(file.path(base_env_folder, "unique_background_10x_obs_all_env.csv"))) %>% 
  select(all_of(env_var)) %>% 
  drop_na()

#Background 20x
bg20_file <- file.path(base_env_folder, "unique_crabeater_obs_all_env.csv")
mod_bg_20 <- read_csv(bg20_file) %>% 
  #Selecting observations for the Indian sector during the weaning period
  filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
  bind_rows(read_csv(file.path(base_env_folder, "unique_background_20x_obs_all_env.csv"))) %>% 
  select(all_of(env_var)) %>% 
  drop_na()

# Data transformation -----------------------------------------------------
#Defining function
prep_train <- function(da){
  set.seed(42)
  #Creating recipe
  recipe_bg <- recipe(presence ~ ., da %>% select(-c(xt_ocean, yt_ocean))) %>% 
    step_center(all_predictors()) %>% 
    step_scale(all_predictors())
  #Dividing dataset (training and testing) - Keep training only
  split_bg <- initial_split(da, prop = 0.75, strata = "presence") 
  train_bg <- training(split_bg) %>% 
    arrange(desc(presence))
  #Prepping recipe
  prep_bg <- prep(recipe_bg, training = train_bg)
  #Applying recipe to training data
  baked_bg <- bake(prep_bg, new_data = train_bg)
  out <- list(train_split = train_bg,
              baked_train = baked_bg)
  return(out)
}

#Baked training data
training_bg5 <- prep_train(mod_bg_5) 
training_bg10 <- prep_train(mod_bg_10) 
training_bg20 <- prep_train(mod_bg_20) 


# Using SDM tune formatting -----------------------------------------------
sdm_format <- function(da_original, da_train){
  coords <- da_original %>% 
    arrange(desc(presence)) %>% 
    select(xt_ocean, yt_ocean) %>% 
    rename("x" = "xt_ocean", "y" = "yt_ocean")
  sdmt_bg <- SWD(species = "Crabeater seals",
                  coords = coords,
                  data = da_train %>% select(-presence),
                  pa = da_train$presence)
  return(sdmt_bg)
}

sdmt_bg5 <- sdm_format(training_bg5$train_split, training_bg5$baked_train)
sdmt_bg10 <- sdm_format(training_bg10$train_split, training_bg10$baked_train)
sdmt_bg20 <- sdm_format(training_bg20$train_split, training_bg20$baked_train)


# Modelling ---------------------------------------------------------------
default_model_bg5 <- train(method = "Maxent", data = sdmt_bg5, iter = 5000)
cat("Training auc background 5x: ", auc(default_model_bg5))
plotROC(default_model_bg5)

default_model_bg10 <- train(method = "Maxent", data = sdmt_bg10, iter = 5000)
cat("Training auc background 10x: ", auc(default_model_bg10))
plotROC(default_model_bg10)

default_model_bg20 <- train(method = "Maxent", data = sdmt_bg20, iter = 5000)
cat("Training auc background 20x: ", auc(default_model_bg20))
plotROC(default_model_bg20)
