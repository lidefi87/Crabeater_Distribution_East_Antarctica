Calculating threshold for individual SDM algorithms
================
Denisse Fierro Arcos
2023-12-14

- <a href="#calculating-threshold-for-individual-sdm-algorithms"
  id="toc-calculating-threshold-for-individual-sdm-algorithms">Calculating
  threshold for individual SDM algorithms</a>
  - <a href="#loading-libaries" id="toc-loading-libaries">Loading
    libaries</a>
  - <a
    href="#loading-models-and-predicting-presence-access-om2-01-matching-observations"
    id="toc-loading-models-and-predicting-presence-access-om2-01-matching-observations">Loading
    models and predicting presence (ACCESS-OM2-01 matching observations)</a>
    - <a href="#calculating-thresholds-for-each-sdm"
      id="toc-calculating-thresholds-for-each-sdm">Calculating thresholds for
      each SDM</a>
  - <a href="#loading-models-and-predicting-presence-access-om2-01-full-set"
    id="toc-loading-models-and-predicting-presence-access-om2-01-full-set">Loading
    models and predicting presence (ACCESS-OM2-01 full set)</a>
    - <a href="#calculating-thresholds-for-each-sdm-1"
      id="toc-calculating-thresholds-for-each-sdm-1">Calculating thresholds
      for each SDM</a>
  - <a href="#loading-models-and-predicting-presence-observations"
    id="toc-loading-models-and-predicting-presence-observations">Loading
    models and predicting presence (Observations)</a>
    - <a href="#calculating-thresholds-for-each-sdm-2"
      id="toc-calculating-thresholds-for-each-sdm-2">Calculating thresholds
      for each SDM</a>

# Calculating threshold for individual SDM algorithms

We calculated max SSS threshold using an adapted form of the
`thresholds` function from the `SDMTune` package to find the probability
above which the crabeaters are considered to occur. We will use these
threshold to create binary maps of presence that we will use as a
measure of uncertainty. Any probabilities above or equal to the
threshold will be classified as presence (value of 1), while any values
below the threshold were considered absence (value of 0).

## Loading libaries

``` r
library(tidyverse)
library(mgcv)
library(SDMtune)
library(randomForest)
library(gbm)
library(prg)
library(pROC)
source("Scripts/05_SDMs/useful_functions.R")
```

## Loading models and predicting presence (ACCESS-OM2-01 matching observations)

``` r
mod_match_obs <- read_csv("Environmental_Data/mod-match-obs_env_pres_bg_20x_Indian_weaning.csv") %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month))

#Preparing training data and testing data for GAM
mod_match_obs <- prep_data(mod_match_obs, "month", split = T)

#Applying SWD format for all other algorithms
mod_data_sdm <- mod_match_obs$baked_test %>% 
  select(!year) %>% 
  sdm_format() 
```

### Calculating thresholds for each SDM

``` r
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
```

## Loading models and predicting presence (ACCESS-OM2-01 full set)

``` r
mod_full <- read_csv("Environmental_Data/model_env_pres_bg_20x_Indian_weaning.csv") %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month))

#Preparing training data and testing data for GAM
mod_full <- prep_data(mod_full, "month", split = T)

#Applying SWD format for all other algorithms
mod_full_data_sdm <- mod_full$baked_test %>% 
  select(!year) %>% 
  sdm_format() 
```

### Calculating thresholds for each SDM

``` r
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
```

## Loading models and predicting presence (Observations)

``` r
obs <- read_csv("Environmental_Data/obs_env_pres_bg_20x_Indian_weaning.csv") %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month))

#Preparing training data and testing data for GAM
obs <- prep_data(obs, "month", split = T)

#Applying SWD format for all other algorithms
obs_data_sdm <- obs$baked_test %>% 
  select(!year) %>% 
  sdm_format() 
```

### Calculating thresholds for each SDM

``` r
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
```
