# Loading libraries
library(tidyverse)
library(rsample)
library(dplyr, lib.loc = '~/R/rocker-rstudio/4.2')
library(recipes)
library(caret)

#Loading data
crab_ind <- read_csv("Environmental_Data/Env_obs/Obs_BG_5x_Indian_weaning_LowVIF.csv") %>% 
  #Removing any rows with NA values
  drop_na()

#Selecting variables (response and predictors) to be included in the model
model_data <- crab_ind %>% 
  select(presence:dist_ice_edge_km)

#Calculating weights to be applied to presence column
wt <- model_data %>% 
  count(presence) %>% 
  pivot_wider(names_from = presence, values_from = n) %>% 
  mutate(weight = `1`/`0`) %>% 
  pull(weight)

#Getting a list of numeric predictive variables that are NOT presence
env_var <-  model_data %>% 
  select(where(is.numeric) & !presence) %>% 
  names()

#Building formulas for GAM
#Full model
gam_fm <- paste("presence ~", paste(paste0("s(", env_var, ")"), collapse = " + "))

#Model scopes - to be tested
scope <- list(
  "bottom_slope_deg" = ~1 + bottom_slope_deg + s(bottom_slope_deg, df = 2) + s(bottom_slope_deg, df = 3) + s(bottom_slope_deg, df = 4) + s(bottom_slope_deg, df = 5),
  "dist_shelf_km" = ~1 + dist_shelf_km + s(dist_shelf_km, df = 2) + s(dist_shelf_km, df = 3) + s(dist_shelf_km, df = 4) + s(dist_shelf_km, df = 5),
  "dist_coast_km" = ~1 + dist_coast_km + s(dist_coast_km, df = 2) + s(dist_coast_km, df = 3) + s(dist_coast_km, df = 4) + s(dist_coast_km, df = 5),
  "depth_m" = ~1 + depth_m + s(depth_m, df = 2) + s(depth_m, df = 3) + s(depth_m, df = 4) + s(depth_m, df = 5),
  "SIC" = ~1 + SIC + s(SIC, df = 2) + s(SIC, df = 3) + s(SIC, df = 4) + s(SIC, df = 5),
  "SST_degC" = ~1 + SST_degC + s(SST_degC, df = 2) + s(SST_degC, df = 3) + s(SST_degC, df = 4) + s(SST_degC, df = 5),
  "lt_pack_ice" = ~1 + lt_pack_ice + s(lt_pack_ice, df = 2) + s(lt_pack_ice, df = 3) + s(lt_pack_ice, df = 4) + s(lt_pack_ice, df = 5),
  "dist_ice_edge_km" = ~1 + dist_ice_edge_km + s(dist_ice_edge_km, df = 2) + s(dist_ice_edge_km, df = 3) + s(dist_ice_edge_km, df = 4) + s(dist_ice_edge_km, df = 5)
)


#Splitting data - Training and testing
set.seed(42)
#Ensuring training and test data keep similar obs/background ratios
split <- initial_split(crab_ind, prop = 0.75, strata = "presence")
train <- training(split)
test <- testing(split)
#Getting weights for train and test sets
train_wt <- if_else(train$presence == 0, wt, 1)

#Preparing train data 
baked_train <- recipe(presence ~., data = train %>% select(presence:dist_ice_edge_km, month)) %>% 
  #Standardising data
  step_center(all_numeric(), -all_outcomes()) %>% 
  step_scale(all_numeric(), -all_outcomes()) %>% 
  #Applying recipe to training dataset
  prep(train) %>% 
  bake(train)


# GAM step-wise -----------------------------------------------------------
library(gam)
gam.object <- gam(formula = presence ~s(SIC), data = baked_train,
                  family = binomial(link = "cloglog"),
                  weights = train_wt, select = T)
step.object <- step.Gam(gam.object, scope = scope, direction = "both")
step.object_for <- step.Gam(gam.object, scope = scope, direction = "forward")

gam.object.b <- gam(presence ~ s(bottom_slope_deg, df = 5) + s(dist_shelf_km, df = 5) + s(dist_coast_km, df = 5) + s(depth_m, df = 5) + s(SIC, df = 5) + s(SST_degC, df = 5) + s(lt_pack_ice, df = 5) + s(dist_ice_edge_km, df = 5),
                    data = baked_train,
                    family = binomial(link = "cloglog"),
                    weights = train_wt)
step.object_back <- step.Gam(gam.object.b, scope = scope, direction = "backward")

mean_model <- read_csv("Environmental_Data/ACCESS-OM2-01/All_values_ACCESS-OM2-01_env_vars.csv")

mean_model <- mean_model %>% 
  select(!c(xt_ocean, yt_ocean)) %>%
  scale() %>% 
  as.data.frame() %>% 
  cbind(mean_model %>% select(c(xt_ocean, yt_ocean))) %>% 
  drop_na()

res_back <- mean_model %>% 
  mutate(pred = as.vector(predict(step.object_back,
                                  mean_model %>% select(bottom_slope_deg:dist_ice_edge_km), 
                                  type = "response"))) 

library(stars)
library(cmocean)
ras_pred <- res_back %>%
  select(xt_ocean:pred) %>% 
  st_as_stars()


ggplot() +
  geom_stars(data = ras_pred) +
  scale_fill_cmocean(name = "ice", direction = -1, 
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(0, 1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Crabeater seal distribution during weaning season (observations)",
       x = "Longitude",
       y = "Latitude",
       fill = "Probability")#+



res_back %>% 
  ggplot(aes(xt_ocean, yt_ocean, color = pred))+
  geom_point()

res_back %>% 
  write_csv("obs-pred_mean_env_best_gam_model.csv")


# pres <- res_back %>% 
#   ggplot(aes(xt_ocean, yt_ocean, color = presence))+
#   geom_point()
# 
# pred <- res_back %>% 
#   ggplot(aes(xt_ocean, yt_ocean, color = pred))+
#   geom_point()

#Applying GAM
library(mgcv)
library(splines)
mod_gam <- gam(formula = presence ~ ns(dist_shelf_km, df = 4) + dist_coast_km + ns(depth_m, df = 5) + ns(SIC, df = 5) + ns(SST_degC, df = 5) + ns(dist_ice_edge_km, df = 5), 
               data = baked_train, family = binomial(link = "cloglog"), 
               weights = train_wt, method = "REML")
summary(mod_gam)


gam_fm <- paste("presence ~", paste0(env_var, collapse = " + "))
cv_mod <- train(as.formula(gam_fm), data = baked_train,
                         family = binomial(link = "cloglog"),
                         weights = train_wt, method = "gam",
                         trControl = trainControl(method = "cv", number = 10))
summary(cv_mod)


#Checking predictions
res <- test %>%  
  mutate(pred = predict(mod_gam, test, type = "response"))

pres <- res %>% 
  ggplot(aes(xt_ocean, yt_ocean, color = presence))+
  geom_point()

pred <- res %>% 
  ggplot(aes(xt_ocean, yt_ocean, color = pred))+
  geom_point()

mean_model <- read_csv("Cleaned_Data/Env_obs/All_values_ACCESS-OM2-01_env_vars.csv")


pred_mean_model <- mean_model %>%  
  mutate(pred = as.vector(predict(mod_gam, 
                        mean_model %>% select(bottom_slope_deg:dist_ice_edge_km),
                        type = "response")))

pred_mean_model %>% 
  ggplot(aes(xt_ocean, yt_ocean, color = pred))+
  geom_point()

library(patchwork)
pres/pred



# crab_ind %>% 
#   select(bottom_slope_deg:dist_ice_edge_km) %>% 
#   scale() %>% 
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val))+geom_histogram()+facet_wrap(.~var, scales = "free")

