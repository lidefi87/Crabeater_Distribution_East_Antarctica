# Loading libraries
library(tidyverse)
library(mgcv)
library(rsample)
library(recipes)
library(caret)

#Loading data
crab_ind <- read_csv("Cleaned_Data/Obs_BG_5x_Indian_weaning_LowVIF.csv") %>% 
  #Removing any rows with NA values
  drop_na()

#Selecting variables (response and predictors) to be included in the model
model_data <- crab_ind %>% 
  select(presence:dist_ice_edge_km, month) %>% 
  mutate(month = factor(month))

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

#Model scopes
scope <- list()
for(var in env_var){
  df_list <- 2:5 %>% map(\(x) paste0("s(", var, ", df = ", x, ")")) %>% unlist()
  scope[var] <- paste0("~1 + ", var, " + ", paste(df_list, collapse = " + "))
  }


scope <- list(
  "bottom_slope_deg" = ~1 + bottom_slope_deg + s(bottom_slope_deg, df = 2) + s(bottom_slope_deg, df = 3) + s(bottom_slope_deg, df = 4) + s(bottom_slope_deg, df = 5),
  "dist_shelf_km" = ~1 + dist_shelf_km + s(dist_shelf_km, df = 2) + s(dist_shelf_km, df = 3) + s(dist_shelf_km, df = 4) + s(dist_shelf_km, df = 5),
  "dist_coast_km" = ~1 + dist_coast_km + s(dist_coast_km, df = 2) + s(dist_coast_km, df = 3) + s(dist_coast_km, df = 4) + s(dist_coast_km, df = 5),
  "depth_m" = ~1 + depth_m + s(depth_m, df = 2) + s(depth_m, df = 3) + s(depth_m, df = 4) + s(depth_m, df = 5),
  "SIC" = ~1 + SIC + s(SIC, df = 2) + s(SIC, df = 3) + s(SIC, df = 4) + s(SIC, df = 5),
  "bottom_temp_degC" = ~1 + bottom_temp_degC + s(bottom_temp_degC, df = 2) + s(bottom_temp_degC, df = 3) + s(bottom_temp_degC, df = 4) + s(bottom_temp_degC, df = 5),
  "SSS_PSU" = ~1 + SSS_PSU + s(SSS_PSU, df = 2) + s(SSS_PSU, df = 3) + s(SSS_PSU, df = 4) + s(SSS_PSU, df = 5),
  "vel_lat_surf_msec" = ~1 + vel_lat_surf_msec + s(vel_lat_surf_msec, df = 2) + s(vel_lat_surf_msec, df = 3) + s(vel_lat_surf_msec, df = 4) + s(vel_lat_surf_msec, df = 5),
  "vel_lat_bottom_msec" = ~1 + vel_lat_bottom_msec + s(vel_lat_bottom_msec, df = 2) + s(vel_lat_bottom_msec, df = 3) + s(vel_lat_bottom_msec, df = 4) + s(vel_lat_bottom_msec, df = 5),
  "vel_lon_surf_msec" = ~1 + vel_lon_surf_msec + s(vel_lon_surf_msec, df = 2) + s(vel_lon_surf_msec, df = 3) + s(vel_lon_surf_msec, df = 4) + s(vel_lon_surf_msec, df = 5),
  "vel_lon_bottom_msec" = ~1 + vel_lon_bottom_msec + s(vel_lon_bottom_msec, df = 2) + s(vel_lon_bottom_msec, df = 3) + s(vel_lon_bottom_msec, df = 4) + s(vel_lon_bottom_msec, df = 5),
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

#Applying GAM
mod_gam <- gam(formula = as.formula(gam_fm), data = baked_train,
                     family = binomial(link = "cloglog"),
                     weights = train_wt, select = T,
                     method = "REML")

summary(mod_gam)


gam_fm_no_month <- paste("presence ~", paste0(env_var, collapse = " + "))
cv_mod_no_month <- train(as.formula(gam_fm_no_month), data = baked_train,
                         family = binomial(link = "cloglog"),
                         weights = train_wt, method = "gam",
                         trControl = trainControl(method = "cv", number = 10))

names_no_bottom <- model_data %>% 
  select(where(is.numeric) & !contains("bottom") & !"presence") %>% 
  names()
  
gam_fm_no_bottom <- paste("presence ~", paste0("s(", names_no_bottom[-1], ")", collapse = " + "))

mod_gam_no_bottom <- gam(formula = as.formula(gam_fm_no_bottom), data = baked_train,
                         family = binomial(link = "cloglog"),
                         weights = train_wt, select = T,
                         method = "REML")
summary(mod_gam_no_bottom)



library(gam)
gam.object <- gam(formula = presence ~s(SIC), data = baked_train,
                  family = binomial(link = "cloglog"),
                  weights = train_wt, select = T)
step.object <- step.Gam(gam.object, scope = scope, direction = "both")
step.object_for <- step.Gam(gam.object, scope = scope, direction = "forward")

gam.object.b <- gam(presence ~ s(bottom_slope_deg, df = 5) + s(dist_shelf_km, df = 5) + s(dist_coast_km, df = 5) + s(depth_m, df = 5) + s(SIC, df = 5) + s(bottom_temp_degC, df = 5) + s(SSS_PSU, df = 5) + s(vel_lat_surf_msec, df = 5) + s(vel_lat_bottom_msec, df = 5) + s(vel_lon_surf_msec, df = 5) + s(vel_lon_bottom_msec, df = 5) + s(lt_pack_ice, df = 5) + s(dist_ice_edge_km, df = 5),
                    data = baked_train,
                    family = binomial(link = "cloglog"),
                    weights = train_wt)
step.object_back <- step.Gam(gam.object.b, scope = scope, direction = "backward")


cv_mod <- train(presence ~., data = baked_train,
                family = binomial(link = "cloglog"),
                weights = train_wt, method = "gam",
                trControl = trainControl(method = "cv", number = 10))

cv_mod$finalModel
summary(cv_mod)



#Checking predictions
res <- test %>%  
  mutate(pred = predict(mod_gam, test, type = "response"))

res_back <- test %>% 
  mutate(pred = predict(step.object_back, test, type = "response")) 

pres <- res_back %>% 
  ggplot(aes(xt_ocean, yt_ocean, color = presence))+
  geom_point()

pred <- res_back %>% 
  ggplot(aes(xt_ocean, yt_ocean, color = pred))+
  geom_point()




res %>% 
  ggplot(aes(xt_ocean, pred))+
  geom_point()


# crab_ind %>% 
#   select(bottom_slope_deg:dist_ice_edge_km) %>% 
#   scale() %>% 
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val))+geom_histogram()+facet_wrap(.~var, scales = "free")

