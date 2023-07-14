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
  select(presence:dist_ice_edge_km)

#Calculating weights to be applied to presence column
wt <- model_data %>% 
  group_by(presence) %>% 
  count() %>% 
  pivot_wider(names_from = presence, values_from = n) %>% 
  mutate(weight = `1`/`0`) %>% 
  pull(weight)

#Getting a list of prediction variables
env_var <- names(model_data)[-1]

#Building formulas for GAM
#Full model
gam_fm <- paste("presence ~", paste(paste0("s(", env_var, ")"), collapse = " + "))

#Model scopes
scope <- list()
for(var in env_var){
  scope[var] <- paste0("~1 + ", var, " + poly(", var, ", 2)")
}

#Splitting data - Training and testing
set.seed(42)
#Ensuring training and test data keep similar obs/background ratios
split <- initial_split(model_data, prop = 0.75, strata = "presence")
train <- training(split)
test <- testing(split)
#Getting weights for train and test sets
train_wt <- if_else(train$presence == 0, wt, 1)

#Preparing train data 
baked_train <- recipe(presence ~., data = train) %>% 
  #Standardising data
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) %>% 
  #Applying recipe to training dataset
  prep(train) %>% 
  bake(train)

#Applying GAM
mod_gam <- gam(formula = as.formula(gam_fm), data = baked_train,
                     family = binomial(link = "cloglog"),
                     weights = train_wt,
                     method = "REML")


cv_mod <- train(presence ~., data = baked_train,
                family = binomial(link = "cloglog"),
                weights = train_wt, method = "gam",
                trControl = trainControl(method = "cv", number = 10))



#Checking predictions
res <- test %>% 
  select(presence) %>% 
  mutate(pred = predict(mod_gam, test, type = "response"))

res %>% 
  ggplot(aes(presence, pred))+
  geom_point()


# crab_ind %>% 
#   select(bottom_slope_deg:dist_ice_edge_km) %>% 
#   scale() %>% 
#   as.data.frame() %>% 
#   pivot_longer(everything(), names_to = "var", values_to = "val") %>% 
#   ggplot(aes(x = val))+geom_histogram()+facet_wrap(.~var, scales = "free")

