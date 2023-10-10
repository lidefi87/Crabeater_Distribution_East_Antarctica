#Useful functions to be reused in SDM 
library(tidyverse)
library(rsample)
library(recipes)

#Defining function to prepare data for training and testing
prep_data <- function(da, cat_vars, split = T){
  
  #Creating recipe - steps to be followed when pre-processing data
  recipe <- recipe(presence ~ ., da) %>%
    #Coordinates (xt_ocean and yt_ocean) are not pre-processed
    add_role(xt_ocean, yt_ocean, new_role = "coords") %>% 
    #Categorical data
    add_role(all_of(cat_vars), new_role = "cat") %>% 
    #Scaled and center all predictors - exclude coordinates
    step_center(all_predictors(), -c(has_role("coords"), has_role("cat"), "year")) %>% 
    step_scale(all_predictors(), -c(has_role("coords"), has_role("cat"), "year"))
  
  #Dividing dataset (training and testing) - If "split" set to True
  if(split){
    #Setting seed for reproducible results
    set.seed(42)
    #Splitting into testing and training
    split <- initial_split(da, prop = 0.75, strata = "presence") 
    train <- training(split) %>% 
      arrange(desc(presence))
    test <- testing(split) %>% 
      arrange(desc(presence))
    #Applying recipe to training data
    baked_train <- prep(recipe, training = train) %>% 
      bake(new_data = train)
    #Applying recipe to testing data
    baked_test <- prep(recipe, new_data = test) %>% 
      bake(new_data = test)
    #Creating a list with training and testing data
    out <- list(baked_train = baked_train,
                baked_test = baked_test)
  }else{
    #Applying recipe to training data
    out <- prep(recipe, training = da) %>% 
      bake(new_data = da)
  }
  
  #Return list as a result
  return(out)
}



# library(caret)
# library(vip)
# set.seed(42)
# mod1 <- train(y ~ x, data = data, method = 'lm', 
#              trControl = trainControl(method = "cv", number = 10))
# 
# #To check variable importance
# vip(mod1, num_features = '#', method = "model")
# 
# #Produces RMSE - how far (on average is an estimate from the actual value when model is applied to unseen data)
# 
# #Extract out of sample performance measures
# summary(resamples(list(model1 = mod1,
#                        model2 = mod2, ....)))





