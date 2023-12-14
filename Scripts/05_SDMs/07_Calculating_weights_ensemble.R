###############################################################################
# Calculating weights to create SDM ensemble mean
#
# Author: Denisse Fierro Arcos
# Date: 2023-12-14
# 
# Since not all SDMs performed equally well, the final ensemble will be 
# calculated as a weighted mean of the four SDMs used in this project: GAMs, 
# Maxent, BRTs and RFs. The weighting will depend on three model performance 
# metrics: area under the the receiver operating curve (AUC ROC), area under the
# precision-recall gain curve (AUC PRG) and the Pearson correlation between the
# model predictions and the testing dataset.

# Loading libraries -------------------------------------------------------
library(tidyverse)


# Loading model evaluation results ----------------------------------------
model_eval_path <- "SDM_outputs/model_evaluation.csv"
model_eval <- read_csv(model_eval_path) 


# Calculating weights -----------------------------------------------------
model_eval <- model_eval %>% 
  #Weights are based on the sum of three metrics
  rowwise() %>% 
  mutate(sum_metric_model = sum(across(auc_roc:pear_cor))) %>% 
  #Add values per source of environmental data
  group_by(env_trained) %>% 
  mutate(metric_group = sum(sum_metric_model),
         weights = sum_metric_model/metric_group)


# Saving weights ----------------------------------------------------------
model_eval %>% 
  write_csv("SDM_outputs/model_evaluation.csv")
  
