#Useful functions to be reused in SDM 
library(tidyverse)
library(rsample)
library(recipes)
library(pROC)

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


#Defining function to prepare data for prediction
prep_pred <- function(da, cat_vars){
  
  #Creating recipe - steps to be followed when pre-processing data
  recipe <- recipe(1 ~ ., da) %>%
    #Coordinates (xt_ocean and yt_ocean) are not pre-processed
    add_role(xt_ocean, yt_ocean, new_role = "coords") %>% 
    #Categorical data
    add_role(all_of(cat_vars), new_role = "cat") %>% 
    #Scaled and center all predictors - exclude coordinates
    step_center(all_predictors(), -c(has_role("coords"), has_role("cat"))) %>% 
    step_scale(all_predictors(), -c(has_role("coords"), has_role("cat")))
  
    #Applying recipe to training data
    out <- prep(recipe, training = da) %>% 
      bake(new_data = da)
  
  #Return list as a result
  return(out)
}


# Defining function to calculate downweights
down_weights <- function(da){
  weight <- da %>% 
    #Count presences and background points
    count(presence) %>% 
    pivot_wider(names_from = presence, values_from = n) %>% 
    #Calculate downweight by dividing total presences by background points
    mutate(weight = `1`/`0`) %>% 
    pull(weight)
  
  #Creating a vector with weights to be applied to presence column
  weights <- model_data %>% 
    mutate(weight = case_when(presence == 0 ~ weight,
                              T ~ presence)) %>% 
    pull(weight)
  
  #Return downweights
  return(weights)
}

#Formatting data for SDM tune
sdm_format <- function(data){
  coords <- data %>% 
    arrange(desc(presence)) %>% 
    select(xt_ocean, yt_ocean) %>% 
    rename("x" = "xt_ocean", "y" = "yt_ocean") %>% 
    as.data.frame()
  sdmt_bg <- SWD(species = "Crabeater seals",
                 coords = coords,
                 data = data %>% 
                   select(-c(presence, xt_ocean, yt_ocean)) %>%
                   as.data.frame(),
                 pa = data$presence)
  return(sdmt_bg)
}


# Defining function to calculate variable importance for GAM
# This function comes from the SDMTune package and it has been adapted to work
# with GAMs
compute_permutation_gam <- function(model, model_auc, vars, data_origin, 
                                #Use same default value as SDMTune
                                permut = 10){
  #Create matrix to store results
  permuted_auc <- matrix(nrow = permut, ncol = length(vars))
  #For reproducibility
  set.seed(25)
  
  #Loop through variables
  for (j in seq_along(vars)) {
    #Loop through number of permutations
    for (i in seq_len(permut)) {
      #Resample variable and replace in original data
      data <- sample(pull(data_origin[, vars[j]]))
      train_copy <- data_origin
      train_copy[, vars[j]] <- data
      
      #Calculate AUC ROC from resampled data
      #Predict with resampled data
      pred <- predict(model, train_copy, type = "response")
      #Calculate AUC ROC
      auc_roc <- pROC::roc(data_origin$presence, pred) %>% 
        pROC::auc() %>% 
        as.numeric()
      #Add results to matrix
      permuted_auc[i, j] <- auc_roc
    }}
  
  #If multiple permutations request calculated SD and mean of AUC
  if (permut > 1) {
    sd_auc <- apply(permuted_auc, 2, sd)
    permuted_auc <- apply(permuted_auc, 2, mean)
  }
  
  #Calculate variable importance
  perm_imp <- pmax(0, (model_auc - permuted_auc))
  perm_imp <- 100 * perm_imp / sum(perm_imp)
  perm_imp <- round(perm_imp, 1)
  
  #Save results as data frame
  output <- data.frame(Variable = vars, Permutation_importance = perm_imp,
                       stringsAsFactors = FALSE)
  if (permut > 1)
    output$sd <- round(sd_auc, 3)
  
  output
}


# Defining function to plot variable importance
# This function comes from the SDMTune package to standardise all plots
plotVarImp_gam <- function(df, color = "grey"){
  df <- df[order(df[, 2]), ]
  df[, 2] <- df[, 2] / 100
  df[, 1] <- factor(df[, 1], levels = df[, 1])
  y_name <- colnames(df)[2]
  
  ggplot(df, aes(x = .data$Variable, y = .data[[y_name]])) +
    labs(x = "", y = sub("_", " ", y_name)) +
    scale_y_continuous(labels = scales::percent) +
    geom_bar(position = "dodge",
                      stat = "identity",
                      fill = color) +
    coord_flip() +
    theme_minimal() +
    theme(text = element_text(colour = "#666666"))
}


# Defining function to plot marginal responses in GAM
# This function comes from the SDMTune packages and has been adapted to
# plot GAM results
plotResponse_gam <- function(model,
                         data_origin,
                         var, nested_by = NULL,
                         type = "response",
                         only_presence = T,
                         fun = mean,
                         rug = T,
                         color = "#4bc0c0") {
  
  #Dividing into presence and absence
  p <- data_origin[data_origin$presence == 1,]
  a <- data_origin[data_origin$presence == 0,]
  
  #Getting presence data only
  df <- p
  
  #Identifying categorical and continuous variables
  cont_vars <- names(Filter(is.numeric, df))
  cat_vars <- names(Filter(is.factor, df))
  
  p_rug <- data.frame(x = p[[var]])
  a_rug <- data.frame(x = a[[var]])
  
  #Loop through variables
  if (var %in% cat_vars) {
    categ <- as.numeric(levels(df[[var]]))
    n_rows <- length(categ)
  } else if (!is.null(nested_by)) {
    nest_cat <- as.numeric(levels(df[[nested_by]]))
    n_rows <- 200
    p_rug$nested_by <- p[[nested_by]]
    a_rug$nested_by <- a[[nested_by]]
  } else {
    n_rows <- 100
  }
  
  #Define function to get data to plot marginal responses
  get_plot_data <- function(model, data_origin, var, nested_by = NULL, 
                            nest_cat = NULL, df, cont_vars, cat_vars, n_rows, 
                            fun, categ) {
    
    data <- data.frame(matrix(NA, nrow = 1, ncol = ncol(df)))
    colnames(data) <- colnames(df)
    data[cont_vars] <- apply(df[cont_vars], 2, fun)
    if (is.null(nested_by)){
      data[cat_vars] <- as.factor(names(which.max(table(df[[cat_vars]]))))
    } else {
      data[nested_by] <- as.factor(nest_cat)
    }
    data <- do.call("rbind", replicate(n_rows, data, simplify = FALSE))
    
    if (var %in% cont_vars) {
      var_min <- min(data_origin[[var]], na.rm = T)
      var_max <- max(data_origin[[var]], na.rm = T)
      data[var] <- seq(var_min, var_max, length.out = n_rows)
      for (c in cat_vars) {
        levels(data[, c]) <- levels(df[[c]])
      }
    } else {
      data[var] <- factor(categ)
    }
    
    pred <- predict(model, data, type = "response")
    
    data.frame(x = data[[var]],
               y = pred)
  }
  
  #Getting data to plot marginal responses
  if (!is.null(nested_by)){
    plot_data <- data.frame()
    for (n in nest_cat){
      p_data <- get_plot_data(model, data_origin[data_origin$month == n,], 
                              var, nested_by, n, df,
                              cont_vars, cat_vars, n_rows, fun, categ)
      p_data$nested_by <- n
      plot_data <- rbind(plot_data, p_data)
    }
    my_plot <- ggplot(plot_data, aes(x = .data$x, y = .data$y)) +
      geom_line(colour = color)+
      facet_grid(~nested_by, scales = "free_y")
  } else {
    plot_data <- get_plot_data(model, data_origin, var, nested_by = NULL, 
                               nest_cat = NULL, df, cont_vars, cat_vars, n_rows,
                               fun, categ)
    
    if (var %in% cont_vars) {
      my_plot <- ggplot(plot_data, aes(x = .data$x, y = .data$y)) +
        geom_line(colour = color)
      
    } else {
      my_plot <- ggplot(plot_data, aes(x = .data$x, y = .data$y)) +
        geom_bar(stat = "identity", fill = color)
    }
  }
  
  my_plot <- my_plot +
    labs(x = var, y = "Probability of presence") +
    theme_minimal() +
    theme(text = element_text(colour = "#666666"))
  
  if (rug == TRUE & var %in% cont_vars) {
    my_plot <- my_plot +
      geom_rug(data = p_rug, inherit.aes = FALSE, aes(.data$x),
               sides = "t", color = "#4C4C4C") +
      geom_rug(data = a_rug, inherit.aes = FALSE, aes(.data$x),
               sides = "b", color = "#4C4C4C")
  }
  
  my_plot
  
}


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


confMatrix_adap_ensemble <- function(data,
                             pred) {
  
  n_p <- sum(data$presence == 1)
  n_a <- sum(data$presence == 0)
  p_pred <- pred[1:n_p]
  a_pred <- pred[(n_p + 1):(n_p + n_a)]
  
  th <- sort(unique(pred))
  th <- c(0, th, 1)
  
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


thresholds_adap_ensemble <- function(data,
                             pred) {
  
  n_pres <- sum(data$presence == 1)
  
  cm_train <- confMatrix_adap_ensemble(data, pred)
  tpr <- cm_train$tp / (cm_train$tp + cm_train$fn)
  tnr <- cm_train$tn / (cm_train$fp + cm_train$tn)
  fpr <- cm_train$fp / (cm_train$fp + cm_train$tn)
  
  mtp <- min(pred)
  ess <- cm_train$th[which.min(abs(tpr - tnr))]
  mss <- cm_train$th[which.max(tpr + tnr)]
  
  ths <- c(mtp, ess, mss)
  rownames <- c("Minimum training presence",
                "Equal training sensitivity and specificity",
                "Maximum training sensitivity plus specificity")
  colnames <- c("Threshold",
                "Value",
                "Fractional predicted area",
                "Training omission rate")
  
  or_train <- vector(mode = "numeric", length = length(ths))
  fpa <- vector(mode = "numeric", length = length(ths))
  
  for (i in seq_along(ths)) {
    index <- which.min(abs(cm_train$th - ths[i]))
    or_train[i] <- cm_train[index, ]$fn / n_pres
    fpa[i] <- fpr[index]
  }
  
  output <- data.frame(th = rownames, val = ths, fpa = fpa, or = or_train,
                       stringsAsFactors = FALSE)
  
  colnames(output) <- colnames
  
  output
}




