Variable importance plots
================
Denisse Fierro Arcos
2024-03-28

- <a href="#variable-importance-by-species-distribution-model-sdms"
  id="toc-variable-importance-by-species-distribution-model-sdms">Variable
  importance by Species Distribution Model (SDMs)</a>
  - <a href="#loading-libraries" id="toc-loading-libraries">Loading
    libraries</a>
  - <a href="#loading-ggplot2-objects"
    id="toc-loading-ggplot2-objects">Loading <code>ggplot2</code>
    objects</a>
  - <a href="#bar-plot-for-models-trained-with-access-om2-01-full-set"
    id="toc-bar-plot-for-models-trained-with-access-om2-01-full-set">Bar
    plot for models trained with ACCESS-OM2-01 full set</a>
  - <a
    href="#optional-create-a-single-plot-with-sdms-trained-with-full-access-om2-01-environmental-data"
    id="toc-optional-create-a-single-plot-with-sdms-trained-with-full-access-om2-01-environmental-data">Optional:
    Create a single plot with SDMs trained with full ACCESS-OM2-01
    environmental data</a>
  - <a href="#bar-plot-for-models-trained-with-access-om2-01-reduced-set"
    id="toc-bar-plot-for-models-trained-with-access-om2-01-reduced-set">Bar
    plot for models trained with ACCESS-OM2-01 reduced set</a>
  - <a
    href="#optional-create-a-single-plot-with-sdms-trained-with-reduced-access-om2-01-environmental-data"
    id="toc-optional-create-a-single-plot-with-sdms-trained-with-reduced-access-om2-01-environmental-data">Optional:
    Create a single plot with SDMs trained with reduced ACCESS-OM2-01
    environmental data</a>
  - <a href="#bar-plot-for-models-trained-with-observations"
    id="toc-bar-plot-for-models-trained-with-observations">Bar plot for
    models trained with observations</a>
  - <a
    href="#optional-create-a-single-plot-with-sdms-trained-with-remotely-sensed-environmental-data"
    id="toc-optional-create-a-single-plot-with-sdms-trained-with-remotely-sensed-environmental-data">Optional:
    Create a single plot with SDMs trained with remotely-sensed
    environmental data</a>
  - <a href="#multipanel-option---all-results"
    id="toc-multipanel-option---all-results">Multipanel option - all
    results</a>

# Variable importance by Species Distribution Model (SDMs)

We use variable importance plots produced by the modelling scripts
included in this folder. We will combine them into a single figure for
each environmental dataset.

## Loading libraries

``` r
library(tidyverse)
library(cowplot)
```

``` r
knitr::opts_chunk$set(fig.path = "figures/") 
```

## Loading `ggplot2` objects

``` r
ggobj_list <- list.files("../../SDM_outputs/", pattern = "*_var_imp_.*rds", 
                         full.names = T, recursive = T)
```

## Bar plot for models trained with ACCESS-OM2-01 full set

``` r
var_imp_plot <- function(ggobj, mod_name, lims){
  #Access data used in plot
  p <- ggobj$data %>% 
    #arrange by variable importance
    arrange(desc(Permutation_importance)) %>% 
    # rowid_to_column("id") %>% 
    #calculate cumulative sum for importance contribution
    mutate(cum_sum = cumsum(Permutation_importance),
           #assign value of 1 if cum sum is over 50%
           fill = case_when(cum_sum > 0.5 ~ 1, T ~ 0),
           #calculate cumulative sum for this column
           fill = cumsum(fill),
           #assign blue colour if column is 0 or 1, otherwise assign grey
           fill = case_when(fill <= 1 ~ "#004488", T ~ "#bbbbbb")) %>% 
    #Initialise plot
    ggplot(aes(y = Variable, x = Permutation_importance, fill = fill))+
    #Column plot
    geom_col()+
    #use values in fill column as colours
    scale_fill_identity()+
    #Apply b&W theme
    theme_bw()+
    #Change labels
    labs(x = "Permutation importance", title = mod_name)+
    theme(axis.title.y = element_blank(), 
          plot.title = element_text(hjust = 0.5))+
    #Scales to be applied based on maximum values seen in plots
    scale_x_continuous(limits = lims, labels = scales::label_percent())
}
```

``` r
#Get list of files for ACCESS-OM2-01 full set
mod_full <- ggobj_list %>% 
  str_subset(".*mod_full.rds")

#Create empty list to store plots
mod_full_plots <- list()

#Loop through each element of the list
for(i in seq_along(mod_full)){
  #Get model name from file name
  model <- str_remove(basename(mod_full[i]), "_var.*")
  if(model == "Maxent"){
    model <- "MaxEnt"
  }
  #Load ggplot2 object
  ggobj <- readRDS(mod_full[i])
  p <- var_imp_plot(ggobj, model, c(0, .35))
  #If models are on the first row (BRT or RF), remove x axis title
  if(str_detect(model, "BRT|RF")){
    p <- p+
      labs(x = "")
  }
  #Save plots into list
  mod_full_plots[[model]] <- p
}
```

## Optional: Create a single plot with SDMs trained with full ACCESS-OM2-01 environmental data

``` r
#Turn into a single plot
mod_full_plots <- plot_grid(mod_full_plots$RF, mod_full_plots$BRT,
                            mod_full_plots$MaxEnt, mod_full_plots$GAM,
                            nrow = 2,
                            labels = c("A", "B", "C", "D"))

ggsave("../../SDM_outputs/var_imp_mod_full.png", mod_full_plots,
       device = "png", width = 9)
```

    ## Saving 9 x 5 in image

## Bar plot for models trained with ACCESS-OM2-01 reduced set

``` r
#Get list of files for ACCESS-OM2-01 reduced set
mod_match_obs <- ggobj_list %>% 
  str_subset(".*mod_match_obs.rds")

#Create empty list to store plots
mod_match_obs_plots <- list()

#Loop through each element of the list
for(i in seq_along(mod_match_obs)){
  #Get model name from file name
  model <- str_remove(basename(mod_match_obs[i]), "_var.*")
  if(model == "Maxent"){
    model <- "MaxEnt"
  }
  #Load ggplot2 object
  ggobj <- readRDS(mod_match_obs[i])
  #Create plot
  p <- var_imp_plot(ggobj, model, c(0, .4))
  
  #If models are on the first row (BRT or RF), remove x axis title
  if(str_detect(model, "BRT|RF")){
    p <- p+
      labs(x = "")
  }
  #Save plots into list
  mod_match_obs_plots[[model]] <- p
}
```

## Optional: Create a single plot with SDMs trained with reduced ACCESS-OM2-01 environmental data

``` r
#Turn into a single plot
mod_match_obs_plots <- plot_grid(mod_match_obs_plots$RF,
                                 mod_match_obs_plots$BRT,
                                 mod_match_obs_plots$MaxEnt,
                                 mod_match_obs_plots$GAM,
                                 nrow = 2,
                            labels = c("A", "B", "C", "D"))

ggsave("../../SDM_outputs/var_imp_mod_match_obs.png", mod_match_obs_plots,
       device = "png", width = 9)
```

    ## Saving 9 x 5 in image

## Bar plot for models trained with observations

``` r
#Get list of files for observations
obs <- ggobj_list %>% 
  str_subset(".*imp_obs.rds")

#Create empty list to store plots
obs_plots <- list()

#Loop through each element of the list
for(i in seq_along(obs)){
  #Get model name from file name
  model <- str_remove(basename(obs[i]), "_var.*")
  if(model == "Maxent"){
    model <- "MaxEnt"
  }
  #Load ggplot2 object
  ggobj <- readRDS(obs[i])
  #Create plot
  p <- var_imp_plot(ggobj, model, c(0, .75))
  
  #If models are on the first row (BRT or RF), remove x axis title
  if(str_detect(model, "BRT|RF")){
    p <- p+
      labs(x = "")
  }
  #Save plots into list
  obs_plots[[model]] <- p
}
```

## Optional: Create a single plot with SDMs trained with remotely-sensed environmental data

``` r
#Turn into a single plot
obs_plots <- plot_grid(obs_plots$RF, obs_plots$BRT, obs_plots$MaxEnt, 
                       obs_plots$GAM, nrow = 2, 
                       labels = c("A", "B", "C", "D"))

ggsave("../../SDM_outputs/var_imp_obs.png", obs_plots, 
       device = "png", width = 9)
```

    ## Saving 9 x 5 in image

## Multipanel option - all results

``` r
title_mod_full <- ggdraw()+
  draw_label("ACCESS-OM2-01 full set", fontface = "bold", x = 0.6)

plot_mod_full <- plot_grid(title_mod_full, 
                           plot_grid(mod_full_plots$RF, mod_full_plots$BRT, 
                                     mod_full_plots$MaxEnt+labs(x = ""),
                                     mod_full_plots$GAM, nrow = 4, 
                                     labels = "AUTO", label_x = 0.25),
                           ncol = 1, rel_heights = c(0.05, 1))

title_mod_obs <- ggdraw()+
  draw_label("ACCESS-OM2-01 reduced set", fontface = "bold", x = 0.6)

plot_mod_obs <- plot_grid(title_mod_obs, 
                          plot_grid(mod_match_obs_plots$RF, 
                                    mod_match_obs_plots$BRT, 
                                    mod_match_obs_plots$MaxEnt+labs(x = ""),
                                    mod_match_obs_plots$GAM, nrow = 4, 
                                    labels = c("E", "F", "G", "H"),
                                    label_x = 0.25),
                          ncol = 1, rel_heights = c(0.05, 1))

title_obs <- ggdraw()+
  draw_label("Remotely-sensed data", fontface = "bold", x = 0.65)

plot_obs <- plot_grid(title_obs, 
                      plot_grid(obs_plots$RF, obs_plots$BRT, 
                                obs_plots$MaxEnt+labs(x = ""), obs_plots$GAM, 
                                nrow = 4, labels = c("I", "J", "K", "L"),
                                label_x = 0.25),
                      ncol = 1, rel_heights = c(0.05, 1))

plot_grid(plot_mod_full, plot_mod_obs, plot_obs, ncol = 3)

ggsave("../../SDM_outputs/var_imp_all.pdf", height = 12, width = 11, dpi = 320, 
       bg = "white")
```
