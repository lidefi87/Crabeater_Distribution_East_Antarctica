Maxent trained with observations
================
Denisse Fierro Arcos
2023-11-02

- <a href="#maxent-via-sdmtune" id="toc-maxent-via-sdmtune">Maxent via
  SDMtune</a>
  - <a href="#loading-libraries" id="toc-loading-libraries">Loading
    libraries</a>
  - <a href="#setting-up-notebook" id="toc-setting-up-notebook">Setting up
    notebook</a>
  - <a href="#loading-environmental-data-from-observations"
    id="toc-loading-environmental-data-from-observations">Loading
    environmental data from observations</a>
    - <a href="#splitting-data-into-testing-and-training"
      id="toc-splitting-data-into-testing-and-training">Splitting data into
      testing and training</a>
  - <a href="#loading-mean-environmental-conditions-from-observations"
    id="toc-loading-mean-environmental-conditions-from-observations">Loading
    mean environmental conditions from observations</a>
  - <a href="#loading-layers-for-plotting"
    id="toc-loading-layers-for-plotting">Loading layers for plotting</a>
  - <a href="#modelling" id="toc-modelling">Modelling</a>
  - <a href="#variable-importance" id="toc-variable-importance">Variable
    importance</a>
  - <a href="#jacknife-tests" id="toc-jacknife-tests">Jacknife tests</a>
    - <a href="#plotting-jacknife-results"
      id="toc-plotting-jacknife-results">Plotting Jacknife results</a>
  - <a href="#auc-curves" id="toc-auc-curves">AUC curves</a>
  - <a href="#true-skill-statistic-tss"
    id="toc-true-skill-statistic-tss">True Skill Statistic (TSS)</a>
  - <a href="#model-report" id="toc-model-report">Model report</a>
  - <a href="#reducing-model-variables"
    id="toc-reducing-model-variables">Reducing model variables</a>
  - <a href="#training-model-with-reduced-variables"
    id="toc-training-model-with-reduced-variables">Training model with
    reduced variables</a>
  - <a href="#variable-importance-1" id="toc-variable-importance-1">Variable
    importance</a>
  - <a href="#model-report-1" id="toc-model-report-1">Model report</a>
  - <a href="#auc-curves-1" id="toc-auc-curves-1">AUC curves</a>
  - <a href="#true-skill-statistic-tss-1"
    id="toc-true-skill-statistic-tss-1">True Skill Statistic (TSS)</a>
  - <a href="#performance-metrics" id="toc-performance-metrics">Performance
    metrics</a>
  - <a href="#predictions" id="toc-predictions">Predictions</a>
  - <a href="#plotting-predictions" id="toc-plotting-predictions">Plotting
    predictions</a>

# Maxent via SDMtune

MaxEnt is one of the most widely used species distribution model
algorithm.

In this project, we will use MaxEnt as one of the models to be
considered in our Species Distribution Model ensemble to estimate the
distribution of crabeater seals in the recent past.

## Loading libraries

``` r
library(tidyverse)
library(SDMtune)
library(stars)
library(sf)
library(cmocean)
library(cowplot)
library(prg)
source("useful_functions.R")
```

## Setting up notebook

Selecting an output folder for GAM results exists and getting a list of
data files.

``` r
#Location of folder for outputs
out_folder <- "../../SDM_outputs/Maxent/Obs"
#If folder does not exist, create one
if(!dir.exists(out_folder)){
  dir.create(out_folder, recursive = T)
}

#Get path to files containing data
file_list <- list.files("../../Environmental_Data/", pattern = "Indian", 
                        full.names = T)
```

## Loading environmental data from observations

We will use the datasets created in the notebook
`02_Merging_background_presence_data.Rmd` located within the
`Scripts/05_SDMs` folder. These datasets include the crabeater seal
observations, background points, and environmental data.

We will also define categorical and continuous explanatory variables.

``` r
#Loading data
obs_data <- read_csv(str_subset(file_list, "/obs.*VIF")) %>% 
  select(!c(sector, zone, season_year:decade)) %>% 
  #Setting month as factor and ordered factor
  mutate(month = as.factor(month)) %>% 
  drop_na()
```

    ## Rows: 32512 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (4): sector, zone, season_year, life_stage
    ## dbl (13): year, yt_ocean, xt_ocean, month, decade, presence, bottom_slope_de...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
covars <- str_subset(names(obs_data), "presence|_ocean", negate = T)
```

### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations. We will
then transform the data into SWD (“samples with data”) format, which is
the required format for inputs used in the `SDMtune` library.

``` r
#Getting training data
obs_data <- prep_data(obs_data, "month", split = F)

#Applying SWD format to model data
model_data <- obs_data %>% 
  select(!year) %>% 
  sdm_format() %>% 
  trainValTest(test = 0.25, only_presence = T, seed = 42)
```

## Loading mean environmental conditions from observations

This dataset includes the mean environmental conditions per month
(November and December) over the entire period of study (1981 to 2013).

``` r
mean_obs <- read_csv("../../Environmental_Data/Env_obs/All_values_month_Obs_env_vars.csv") %>% 
  mutate(month = as.factor(month)) %>% 
  #Drop variables with high multicollinearity
  select(ends_with("_ocean")|any_of(covars))
```

    ## Rows: 730244 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (11): yt_ocean, xt_ocean, bottom_slope_deg, dist_shelf_km, dist_coast_km...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
mean_obs_baked <- prep_pred(mean_obs, "month")
```

## Loading layers for plotting

We will extract this layer from the `rnaturalearth` package. We will
then reproject this layer to South Polar Stereographic (`EPSG 3976`).

``` r
#Loading layer
antarctica <- rnaturalearth::ne_countries(continent = "Antarctica",
                                          returnclass = "sf") %>% 
  #Transforming to South Polar Stereographic
  st_transform(3976)
```

## Modelling

MaxEnt has different feature classes (`fc`, otherwise known as
restrictions) available for modelling. These `fc` include:  
- `l` - lineal,  
- `q` - quadratic,  
- `p` - product,  
- `t` - threshold,  
- `h` - hinge  
and any possible combination of these 5 features.

Regularisation (`reg`) refers to *L1 regularisation* also known as
*Lasso (Least Absolute Shrinkage and Selection Operator) regression*.
This involves adding an absolute value of magnitude as a penalty term to
the loss function. It is used to prevent overfitting. In MaxEnt a `reg`
value lower than 1 results in a outputs that fit closer to presence
data. The risk of using values that are too small is a model that
overfits and therefore does not generalised well. While, `reg` values
larger than 1 result in less localised predictions, producing smoother
or more diffuse distributions.

Here, we use the `SDMtune` library to test various value combinations
for regularisation, feature classes and number of iterations. We will
identify the “best model” using the `AUC` for the testing dataset.

``` r
#Train model
default_model <- train(method = "Maxent", data = model_data[[1]])

# Define the hyperparameters to test
hyp_parm <- list(reg = seq(0.5, 5, 0.5),
                 #Feature classes
                 fc = c("lq", "lh", "lqp", "lqph", "lqpht"),
                 #Number of iterations
                 iter = c(500, 1000, 1500))

# Test all the possible combinations with gridSearch
gs_obs <- optimizeModel(default_model, hypers = hyp_parm, metric = "auc", 
                        test = model_data[[2]], seed = 42)

#Check best performing models based on AUC
gs_obs@results %>% 
  #Adding index as column to identify best model easily
  rownames_to_column("index") %>% 
  #Arranging results by AUC from testing data (descending order)
  arrange(-test_AUC) %>% 
  #Showing only the top 5 models
  head(n = 5)

#Best model based on test AUC and smallest AUC difference between train and test
best_max_obs <- gs_obs@models[[1]]

out <- file.path(out_folder, "initial_Maxent_model")
#If folder does not exist, create one
if(!dir.exists(out)){
  dir.create(out, recursive = T)
}

#Save best performing model to disk to avoid having to tune model again
best_max_obs %>% 
  saveRDS(file.path(out, "best_maxent_obs_grid.rds"))
```

## Variable importance

We will check the contribution of each environmental variable to the
model performance.

``` r
#Calculating importance
var_imp_best <- varImp(best_max_obs) 
```

    ## Variable importance  ■■■■■                             12% | ETA: 48s - 00:00:6…Variable importance  ■■■■■■■■■                         25% | ETA: 38s - 00:00:1…Variable importance  ■■■■■■■■■■■■                      38% | ETA: 31s - 00:00:1…Variable importance  ■■■■■■■■■■■■■■■■                  50% | ETA: 25s - 00:00:2…Variable importance  ■■■■■■■■■■■■■■■■■■■■              62% | ETA: 17s - 00:00:2…Variable importance  ■■■■■■■■■■■■■■■■■■■■■■■           75% | ETA: 11s - 00:00:3…Variable importance  ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% | ETA:  5s - 00:00:3…Variable importance  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s - 00:00:4…

``` r
#Plotting results
p <- var_imp_best %>% 
  plotVarImp()

p
```

![](04c_MaxEnt_obs_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Sea ice concentration (`SIC`) is the variable with the highest
importance in this model, and its about twice as important as the sea
surface temperature (`SST_degC`), which was the second most important
variable.

Both `month` and the slope of the seafloor (`bottom_slope_deg`)
contribute less than 5% each to the model, so we may be able to drop
these without affecting overall model performance.

## Jacknife tests

We can now check which variables contributed the most to this model.

``` r
jk_obs <- doJk(best_max_obs, metric = "auc", test = model_data[[2]])
```

    ## Loading required namespace: rJava

    ## Jk Test  ■■■                                6% | ETA: 11m - 00:00:43.4Jk Test  ■■■■■                             12% | ETA:  5m - 00:00:43.9Jk Test  ■■■■■■■                           19% | ETA:  6m - 00:01:27.2Jk Test  ■■■■■■■■■                         25% | ETA:  4m - 00:01:27.9Jk Test  ■■■■■■■■■■                        31% | ETA:  5m - 00:02:3.4 Jk Test  ■■■■■■■■■■■■                      38% | ETA:  3m - 00:02:4.2Jk Test  ■■■■■■■■■■■■■■                    44% | ETA:  3m - 00:02:40 Jk Test  ■■■■■■■■■■■■■■■■                  50% | ETA:  3m - 00:02:41.1Jk Test  ■■■■■■■■■■■■■■■■■■                56% | ETA:  2m - 00:03:1.3 Jk Test  ■■■■■■■■■■■■■■■■■■■■              62% | ETA:  2m - 00:03:2.9Jk Test  ■■■■■■■■■■■■■■■■■■■■■■            69% | ETA:  2m - 00:03:26.2Jk Test  ■■■■■■■■■■■■■■■■■■■■■■■           75% | ETA:  1m - 00:03:27.6Jk Test  ■■■■■■■■■■■■■■■■■■■■■■■■■         81% | ETA:  1m - 00:04:0.2 Jk Test  ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% | ETA: 34s - 00:04:0.8Jk Test  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■     94% | ETA: 18s - 00:04:31 Jk Test  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s - 00:04:33.8

``` r
jk_obs
```

    ##           Variable Train_AUC_without Train_AUC_withonly Test_AUC_without
    ## 1            month         0.7345890          0.5151006        0.6889295
    ## 2 bottom_slope_deg         0.7330105          0.5260331        0.6933127
    ## 3    dist_coast_km         0.7221313          0.5609909        0.6785276
    ## 4          depth_m         0.7219998          0.5725372        0.6822852
    ## 5              SIC         0.7048972          0.6117493        0.6557080
    ## 6         SST_degC         0.7125643          0.6008730        0.6756056
    ## 7      lt_pack_ice         0.7294782          0.5624657        0.6855787
    ## 8 dist_ice_edge_km         0.7127247          0.6159980        0.6705409
    ##   Test_AUC_withonly
    ## 1         0.4988169
    ## 2         0.4781427
    ## 3         0.5351563
    ## 4         0.5293081
    ## 5         0.6199659
    ## 6         0.5813533
    ## 7         0.5497532
    ## 8         0.5960851

### Plotting Jacknife results

``` r
plotJk(jk_obs, type = "train", ref = SDMtune::auc(best_max_obs))
```

![](04c_MaxEnt_obs_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

When used on their own, `SIC`, `SST_degC` and distance to the sea ice
edge (`dist_ice_edge_km`) contribute the most to model performance.
Their removal affects performance negatively. While `month` and
`bottom_slope_deg` contribute the least and their removal has almost no
impact.

``` r
plotJk(jk_obs, type = "test", ref = SDMtune::auc(best_max_obs, 
                                                 test = model_data[[2]]))
```

![](04c_MaxEnt_obs_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Results from the training dataset are pretty similar. We can check if we
are able to remove `bottom_slope_deg` and `month` from the model without
affecting performance.

## AUC curves

We will calculate AUC curves, so we can compare to the simplified models
we will test.

``` r
plotROC(best_max_obs, test = model_data[[2]])
```

    ## Warning: The following aesthetics were dropped during statistical transformation: m and
    ## d.
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in
    ##   the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
    ##   variable into a factor?

![](04c_MaxEnt_obs_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## True Skill Statistic (TSS)

``` r
tss(best_max_obs)
```

    ## [1] 0.3516164

This model is not a good performer.

## Model report

Before moving onto testing a new model, we will save a report with the
information shown above.

``` r
modelReport(best_max_obs, type = "cloglog", 
            folder = file.path(out_folder, "initial_Maxent_model"),
            test = model_data[[2]], response_curves = T, only_presence = T, 
            jk = T)
```

## Reducing model variables

As discussed above, there are some variables that do not contribute much
to the model performance, such as `month` and the slope of the seafloor
(`bottom_slope_deg`). We will check here if we are able to exclude them
without negatively impacting performance.

We will use 10% as a threshold, but we will only exclude variables if
they do not decrease model performance.

``` r
simple_model <- reduceVar(best_max_obs, metric = "auc", test = model_data[[2]],
                          th = 10, permut = 10, use_jk = T)
```

    ## ✔ The variables bottom_slope_deg and month have been removed

``` r
simple_model
```

    ## 

    ## ── Object of class: <SDMmodel> ──

    ## 

    ## Method: Maxent

    ## 

    ## ── Hyperparameters

    ## • fc: "lqpht"

    ## • reg: 1

    ## • iter: 1500

    ## 

    ## ── Info

    ## • Species: Crabeater seals

    ## • Presence locations: 1366

    ## • Absence locations: 29790

    ## 

    ## ── Variables

    ## • Continuous: "dist_coast_km", "depth_m", "SIC", "SST_degC", "lt_pack_ice", and
    ## "dist_ice_edge_km"

    ## • Categorical: NA

As suspected, `month` and the slope of the seafloor (`bottom_slope_deg`)
can be excluded from our model.Now, we will fitted the model again to
find the best hyperparameter values.

## Training model with reduced variables

``` r
#Select variables
model_data_simple <- obs_data %>% 
  select(!c(year, month, bottom_slope_deg)) %>% 
  sdm_format() %>% 
  trainValTest(test = 0.25, only_presence = T, seed = 42)

#Train model
default_model <- train(method = "Maxent", data = model_data_simple[[1]])

# Test all the possible hyper parameter combinations with gridSearch
gs_mod <- optimizeModel(default_model, hypers = hyp_parm, metric = "auc", 
                        test = model_data_simple[[2]], seed = 42)

#Check best performing models based on AUC
gs_mod@results %>% 
  #Adding index as column to identify best model easily
  rownames_to_column("index") %>% 
  #Arranging results by AUC from testing data (descending order)
  arrange(-test_AUC) %>% 
  #Showing only the top 5 models
  head(n = 5)

#Best model based on test AUC and smallest AUC difference between train and test
best_max_mod <- gs_mod@models[[1]]

out <- file.path(out_folder, "reduced_Maxent_model")
#If folder does not exist, create one
if(!dir.exists(out)){
  dir.create(out, recursive = T)
}

best_max_mod %>% 
  saveRDS(file.path(out, "best_red_maxent_model.rds"))
```

The only hyperparameter changed was the number of iterations. We will
save a report for this model.

## Variable importance

``` r
#Calculating variable contribution based on permutations
var_imp_best <- varImp(best_max_mod) 

#Plotting results
p <- var_imp_best %>% 
  plotVarImp()

saveRDS(p, "../../SDM_outputs/Maxent/Maxent_var_imp_obs.rds")
```

![](04c_MaxEnt_obs_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

## Model report

``` r
#Produce report
modelReport(best_max_mod, type = "cloglog", folder = out, 
            test = model_data_simple[[2]], 
            response_curves = T, only_presence = T, jk = T)
```

## AUC curves

We will calculate AUC curves, so we can compare to the simplified models
we will test.

``` r
plotROC(best_max_mod, test = model_data_simple[[2]])
```

    ## Warning: The following aesthetics were dropped during statistical transformation: m and
    ## d.
    ## ℹ This can happen when ggplot fails to infer the correct grouping structure in
    ##   the data.
    ## ℹ Did you forget to specify a `group` aesthetic or to convert a numerical
    ##   variable into a factor?

![](04c_MaxEnt_obs_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

## True Skill Statistic (TSS)

``` r
tss(best_max_mod)
```

    ## [1] 0.356092

## Performance metrics

To be able to compare the performance of this model with the three other
SDM algorithms to be used in the SDM ensemble, we will calculate three
metrics: area under the receiver operating curve ($AUC_{ROC}$), area
under the precisison-recall gain curve ($AUC_{PRG}$) and the Pearson
correlation between the model predictions and the testing dataset.

``` r
#Predicting values using testing dataset
pred <- predict(best_max_mod, model_data_simple[[2]]@data, type = "cloglog")

#AUC ROC
auc_roc <- SDMtune::auc(best_max_mod, model_data_simple[[2]])

#AUC PRG
auc_prg <- create_prg_curve(model_data_simple[[2]]@pa, pred) %>% 
  calc_auprg()

#Pearson correlation
cor <- cor(pred, model_data_simple[[2]]@pa)

print(c(paste0("AUC ROC: ", round(auc_roc, 3)),
        paste0("AUC PRG: ", round(auc_prg, 3)),
        paste0("Pearson correlation: ", round(cor, 3))))
```

Saving model evaluation results.

``` r
#Load model evaluation data frame and add results
model_eval_path <- "../../SDM_outputs/model_evaluation.csv"
read_csv(model_eval_path) %>% 
  bind_rows(data.frame(model = "Maxent", env_trained = "observations", 
                       auc_roc = auc_roc, auc_prg = auc_prg, 
                       pear_cor = cor)) %>% 
  write_csv(model_eval_path)
```

## Predictions

We will use this reduced Maxent model to predict the distribution of
crabeater seals.

``` r
pred_obs <- mean_obs_baked %>% 
  drop_na() %>% 
  mutate(pred = as.vector(predict(best_max_mod,
                                  data = drop_na(mean_obs_baked),
                                  type = "cloglog")))

pred_obs_ras <- pred_obs %>% 
  #Select relevant variables only
  select(xt_ocean, yt_ocean, pred, month) %>% 
  right_join(mean_obs_baked %>%
  #Select relevant variables only
  select(xt_ocean, yt_ocean, month)) %>% 
  #Set dimensions
  st_as_stars(dims = c("xt_ocean", "yt_ocean", "month")) %>% 
  #Ensuring month dimension is shown correctly
  st_set_dimensions("month", values = c(11, 12)) %>%
  #Set CRS
  st_set_crs(4326) %>% 
  #Transform to South Pole stereographic
  st_transform(crs = st_crs(3976))

#Saving outputs
#Data frame
pred_obs %>% 
  write_csv(file.path(out, "mean_pred_obs.csv"))
#Saving as R dataset so it can be easily open with readRDS
saveRDS(pred_obs_ras,
        file.path(out, "mean_pred_obs_raster.rds"))
```

## Plotting predictions

``` r
#Plotting November distribution
#Prepping data
nov <- pred_obs_ras %>% 
  slice(index = 1, along = "month") 

#Plotting
nov_plot <- ggplot()+
  geom_stars(data = nov)+
  geom_sf(data = antarctica)+
  lims(x = c(0, 4000000))+
  #Set colour palette
  scale_fill_cmocean(name = "haline", direction = -1, 
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(0, 1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "November",
       x = "Longitude",
       y = "Latitude")

dec <- pred_obs_ras %>% 
  slice(index = 2, along = "month") 

dec_plot <- ggplot() +
  geom_stars(data = dec) +
  geom_sf(data = antarctica)+
  lims(x = c(0, 4000000))+
  scale_fill_cmocean(name = "haline", direction = -1, 
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(0, 1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "December",
       x = "Longitude",
       y = " ",
       fill = "Probability")

#Get legend
legend <- get_legend(dec_plot)

#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_obs <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Mean crabeater seal distribution\n(observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
final <- plot_grid(title, plot_obs, ncol = 1, 
          rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_obs.png"), plot = final, 
       device = "png", bg = "white", width = 8.75, height = 7)
```
