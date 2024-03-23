Generalised Additive Model
================
Denisse Fierro Arcos
2023-09-22

- <a href="#generalised-additive-model-gam"
  id="toc-generalised-additive-model-gam">Generalised Additive Model
  (GAM)</a>
  - <a href="#loading-libraries" id="toc-loading-libraries">Loading
    libraries</a>
  - <a href="#setting-up-notebook" id="toc-setting-up-notebook">Setting up
    notebook</a>
    - <a href="#loading-mean-environmental-conditions-from-access-om2-01"
      id="toc-loading-mean-environmental-conditions-from-access-om2-01">Loading
      mean environmental conditions from ACCESS-OM2-01</a>
    - <a href="#loading-mean-environmental-conditions-from-observations"
      id="toc-loading-mean-environmental-conditions-from-observations">Loading
      mean environmental conditions from observations</a>
  - <a href="#loading-layers-for-plotting"
    id="toc-loading-layers-for-plotting">Loading layers for plotting</a>
  - <a
    href="#loading-environmental-data-from-access-om2-01-and-setting-up-variables"
    id="toc-loading-environmental-data-from-access-om2-01-and-setting-up-variables">Loading
    environmental data from ACCESS-OM2-01 and setting up variables</a>
  - <a href="#environmental-variables-matching-observations"
    id="toc-environmental-variables-matching-observations">Environmental
    variables matching observations</a>
    - <a href="#building-gam-formula" id="toc-building-gam-formula">Building
      GAM formula</a>
    - <a href="#splitting-data-into-testing-and-training"
      id="toc-splitting-data-into-testing-and-training">Splitting data into
      testing and training</a>
    - <a href="#modelling" id="toc-modelling">Modelling</a>
    - <a href="#comparing-models" id="toc-comparing-models">Comparing
      models</a>
    - <a href="#performance-metrics" id="toc-performance-metrics">Performance
      metrics</a>
    - <a href="#calculating-variable-importance"
      id="toc-calculating-variable-importance">Calculating variable
      importance</a>
    - <a href="#saving-marginal-gam-plots"
      id="toc-saving-marginal-gam-plots">Saving marginal GAM plots</a>
    - <a href="#predictions" id="toc-predictions">Predictions</a>
  - <a href="#all-environmental-variables-available-in-access-om2-01"
    id="toc-all-environmental-variables-available-in-access-om2-01">All
    environmental variables available in ACCESS-OM2-01</a>
    - <a href="#building-gam-formula-1"
      id="toc-building-gam-formula-1">Building GAM formula</a>
    - <a href="#splitting-data-into-testing-and-training-1"
      id="toc-splitting-data-into-testing-and-training-1">Splitting data into
      testing and training</a>
    - <a href="#modelling-1" id="toc-modelling-1">Modelling</a>
    - <a href="#comparing-models-1" id="toc-comparing-models-1">Comparing
      models</a>
    - <a href="#performance-metrics-1"
      id="toc-performance-metrics-1">Performance metrics</a>
    - <a href="#calculating-variable-importance-1"
      id="toc-calculating-variable-importance-1">Calculating variable
      importance</a>
    - <a href="#saving-marginal-gam-plots-1"
      id="toc-saving-marginal-gam-plots-1">Saving marginal GAM plots</a>
    - <a href="#predictions-1" id="toc-predictions-1">Predictions</a>
  - <a href="#environmental-variables-from-observations"
    id="toc-environmental-variables-from-observations">Environmental
    variables from observations</a>
    - <a href="#building-gam-formula-2"
      id="toc-building-gam-formula-2">Building GAM formula</a>
    - <a href="#splitting-data-into-testing-and-training-2"
      id="toc-splitting-data-into-testing-and-training-2">Splitting data into
      testing and training</a>
    - <a href="#modelling-2" id="toc-modelling-2">Modelling</a>
    - <a href="#comparing-models-2" id="toc-comparing-models-2">Comparing
      models</a>
    - <a href="#performance-metrics-2"
      id="toc-performance-metrics-2">Performance metrics</a>
    - <a href="#calculating-variable-importance-2"
      id="toc-calculating-variable-importance-2">Calculating variable
      importance</a>
    - <a href="#saving-marginal-gam-plots-2"
      id="toc-saving-marginal-gam-plots-2">Saving marginal GAM plots</a>
    - <a href="#predictions-2" id="toc-predictions-2">Predictions</a>
  - <a href="#differences-across-sources-of-environmental-data"
    id="toc-differences-across-sources-of-environmental-data">Differences
    across sources of environmental data</a>
- <a href="#saving-model-evaluation-results"
  id="toc-saving-model-evaluation-results">Saving model evaluation
  results</a>

# Generalised Additive Model (GAM)

GAMs are commonly used in species distribution modelling because of
their flexibility. GAMs use smoothing functions to capture non-linear
relationships between the dependent variable and the predictors (i.e.,
independent variables).

In this project, we will use GAMs as one of the models to be considered
in our Species Distribution Model ensemble to estimate the distribution
of crabeater seals in the recent past.

## Loading libraries

``` r
library(tidyverse)
library(mgcv)
library(stars)
library(cmocean)
library(cowplot)
library(prg)
library(pROC)
source("useful_functions.R")
```

## Setting up notebook

Selecting an output folder for GAM results exists and getting a list of
data files.

``` r
#Location of folder for outputs
out_folder <- "../../SDM_outputs/GAM/Mod_match_obs"
#If folder does not exist, create one
if(!dir.exists(out_folder)){
  dir.create(out_folder, recursive = T)
}

#Get path to files containing data
file_list <- list.files("../../Environmental_Data/", pattern = "Indian", full.names = T)
```

### Loading mean environmental conditions from ACCESS-OM2-01

This dataset includes the mean environmental conditions per month
(November and December) over the entire period of study (1981 to 2013).

``` r
mean_model <- read_csv("../../Environmental_Data/ACCESS-OM2-01/All_values_month_ACCESS-OM2-01_env_vars.csv") %>% 
  mutate(month = factor(month)) %>% 
  #Drop variables with high multicollinearity
  select(!c(freez_pot_Wm2, bottom_sal_psu, SIT_m))
```

    ## Rows: 730244 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (21): yt_ocean, xt_ocean, bottom_slope_deg, dist_shelf_km, dist_coast_km...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of categorical variables
cat_vars <- c("month")

mean_model_baked <- prep_pred(mean_model, cat_vars)
```

### Loading mean environmental conditions from observations

This dataset includes the mean environmental conditions per month
(November and December) over the entire period of study (1981 to 2013).

``` r
mean_obs <- read_csv("../../Environmental_Data/Env_obs/All_values_month_Obs_env_vars.csv") %>% 
  mutate(month = factor(month))
```

    ## Rows: 730244 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (11): yt_ocean, xt_ocean, bottom_slope_deg, dist_shelf_km, dist_coast_km...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of categorical variables
cat_vars <- c("month")

mean_obs_baked <- prep_pred(mean_obs, cat_vars)
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

## Loading environmental data from ACCESS-OM2-01 and setting up variables

We will use the datasets created in the notebook
`02_Merging_background_presence_data.Rmd` located within the
`Scripts/05_SDMs` folder. These datasets include the crabeater seal
observations, background points, and environmental data.

We will also define categorical and continuous explanatory variables.

## Environmental variables matching observations

First, we will look only at the variables with no multicollinearity.
This means that sea surface temperature (`SST`) is excluded even though
this variable is available in the observational dataset.

The variable `month` will be included as an ordinal factor in our
analysis.

``` r
#Loading data
mod_match_obs <- read_csv(str_subset(file_list, "match")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month))
```

    ## Rows: 32368 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (13): year, month, xt_ocean, yt_ocean, presence, bottom_slope_deg, dist_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of covariates
covars <- str_subset(names(mod_match_obs), "presence|_ocean", negate = T)
```

### Building GAM formula

We will run an initial GAM using all available covariates. Sea ice
related variables will be fitted individually to each `month` as we
expect sea ice to vary in extent between November and December because
sea ice should be retreating during this time.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)
```

### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations.

``` r
#Getting training data
mod_match_obs_split <- prep_data(mod_match_obs, cat_vars)

#Selecting model data
model_data <- mod_match_obs_split$baked_train %>%
  select(all_of(covars) | "presence")
```

### Modelling

Background data (`presence` == `0`) will be down-weighted because they
do not truly represent absences. The down-weighting applied represents
the proportion of presence to background points. In this way, the sum of
the weighted background points will be the same as the total of
presences.

We will run the GAM using all covariates available and we will
downweight background points so they have the same effect as all
presences.

``` r
#Calculating downweights
weights <- down_weights(model_data)

#Full model
full_model_gam <- gam(formula = as.formula(full_model),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(full_model_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + 
    ##     s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept) -0.410793   0.192413  -2.135   0.0328 *
    ## month12      0.007132   0.224919   0.032   0.9747  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(year)                     0.0503519      9  0.055 0.290662    
    ## s(bottom_slope_deg)         0.0002205      9  0.000 0.619763    
    ## s(dist_shelf_km)            0.3905267      9  0.652 0.160622    
    ## s(dist_coast_km)            0.0006014      9  0.000 0.347695    
    ## s(depth_m)                  1.5650818      9  7.987 0.001557 ** 
    ## s(SIC):month11              1.4633195      9  3.408 0.053837 .  
    ## s(SIC):month12              6.0886320      9 50.128  < 2e-16 ***
    ## s(lt_pack_ice):month11      2.4477928      9 13.531 0.000133 ***
    ## s(lt_pack_ice):month12      0.2393411      9  0.312 0.231695    
    ## s(dist_ice_edge_km):month11 1.8171975      9 10.113 0.000734 ***
    ## s(dist_ice_edge_km):month12 4.1482311      9 38.112  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0385   Deviance explained = 3.57%
    ## -REML = 1864.9  Scale est. = 1         n = 24276

Removing variables with effects reduced to near zero (`edf` \< 0.001) or
that are non-significant: `year`, `dist_shelf_km` (distance to
continental shelf), `bottom_slope_deg` (seafloor slope), and
`dist_coast_km` (distance to coastline in km).

``` r
#Defining new formula without non-significant covariates
significant_only <-  presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

#Applying new formula
significant_only_gam <- gam(formula = as.formula(significant_only),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

summary(significant_only_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept) -0.40952    0.19646  -2.085   0.0371 *
    ## month12      0.01067    0.22804   0.047   0.9627  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  1.4570      9  7.126 0.003461 ** 
    ## s(SIC):month11              1.5531      9  3.918 0.042621 *  
    ## s(SIC):month12              6.0721      9 50.062  < 2e-16 ***
    ## s(lt_pack_ice):month11      2.5301      9 15.019 6.87e-05 ***
    ## s(lt_pack_ice):month12      0.2687      9  0.366 0.219417    
    ## s(dist_ice_edge_km):month11 1.8054      9  9.851 0.000902 ***
    ## s(dist_ice_edge_km):month12 4.1628      9 38.512  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0381   Deviance explained = 3.54%
    ## -REML = 1864.9  Scale est. = 1         n = 24276

### Comparing models

The Akaike’s Information Criterion (AIC) will be calculated to assess
model performance. We will also check the $r^{2}$ values for both
models.

``` r
#Putting everything in a data frame
sum_gam <- data.frame(model = c("full_model", "significant_only"),
                      AIC = c(AIC(full_model_gam), AIC(significant_only_gam)),
                      rsq = c(summary(full_model_gam)$r.sq, summary(significant_only_gam)$r.sq))
#Checking results
sum_gam
```

    ##              model      AIC        rsq
    ## 1       full_model 1864.132 0.03852279
    ## 2 significant_only 1863.732 0.03814780

These models performed almost the same, so we will keep the most
parsimonious model (i.e., `significant_only`), which has a slightly
lower AIC and slightly higher $r^{2}$. We will include `SST` to this
model and check the effect of `SST` on model performance. The smoothing
for this variable will vary by month as we have done with sea ice
related variables.

``` r
significant_only_SST <-  presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, by = month)

#Applying new formula
significant_only_SST_gam <- gam(formula = as.formula(significant_only_SST),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Check results
summary(significant_only_SST_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)  -0.9736     0.4750  -2.050   0.0404 *
    ## month12       0.4931     0.4842   1.018   0.3084  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  1.4753228      9  7.144 0.003237 ** 
    ## s(SIC):month11              0.0005153      9  0.000 0.584270    
    ## s(SIC):month12              5.6627755      9 42.322  < 2e-16 ***
    ## s(lt_pack_ice):month11      2.6893674      9 17.068 2.95e-05 ***
    ## s(lt_pack_ice):month12      0.8968738      9  1.729 0.110023    
    ## s(dist_ice_edge_km):month11 1.8627920      9 11.365 0.000428 ***
    ## s(dist_ice_edge_km):month12 3.0046447      9 18.396 3.00e-05 ***
    ## s(SST_degC):month11         2.3115850      9  6.574 0.024903 *  
    ## s(SST_degC):month12         5.3547478      9 36.993  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0504   Deviance explained = 4.91%
    ## -REML = 1849.1  Scale est. = 1         n = 24276

The model that includes `SST` explains almost 5% of variability in the
data, which is higher than the previous models. We can see that now sea
ice concentration (`SIC`) for November has become insignificant, and
this is likely because of its high correlation to `SST`. We will remove
`SIC`, and run the model once more to see its effect on deviance
explained.

``` r
significant_only_no_SIC <-  presence ~ month + s(depth_m) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, by = month)

#Applying new formula
significant_only_no_SIC_gam <- gam(formula = as.formula(significant_only_no_SIC),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Check results
summary(significant_only_no_SIC_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(depth_m) + s(lt_pack_ice, by = month) + 
    ##     s(dist_ice_edge_km, by = month) + s(SST_degC, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)  -0.9836     0.4765  -2.064    0.039 *
    ## month12       0.5359     0.4796   1.117    0.264  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  0.8941993      9  8.207 0.001644 ** 
    ## s(lt_pack_ice):month11      2.7126701      9 17.437 2.50e-05 ***
    ## s(lt_pack_ice):month12      0.0002374      9  0.000 0.948834    
    ## s(dist_ice_edge_km):month11 1.8549040      9 11.199 0.000487 ***
    ## s(dist_ice_edge_km):month12 4.9083387      9 27.099 4.61e-06 ***
    ## s(SST_degC):month11         2.3205219      9  6.646 0.024312 *  
    ## s(SST_degC):month12         5.0837684      9 38.320  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0379   Deviance explained = 3.75%
    ## -REML = 1861.3  Scale est. = 1         n = 24276

We will now calculate AIC for the last two models we ran, extract their
$r^{2}$ values and add them to our summary table.

``` r
#Adding results to summary table
sum_gam <- data.frame(model = c("significant_only_SST", "significant_only_no_SIC"),
                      AIC = c(AIC(significant_only_SST_gam), AIC(significant_only_no_SIC_gam)),
                      rsq = c(summary(significant_only_SST_gam)$r.sq, summary(significant_only_no_SIC_gam)$r.sq)) %>% 
  bind_rows(sum_gam, .) %>% 
  #Arrange in increasing order by AIC
  arrange(AIC)

#Checking results
sum_gam
```

    ##                     model      AIC        rsq
    ## 1    significant_only_SST 1843.596 0.05043921
    ## 2 significant_only_no_SIC 1854.366 0.03793804
    ## 3        significant_only 1863.732 0.03814780
    ## 4              full_model 1864.132 0.03852279

We can see that the best performing model based on AIC is the last one
we tested, which includes `SST`, and `SIC`. We can check if the
performance of this model differs significantly from the model that
includes `SST`, but excludes `SIC` (`significant_only_no_SIC` in table
above). To do this, we will performance an ANOVA test.

``` r
anova(significant_only_SST_gam, significant_only_no_SIC_gam, test = "Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, 
    ##     by = month)
    ## Model 2: presence ~ month + s(depth_m) + s(lt_pack_ice, by = month) + 
    ##     s(dist_ice_edge_km, by = month) + s(SST_degC, by = month)
    ##   Resid. Df Resid. Dev      Df Deviance  Pr(>Chi)    
    ## 1     24239     3598.6                               
    ## 2     24248     3642.8 -8.8355  -44.204 1.119e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

We can see these models differ significantly, but since
multicollinearity has a negative effect on GAM predictive ability, we
will use the `significant_only_no_SIC` to predict crabeater
distribution.

### Performance metrics

To be able to compare the performance of this model with the three other
SDM algorithms to be used in the SDM ensemble, we will calculate three
metrics: area under the receiver operating curve ($AUC_{ROC}$), area
under the precisison-recall gain curve ($AUC_{PRG}$) and the Pearson
correlation between the model predictions and the testing dataset.

``` r
#Predicting values using testing dataset
pred <- predict(significant_only_no_SIC_gam, mod_match_obs_split$baked_test, 
                type = "response")

#AUC ROC
auc_roc <- roc(mod_match_obs_split$baked_test$presence, pred) %>% 
  auc() %>% 
  as.numeric()

#AUC PRG
auc_prg <- create_prg_curve(mod_match_obs_split$baked_test$presence, pred) %>% 
  calc_auprg()

#Pearson correlation
cor <- cor(pred, mod_match_obs_split$baked_test$presence)

#Save in data frame for ensemble weighting
model_eval <- data.frame(model = "GAM", env_trained = "mod_match_obs", auc_roc = auc_roc, 
           auc_prg = auc_prg, pear_cor = cor)

print(c(paste0("AUC ROC: ", round(auc_roc, 3)),
        paste0("AUC PRG: ", round(auc_prg, 3)),
        paste0("Pearson correlation: ", round(cor, 3))))
```

### Calculating variable importance

``` r
#Define variables included in best performing model
vars <- c("month", "depth_m", "lt_pack_ice", "dist_ice_edge_km", "SST_degC")

#Calculate variable importance
varimp_mod_match_obs <- compute_permutation_gam(significant_only_no_SIC_gam, 
                                            auc_roc, vars, model_data)

#Plot variable importance
p <- varimp_mod_match_obs %>% 
  plotVarImp_gam()

ggsave(file.path(out_folder, "var_import_mod_match_obs.png"), p, 
       device = "png")
```

### Saving marginal GAM plots

``` r
#Marginal plots where variables are not nested within another variable
unnested <- c("depth_m", "month")

for(v in vars){
  fname <- file.path(out_folder, paste0(v, "_marginal.png"))
  if(v %in% unnested){
    plot <- plotResponse_gam(significant_only_no_SIC_gam, model_data[,c(vars, "presence")],
               v)
  }else{
    plot <- plotResponse_gam(significant_only_no_SIC_gam, model_data[,c(vars, "presence")],
               v, nested_by = "month")
  }
  ggsave(filename = fname, plot = plot, device = "png")
  }
```

### Predictions

We will use the best performing model to predict crabeater seal
distribution using mean monthly environmental conditions obtained from
ACCESS-OM2-01.

``` r
#Prediction monthly crabeater seal distribution
mean_pred_match_obs <- mean_model_baked %>% 
  mutate(pred = as.vector(predict(significant_only_no_SIC_gam, 
                                  mean_model_baked, 
                                  type = "response")))

#Converting from data frame to raster
mean_pred_match_obs_ras <- mean_pred_match_obs %>%
  #Select relevant variables only
  select(xt_ocean, yt_ocean, pred, month) %>% 
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
mean_pred_match_obs %>% 
  write_csv(file.path(out_folder, "mean_pred_match_obs.csv"))
#Saving as R dataset so it can be easily open with readRDS
saveRDS(mean_pred_match_obs_ras,
        file.path(out_folder, "mean_pred_match_obs_raster.rds"))
```

Now that we have predicted the mean distribution of crabeaters based
with our select model, we can plot results for comparison later.

``` r
#Plotting November distribution
#Prepping data
nov <- mean_pred_match_obs_ras %>% 
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

dec <- mean_pred_match_obs_ras %>% 
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
```

    ## Warning in get_plot_component(plot, "guide-box"): Multiple components found;
    ## returning the first one. To return all, use `return_all = TRUE`.

``` r
#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_match_obs <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Mean crabeater seal distribution\n(ACCESS-OM2-01 - simplified)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
final <- plot_grid(title, plot_match_obs, ncol = 1, 
                   rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_match_obs.png"), 
       plot = final, device = "png", bg = "white", width = 8.75, height = 7)
```

## All environmental variables available in ACCESS-OM2-01

We will now use the full range of covariates available in the
ACCESS-OM2-01 model. We will test if the inclusion of additional
environmental variables improve predictive ability. We will follow the
same approach as we did before, first we will test a model with all
covariates and exclude any with non-significant effects.

We tested this dataset for multicollinearity, so we will get the names
of variables with low VIF plus SST as this was found to be an important
covariate before.

``` r
#Location of folder for outputs
out_folder <- "../../SDM_outputs/GAM/Mod_full/"
#If folder does not exist, create one
if(!dir.exists(out_folder)){
  dir.create(out_folder, recursive = T)
}

#Loading data
full_mod <- read_csv(str_subset(file_list, "model_env_pres")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month))
```

    ## Rows: 32368 Columns: 24
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (24): year, month, xt_ocean, yt_ocean, presence, bottom_slope_deg, dist_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of variables to be extracted from main dataset
mod_vars <- read_csv("../../Environmental_Data/ACCESS-OM2-01/Obs_BG_5x_Indian_weaning_LowVIF.csv") %>% 
  select(!c(sector, zone, season_year:decade)) %>% 
  names() %>% 
  append(c("SST_degC", "krill_ggp"))
```

    ## Rows: 12620 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (4): sector, zone, season_year, life_stage
    ## dbl (19): year, yt_ocean, xt_ocean, month, decade, presence, bottom_slope_de...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Extracting relevant columns only
full_mod <- full_mod %>% 
  select(all_of(mod_vars))

#List of environmental covariates to be used in model
covars <- str_subset(names(full_mod), "presence|_ocean", negate = T)
```

### Building GAM formula

We will run an initial GAM using all available covariates. Non-static
variables (i.e., variables that vary with time) will be fitted
individually to each `month` as we expect them to vary between November
and December as we move from spring to summer.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(vel_lon_bottom_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(krill_ggp, by = month)
```

### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations.

``` r
#Getting training data
full_mod_split <- prep_data(full_mod, cat_vars)

#Selecting model data
model_data <- full_mod_split$baked_train %>%
  select(all_of(covars) | "presence")
```

### Modelling

Background data will be down-weighted and then we will run the GAM using
all covariates available.

``` r
#Calculating downweights
weights <- down_weights(model_data)

#Full model
full_model_gam <- gam(formula = as.formula(full_model),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(full_model_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + 
    ##     s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, 
    ##     by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, 
    ##     by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, 
    ##     by = month) + s(vel_lon_bottom_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(krill_ggp, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.0939     0.3032  -3.608 0.000309 ***
    ## month12       0.6101     0.3108   1.963 0.049670 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(year)                        4.622e-04      9  0.000 0.331833    
    ## s(bottom_slope_deg)            9.624e-05      9  0.000 0.673685    
    ## s(dist_shelf_km)               9.285e-01      9 12.797 9.10e-05 ***
    ## s(dist_coast_km)               5.608e-05      9  0.000 0.939201    
    ## s(depth_m)                     9.431e-01      9 15.929 1.77e-05 ***
    ## s(SIC):month11                 1.544e+00      9  4.014 0.030860 *  
    ## s(SIC):month12                 2.527e+00      9  9.431 0.004894 ** 
    ## s(bottom_temp_degC):month11    2.792e+00      9 15.481 0.000158 ***
    ## s(bottom_temp_degC):month12    2.135e-04      9  0.000 0.408188    
    ## s(SSS_psu):month11             1.376e+00      9  5.047 0.010040 *  
    ## s(SSS_psu):month12             4.672e+00      9 40.776  < 2e-16 ***
    ## s(vel_lat_surf_msec):month11   4.348e+00      9 19.805 0.000157 ***
    ## s(vel_lat_surf_msec):month12   1.503e+00      9  5.926 0.013900 *  
    ## s(vel_lat_bottom_msec):month11 8.692e-01      9  6.580 0.005301 ** 
    ## s(vel_lat_bottom_msec):month12 4.252e-05      9  0.000 0.896855    
    ## s(vel_lon_surf_msec):month11   8.795e-01      9  7.018 0.003811 ** 
    ## s(vel_lon_surf_msec):month12   3.899e+00      9 12.138 0.006151 ** 
    ## s(vel_lon_bottom_msec):month11 1.851e-04      9  0.000 0.517805    
    ## s(vel_lon_bottom_msec):month12 1.239e-04      9  0.000 0.602460    
    ## s(lt_pack_ice):month11         2.805e+00      9 18.140 1.60e-05 ***
    ## s(lt_pack_ice):month12         3.041e-01      9  0.439 0.206778    
    ## s(dist_ice_edge_km):month11    2.170e+00      9 18.797 3.46e-06 ***
    ## s(dist_ice_edge_km):month12    3.822e+00      9 23.449 8.59e-06 ***
    ## s(krill_ggp):month11           3.076e+00      9 16.348 7.31e-05 ***
    ## s(krill_ggp):month12           1.425e+00      9 21.241 1.57e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0791   Deviance explained = 7.53%
    ## -REML = 1829.1  Scale est. = 1         n = 24274

This model explains more of the variability in our data than the model
that uses variables for which we also have observations. We will remove
the variables with a low contribution (`edf` $\sim 0.001$). This
includes `year`, bottom slope (`bottom_slope_deg`), distance to the
coastline (`dist_coast_km`), and bottom meridional water velocity
(`vel_lon_bottom_msec`).

We will defined this simplified model and test it below.

``` r
# Simplified model
simpler_model <- presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(krill_ggp, by = month)

#Applying new simplified formula
simpler_model_gam <- gam(formula = as.formula(simpler_model),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(simpler_model_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC, by = month) + 
    ##     s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + 
    ##     s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, 
    ##     by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(krill_ggp, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.0939     0.3032  -3.608 0.000309 ***
    ## month12       0.6101     0.3108   1.963 0.049653 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)               9.286e-01      9 12.798 9.10e-05 ***
    ## s(depth_m)                     9.424e-01      9 15.929 1.77e-05 ***
    ## s(SIC):month11                 1.544e+00      9  4.014 0.030859 *  
    ## s(SIC):month12                 2.527e+00      9  9.431 0.004892 ** 
    ## s(bottom_temp_degC):month11    2.792e+00      9 15.481 0.000158 ***
    ## s(bottom_temp_degC):month12    7.003e-05      9  0.000 0.395940    
    ## s(SSS_psu):month11             1.376e+00      9  5.047 0.010043 *  
    ## s(SSS_psu):month12             4.672e+00      9 40.775  < 2e-16 ***
    ## s(vel_lat_surf_msec):month11   4.348e+00      9 19.805 0.000157 ***
    ## s(vel_lat_surf_msec):month12   1.503e+00      9  5.926 0.013900 *  
    ## s(vel_lat_bottom_msec):month11 8.692e-01      9  6.580 0.005300 ** 
    ## s(vel_lat_bottom_msec):month12 3.000e-05      9  0.000 0.896877    
    ## s(vel_lon_surf_msec):month11   8.794e-01      9  7.018 0.003812 ** 
    ## s(vel_lon_surf_msec):month12   3.899e+00      9 12.138 0.006151 ** 
    ## s(lt_pack_ice):month11         2.805e+00      9 18.140 1.60e-05 ***
    ## s(lt_pack_ice):month12         3.041e-01      9  0.439 0.206781    
    ## s(dist_ice_edge_km):month11    2.170e+00      9 18.797 3.46e-06 ***
    ## s(dist_ice_edge_km):month12    3.822e+00      9 23.450 8.59e-06 ***
    ## s(krill_ggp):month11           3.076e+00      9 16.348 7.31e-05 ***
    ## s(krill_ggp):month12           1.425e+00      9 21.239 1.57e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0791   Deviance explained = 7.53%
    ## -REML = 1829.1  Scale est. = 1         n = 24274

We can see that this simplified model performs almost the same as the
model including all environmental variables available.

We will now explore the influence of `SST` as it was identified as an
influential variable in the previous dataset.

``` r
# Simplified model with SST
simpler_model_SST <- presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, by = month) + s(krill_ggp, by = month)

#Applying new simplified formula
simpler_model_SST_gam <- gam(formula = as.formula(simpler_model_SST),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(simpler_model_SST_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC, by = month) + 
    ##     s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + 
    ##     s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, 
    ##     by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, 
    ##     by = month) + s(krill_ggp, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.1022     0.3012  -3.659 0.000253 ***
    ## month12       0.6393     0.3134   2.040 0.041342 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)               0.8499584      9  5.560 0.006466 ** 
    ## s(depth_m)                     0.9437560      9 16.274 1.37e-05 ***
    ## s(SIC):month11                 1.6431470      9  4.552 0.023214 *  
    ## s(SIC):month12                 4.4878215      9 14.741 0.001410 ** 
    ## s(bottom_temp_degC):month11    2.7434994      9 16.088 0.000105 ***
    ## s(bottom_temp_degC):month12    0.0009083      9  0.001 0.334040    
    ## s(SSS_psu):month11             1.1969284      9  3.500 0.026545 *  
    ## s(SSS_psu):month12             3.9916333      9 14.941 0.001204 ** 
    ## s(vel_lat_surf_msec):month11   4.3845569      9 20.593 0.000105 ***
    ## s(vel_lat_surf_msec):month12   1.4628583      9  5.703 0.014790 *  
    ## s(vel_lat_bottom_msec):month11 0.8651046      9  6.351 0.006056 ** 
    ## s(vel_lat_bottom_msec):month12 0.0002504      9  0.000 0.874792    
    ## s(vel_lon_surf_msec):month11   0.8696302      9  6.413 0.005362 ** 
    ## s(vel_lon_surf_msec):month12   3.5323227      9 12.072 0.004508 ** 
    ## s(lt_pack_ice):month11         2.7748810      9 17.543 2.03e-05 ***
    ## s(lt_pack_ice):month12         0.0005160      9  0.000 0.419663    
    ## s(dist_ice_edge_km):month11    2.1386597      9 17.715 7.34e-06 ***
    ## s(dist_ice_edge_km):month12    2.8423117      9 14.861 0.000381 ***
    ## s(SST_degC):month11            0.0005140      9  0.000 0.387141    
    ## s(SST_degC):month12            5.0550330      9 26.776 2.44e-06 ***
    ## s(krill_ggp):month11           2.9820539      9 14.887 0.000160 ***
    ## s(krill_ggp):month12           0.9247427      9 12.157 0.000180 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0865   Deviance explained = 8.41%
    ## -REML = 1820.2  Scale est. = 1         n = 24274

This model including `SST` and `SIC` explains a higher proportion of the
variability in our data. We will now check the impact of `SIC` on model
performance.

``` r
# Simplified model with SST and no SIC
simpler_model_SST_noSIC <- presence ~ month + s(dist_shelf_km) + s(depth_m) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, by = month) + s(krill_ggp, by = month)

#Applying new simplified formula
simpler_model_SST_noSIC_gam <- gam(formula = as.formula(simpler_model_SST_noSIC),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(simpler_model_SST_noSIC_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(dist_shelf_km) + s(depth_m) + s(bottom_temp_degC, 
    ##     by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, 
    ##     by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, 
    ##     by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, 
    ##     by = month) + s(SST_degC, by = month) + s(krill_ggp, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.4199     0.4632  -3.065  0.00217 **
    ## month12       0.9550     0.4683   2.039  0.04142 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)               8.877e-01      9  7.714 0.001791 ** 
    ## s(depth_m)                     9.406e-01      9 15.221 1.04e-05 ***
    ## s(bottom_temp_degC):month11    2.737e+00      9 16.459 8.73e-05 ***
    ## s(bottom_temp_degC):month12    3.130e-01      9  0.395 0.242807    
    ## s(SSS_psu):month11             1.058e+00      9  2.724 0.047851 *  
    ## s(SSS_psu):month12             4.444e+00      9 41.161  < 2e-16 ***
    ## s(vel_lat_surf_msec):month11   4.377e+00      9 20.621 0.000101 ***
    ## s(vel_lat_surf_msec):month12   1.478e+00      9  5.668 0.015692 *  
    ## s(vel_lat_bottom_msec):month11 8.644e-01      9  6.319 0.006145 ** 
    ## s(vel_lat_bottom_msec):month12 5.443e-05      9  0.000 0.873036    
    ## s(vel_lon_surf_msec):month11   8.657e-01      9  6.180 0.006238 ** 
    ## s(vel_lon_surf_msec):month12   3.851e+00      9 11.778 0.007274 ** 
    ## s(lt_pack_ice):month11         2.881e+00      9 19.065 1.12e-05 ***
    ## s(lt_pack_ice):month12         1.051e-04      9  0.000 0.803357    
    ## s(dist_ice_edge_km):month11    2.165e+00      9 18.521 5.17e-06 ***
    ## s(dist_ice_edge_km):month12    3.304e+00      9 15.004 0.000583 ***
    ## s(SST_degC):month11            1.961e+00      9  5.245 0.028821 *  
    ## s(SST_degC):month12            4.312e+00      9 26.812 4.85e-06 ***
    ## s(krill_ggp):month11           3.264e+00      9 16.018 0.000153 ***
    ## s(krill_ggp):month12           1.726e+00      9 24.993  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0826   Deviance explained = 8.06%
    ## -REML = 1822.5  Scale est. = 1         n = 24274

The predictive ability of the model increased when `SIC` was excluded
and replaced by `SST`. We will calculate the Akaike Information
Criterion (AIC) for all models to compare performance across all models
fitted to the full suite of relevant environmental variables in the
ACCESS-OM2-01 model.

### Comparing models

Once again, we will calculate AIC for each model and obtain the
coefficient of determination ($r^{2}$) to assess model performance.

``` r
full_sum_gam <- data.frame(model = c("full_model", "simpler_model", "simpler_model_SST", "simpler_model_SST_noSIC"),
                      AIC = c(AIC(full_model_gam), AIC(simpler_model_gam), AIC(simpler_model_SST_gam),
                              AIC(simpler_model_SST_noSIC_gam)),
                      rsq = c(summary(full_model_gam)$r.sq, summary(simpler_model_gam)$r.sq,
                              summary(simpler_model_SST_gam)$r.sq, summary(simpler_model_SST_noSIC_gam)$r.sq))

#Arranging models from lowest to highest AIC
full_sum_gam %>% 
  arrange(AIC)
```

    ##                     model      AIC        rsq
    ## 1       simpler_model_SST 1814.849 0.08653505
    ## 2 simpler_model_SST_noSIC 1817.394 0.08255193
    ## 3           simpler_model 1827.808 0.07912673
    ## 4              full_model 1827.815 0.07912691

Even though the model including both `SST` and `SIC` is the best
performing, we will not use it to predict crabeater distribution because
including high correlated predictors can negatively affect the GAM
results. In this case, we will choose the second best performing model,
which includes `SST` only.

### Performance metrics

To be able to compare the performance of this model with the three other
SDM algorithms to be used in the SDM ensemble, we will calculate three
metrics: area under the receiver operating curve ($AUC_{ROC}$), area
under the precisison-recall gain curve ($AUC_{PRG}$) and the Pearson
correlation between the model predictions and the testing dataset.

``` r
#Predicting values using testing dataset
pred <- predict(simpler_model_SST_noSIC_gam, full_mod_split$baked_test, 
                type = "response")

#AUC ROC
auc_roc <- roc(full_mod_split$baked_test$presence, pred) %>% 
  auc() %>% 
  as.numeric()

#AUC PRG
auc_prg <- create_prg_curve(full_mod_split$baked_test$presence, pred) %>% 
  calc_auprg()

#Pearson correlation
cor <- cor(pred, full_mod_split$baked_test$presence)

#Save to data frame
model_eval <- model_eval %>% 
  bind_rows(data.frame(model = "GAM", env_trained = "full_access",
                       auc_roc = auc_roc, auc_prg = auc_prg, pear_cor = cor))

print(c(paste0("AUC ROC: ", round(auc_roc, 3)),
        paste0("AUC PRG: ", round(auc_prg, 3)),
        paste0("Pearson correlation: ", round(cor, 3))))
```

### Calculating variable importance

``` r
#Define variables included in best performing model
vars <- c("month", "dist_shelf_km", "depth_m", "bottom_temp_degC", "SSS_psu",
          "vel_lat_surf_msec", "vel_lat_bottom_msec", "vel_lon_surf_msec", 
          "lt_pack_ice", "dist_ice_edge_km", "SST_degC", "krill_ggp")

#Calculate variable importance
varimp_mod_full <- compute_permutation_gam(simpler_model_SST_noSIC_gam, 
                                            auc_roc, vars, model_data)

#Plot variable importance
p <- varimp_mod_full %>% 
  plotVarImp_gam()

ggsave(file.path(out_folder, "var_import_mod_full.png"), p, 
       device = "png")
```

### Saving marginal GAM plots

``` r
#Marginal plots where variables are not nested within another variable
unnested <- c("depth_m", "month")

for(v in vars){
  fname <- file.path(out_folder, paste0(v, "_marginal.png"))
  if(v %in% unnested){
    plot <- plotResponse_gam(simpler_model_SST_noSIC_gam, model_data[,c(vars, "presence")],
               v)
  }else{
    plot <- plotResponse_gam(simpler_model_SST_noSIC_gam, model_data[,c(vars, "presence")],
               v, nested_by = "month")
  }
  ggsave(filename = fname, plot = plot, device = "png")
}
```

### Predictions

We will use the best performing model to predict crabeater seal
distribution using mean monthly environmental conditions obtained from
ACCESS-OM2-01.

``` r
#Prediction monthly crabeater seal distribution
mean_pred_mod <- mean_model_baked %>% 
  mutate(pred = as.vector(predict(simpler_model_SST_noSIC_gam, 
                                  mean_model_baked, 
                                  type = "response")))

#Converting from data frame to raster
mean_pred_mod_ras <- mean_pred_mod %>%
  #Select relevant variables only
  select(xt_ocean, yt_ocean, pred, month) %>% 
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
mean_pred_mod %>% 
  write_csv(file.path(out_folder, "mean_pred_mod.csv"))
#Saving as R dataset so it can be easily open with readRDS
saveRDS(mean_pred_mod_ras,
        file.path(out_folder, "mean_pred_mod_raster.rds"))
```

Now that we have predicted the mean distribution of crabeaters based
with our select model, we can plot results for comparison later.

``` r
#Plotting November distribution
#Prepping data
nov <- mean_pred_mod_ras %>% 
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

dec <- mean_pred_mod_ras %>% 
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
```

    ## Warning in get_plot_component(plot, "guide-box"): Multiple components found;
    ## returning the first one. To return all, use `return_all = TRUE`.

``` r
#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_mod <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Mean crabeater seal distribution\n(ACCESS-OM2-01 - all variables)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
final <- plot_grid(title, plot_mod, ncol = 1, 
          rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_mod.png"), plot = final, 
       device = "png", bg = "white", width = 8.75, height = 7)
```

## Environmental variables from observations

Finally, we will create a model using all environmental data obtained
from observations.

``` r
#Location of folder for outputs
out_folder <- "../../SDM_outputs/GAM/Obs/"
#If folder does not exist, create one
if(!dir.exists(out_folder)){
  dir.create(out_folder, recursive = T)
}

#Loading data
obs_env_data <- read_csv(str_subset(file_list, "/obs")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month))
```

    ## Rows: 32033 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (13): year, month, xt_ocean, yt_ocean, presence, bottom_slope_deg, dist_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of covariates
covars <- str_subset(names(obs_env_data), "presence|_ocean", negate = T)
```

### Building GAM formula

We will run an initial GAM using all available covariates. Sea ice
related variables will be fitted individually to each `month` as we
expect sea ice to vary in extent between November and December because
sea ice should be retreating during this time.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(SST_degC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)
```

### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations.

``` r
#Getting training data
obs_data_split <- prep_data(obs_env_data, cat_vars)

#Selecting model data
model_data <- obs_data_split$baked_train %>%
  select(all_of(covars) | "presence")
```

### Modelling

Background data will be down-weighted and then we will run the GAM using
all covariates available.

``` r
#Calculating downweights
weights <- down_weights(model_data)

#Full model
full_model_gam <- gam(formula = as.formula(full_model),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(full_model_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + 
    ##     s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(SST_degC, 
    ##     by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.1227     0.3680  -3.051  0.00228 **
    ## month12       0.6459     0.3722   1.735  0.08271 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq p-value    
    ## s(year)                     0.5508983      9  1.006 0.15900    
    ## s(bottom_slope_deg)         0.0002183      9  0.000 0.92376    
    ## s(dist_shelf_km)            0.0002957      9  0.000 0.86746    
    ## s(dist_coast_km)            0.0002146      9  0.000 0.53618    
    ## s(depth_m)                  0.8071627      9  2.272 0.06839 .  
    ## s(SIC):month11              0.4163447      9  0.696 0.15713    
    ## s(SIC):month12              6.1731430      9 80.923 < 2e-16 ***
    ## s(SST_degC):month11         2.4877099      9  9.057 0.00612 ** 
    ## s(SST_degC):month12         4.3254166      9 13.358 0.00391 ** 
    ## s(lt_pack_ice):month11      0.0892455      9  0.093 0.22427    
    ## s(lt_pack_ice):month12      1.9954756      9 10.177 0.00215 ** 
    ## s(dist_ice_edge_km):month11 0.2689607      9  0.346 0.21222    
    ## s(dist_ice_edge_km):month12 0.7673289      9  3.249 0.03400 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.054   Deviance explained = 5.76%
    ## -REML = 1823.1  Scale est. = 1         n = 24024

This model outperforms the GAMs trained with environmental data obtained
from ACCESS-OM2-01. As we have done before, we will remove any variables
that do not contribute much to the model: bottom slope
(`bottom_slope_deg`) and distance to the continental shelf
(`dist_shelf_km`) and the coast (`dist_coast_km`).

``` r
# Simplified model
simple_obs_model <- presence ~ month + s(depth_m) + s(SIC, by = month) + s(SST_degC, by = month)+ s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

#Fitting model
simple_obs_model_gam <- gam(formula = as.formula(simple_obs_model),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(simple_obs_model_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(depth_m) + s(SIC, by = month) + s(SST_degC, 
    ##     by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.1253     0.3662  -3.073  0.00212 **
    ## month12       0.6529     0.3700   1.764  0.07766 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                edf Ref.df Chi.sq p-value    
    ## s(depth_m)                  0.8998      9  2.476 0.06435 .  
    ## s(SIC):month11              0.2977      9  0.417 0.17600    
    ## s(SIC):month12              6.1851      9 80.625 < 2e-16 ***
    ## s(SST_degC):month11         2.4889      9  9.082 0.00593 ** 
    ## s(SST_degC):month12         4.2551      9 12.824 0.00491 ** 
    ## s(lt_pack_ice):month11      0.3613      9  0.480 0.17154    
    ## s(lt_pack_ice):month12      2.0620      9 10.043 0.00268 ** 
    ## s(dist_ice_edge_km):month11 0.1339      9  0.145 0.23894    
    ## s(dist_ice_edge_km):month12 0.7716      9  3.298 0.03317 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0534   Deviance explained = 5.73%
    ## -REML = 1823.2  Scale est. = 1         n = 24024

This simplified model explains almost the same variance as the full
model.

Although multicollinearity was not identified in the observational
dataset, we will assess the effect of keeping either `SST` or `SIC` on
the predictive ability of the GAM.

``` r
# Simplified model no SST
simple_obs_model_noSST <- presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

#Fitting model
simple_obs_model_noSST_gam <- gam(formula = as.formula(simple_obs_model_noSST),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(simple_obs_model_noSST_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.60250    0.08484  -7.102 1.23e-12 ***
    ## month12      0.13757    0.09686   1.420    0.156    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df  Chi.sq p-value    
    ## s(depth_m)                  1.2878519      9   4.272  0.0307 *  
    ## s(SIC):month11              1.6085072      9   4.713  0.0323 *  
    ## s(SIC):month12              6.8157311      9 110.389  <2e-16 ***
    ## s(lt_pack_ice):month11      0.0009058      9   0.001  0.3207    
    ## s(lt_pack_ice):month12      1.9966439      9   8.551  0.0060 ** 
    ## s(dist_ice_edge_km):month11 0.2156642      9   0.257  0.2492    
    ## s(dist_ice_edge_km):month12 0.8017751      9   3.911  0.0236 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0465   Deviance explained = 5.12%
    ## -REML = 1826.1  Scale est. = 1         n = 24024

Removing `SST` from the model results in a decline in the model
predictive ability of about 0.6%. We will now check results if `SIC` is
removed.

``` r
# Simplified model no SIC
simple_obs_model_noSIC <- presence ~ month + s(depth_m) + s(SST_degC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

#Fitting model
simple_obs_model_noSIC_gam <- gam(formula = as.formula(simple_obs_model_noSIC),
    data = model_data,
    family = binomial(link = "cloglog"),
    weights = weights,
    method = "REML", select = T)

#Print summary
summary(simple_obs_model_noSIC_gam)
```

    ## 
    ## Family: binomial 
    ## Link function: cloglog 
    ## 
    ## Formula:
    ## presence ~ month + s(depth_m) + s(SST_degC, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.1392     0.3648  -3.123  0.00179 **
    ## month12       0.6418     0.3792   1.692  0.09058 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  0.7659      9  2.892  0.04438 *  
    ## s(SST_degC):month11         2.4756      9  9.085  0.00731 ** 
    ## s(SST_degC):month12         4.7429      9 17.521  0.00112 ** 
    ## s(lt_pack_ice):month11      0.4403      9  0.646  0.16528    
    ## s(lt_pack_ice):month12      1.5215      9  3.067  0.11075    
    ## s(dist_ice_edge_km):month11 0.3469      9  0.512  0.17225    
    ## s(dist_ice_edge_km):month12 4.5744      9 30.395 3.75e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0372   Deviance explained = 3.98%
    ## -REML = 1851.5  Scale est. = 1         n = 24024

The effect of removing `SIC` is much larger on the model predictive
ability. Finally, we will apply the same model fitted to the ACCESS
data.

We will not apply the same model fitted to the ACCESS data because we
did not use `SIC`.

### Comparing models

We will compare these models using AIC and $r^{2}$. We will summarise
everything in the table below.

``` r
obs_sum_gam <- data.frame(model = c("full_model", "simple_obs_model", "simple_obs_model_noSST",
                                    "simple_obs_model_noSIC"),
                      AIC = c(AIC(full_model_gam), AIC(simple_obs_model_gam), AIC(simple_obs_model_noSST_gam),
                              AIC(simple_obs_model_noSIC_gam)),
                      rsq = c(summary(full_model_gam)$r.sq, summary(simple_obs_model_gam)$r.sq,
                              summary(simple_obs_model_noSST_gam)$r.sq, summary(simple_obs_model_noSIC_gam)$r.sq))

#Rearraging data based on AIC
arrange(obs_sum_gam, AIC)
```

    ##                    model      AIC        rsq
    ## 1       simple_obs_model 1809.655 0.05344687
    ## 2             full_model 1810.246 0.05400046
    ## 3 simple_obs_model_noSST 1812.666 0.04651341
    ## 4 simple_obs_model_noSIC 1840.676 0.03718739

The simple model has the highest AIC, and its $r^{2}$, so we will use
this to predict crabeater distribution.

### Performance metrics

To be able to compare the performance of this model with the three other
SDM algorithms to be used in the SDM ensemble, we will calculate three
metrics: area under the receiver operating curve ($AUC_{ROC}$), area
under the precisison-recall gain curve ($AUC_{PRG}$) and the Pearson
correlation between the model predictions and the testing dataset.

``` r
#Predicting values using testing dataset
pred <- predict(simple_obs_model_gam, obs_data_split$baked_test, 
                type = "response")

#AUC ROC
auc_roc <- roc(obs_data_split$baked_test$presence, pred) %>% 
  auc() %>% 
  as.numeric()

#AUC PRG
auc_prg <- create_prg_curve(obs_data_split$baked_test$presence, pred) %>% 
  calc_auprg()

#Pearson correlation
cor <- cor(pred, obs_data_split$baked_test$presence)

#Save to data frame
model_eval <- model_eval %>% 
  bind_rows(data.frame(model = "GAM", env_trained = "observations", auc_roc = auc_roc, 
                       auc_prg = auc_prg, pear_cor = cor))

print(c(paste0("AUC ROC: ", round(auc_roc, 3)),
        paste0("AUC PRG: ", round(auc_prg, 3)),
        paste0("Pearson correlation: ", round(cor, 3))))
```

### Calculating variable importance

``` r
#Define variables included in best performing model
vars <- c("month", "depth_m", "SIC", "SST_degC", "lt_pack_ice",
          "dist_ice_edge_km")

#Calculate variable importance
varimp_obs <- compute_permutation_gam(simple_obs_model_gam, 
                                            auc_roc, vars, model_data)

#Plot variable importance
p <- varimp_obs %>% 
  plotVarImp_gam()

ggsave(file.path(out_folder, "var_import_obs.png"), p, 
       device = "png")
```

### Saving marginal GAM plots

``` r
#Marginal plots where variables are not nested within another variable
unnested <- c("depth_m", "month")

for(v in vars){
  fname <- file.path(out_folder, paste0(v, "_marginal.png"))
  if(v %in% unnested){
    plot <- plotResponse_gam(simple_obs_model_gam, model_data[,c(vars, "presence")],
               v)
  }else{
    plot <- plotResponse_gam(simple_obs_model_gam, model_data[,c(vars, "presence")],
               v, nested_by = "month")
  }
  ggsave(filename = fname, plot = plot, device = "png")
}
```

### Predictions

We will use the best performing model to predict crabeater seal
distribution using mean monthly environmental conditions obtained from
observations.

``` r
#Prediction monthly crabeater seal distribution
mean_pred_obs <- mean_obs_baked %>% 
  mutate(pred = as.vector(predict(simple_obs_model_gam, 
                                  mean_obs_baked, 
                                  type = "response")))

#Converting from data frame to raster
mean_pred_obs_ras <- mean_pred_obs %>%
  #Select relevant variables only
  select(xt_ocean, yt_ocean, pred, month) %>% 
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
mean_pred_obs %>% 
  write_csv(file.path(out_folder, "mean_pred_obs.csv"))
#Saving as R dataset so it can be easily open with readRDS
saveRDS(mean_pred_obs_ras,
        file.path(out_folder, "mean_pred_obs_raster.rds"))
```

Now that we have predicted the mean distribution of crabeaters based
with our select model, we can plot results for comparison later.

``` r
#Plotting November distribution
#Prepping data
nov <- mean_pred_obs_ras %>% 
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

dec <- mean_pred_obs_ras %>% 
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
```

    ## Warning in get_plot_component(plot, "guide-box"): Multiple components found;
    ## returning the first one. To return all, use `return_all = TRUE`.

``` r
#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_obs <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Mean crabeater seal distribution\n(Observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
final <- plot_grid(title, plot_obs, ncol = 1, 
                      rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_obs.png"), plot = final, 
       device = "png", bg = "white", width = 8.75, height = 7)
```

## Differences across sources of environmental data

Comparing results from ACCESS-OM2-01 (all available variables) trained
GAM with observations trained GAM.

``` r
diff <- mean_pred_mod_ras-mean_pred_obs_ras

#Plotting November distribution
#Prepping data
nov <- diff %>% 
  slice(index = 1, along = "month") 

#Plotting
nov_plot <- ggplot()+
  geom_stars(data = nov)+
  geom_sf(data = antarctica)+
  lims(x = c(0, 4000000))+
  #Set colour palette
  scale_fill_cmocean(name = "curl",
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(-1, 1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "November",
       x = "Longitude",
       y = "Latitude")

dec <- diff %>% 
  slice(index = 2, along = "month") 

dec_plot <- ggplot() +
  geom_stars(data = dec) +
  geom_sf(data = antarctica)+
  lims(x = c(0, 4000000))+
  scale_fill_cmocean(name = "curl",
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(-1, 1)) +
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
```

    ## Warning in get_plot_component(plot, "guide-box"): Multiple components found;
    ## returning the first one. To return all, use `return_all = TRUE`.

``` r
#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_obs_mod <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Differences in crabeater seal distribution\n(Full ACCESS-OM2-01 vs Observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
final <- plot_grid(title, plot_obs_mod, ncol = 1, 
                      rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_comp_pred_obs_mod.png"), plot = final, 
       device = "png", bg = "white", width = 8.75, height = 7)
```

Comparing results from ACCESS-OM2-01 (variables match observations)
trained GAM with observations trained GAM.

``` r
diff_2 <- mean_pred_match_obs_ras-mean_pred_obs_ras

#Plotting November distribution
#Prepping data
nov <- diff_2 %>% 
  slice(index = 1, along = "month") 

#Plotting
nov_plot <- ggplot()+
  geom_stars(data = nov)+
  geom_sf(data = antarctica)+
  lims(x = c(0, 4000000))+
  #Set colour palette
  scale_fill_cmocean(name = "curl",
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(-1, 1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "November",
       x = "Longitude",
       y = "Latitude")

dec <- diff_2 %>% 
  slice(index = 2, along = "month") 

dec_plot <- ggplot() +
  geom_stars(data = dec) +
  geom_sf(data = antarctica)+
  lims(x = c(0, 4000000))+
  scale_fill_cmocean(name = "curl",
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(-1, 1)) +
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
```

    ## Warning in get_plot_component(plot, "guide-box"): Multiple components found;
    ## returning the first one. To return all, use `return_all = TRUE`.

``` r
#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_obs_mod_lim <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Differences in crabeater seal distribution\n(Simplified ACCESS-OM2-01 vs Observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
final <- plot_grid(title, plot_obs_mod_lim, ncol = 1, 
                      rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_comp_pred_obs_mod_lim.png"), plot = final,
       device = "png", bg = "white", width = 8.75, height = 7)
```

# Saving model evaluation results

``` r
model_eval %>% 
  write_csv("../../SDM_outputs/model_evaluation.csv")
```
