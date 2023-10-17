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
    - <a href="#all-environmental-variables-available-in-access-om2-01"
      id="toc-all-environmental-variables-available-in-access-om2-01">All
      environmental variables available in ACCESS-OM2-01</a>
    - <a href="#environmental-variables-from-observations"
      id="toc-environmental-variables-from-observations">Environmental
      variables from observations</a>
    - <a href="#differences" id="toc-differences">Differences</a>

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
source("useful_functions.R")
```

## Setting up notebook

Selecting an output folder for GAM results exists and getting a list of
data files.

``` r
#Location of folder for outputs
out_folder <- "../../SDM_outputs/GAM"
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
  mutate(month = factor(month),
         month_ordered = factor(month, ordered = T)) %>% 
  #Drop variables with high multicollinearity
  select(!c(freez_pot_Wm2, bottom_sal_psu, SIT_m))
```

    ## Rows: 1147094 Columns: 20
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (20): yt_ocean, xt_ocean, bottom_slope_deg, dist_shelf_km, dist_coast_km...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of categorical variables
cat_vars <- c("month", "month_ordered")

mean_model_baked <- prep_pred(mean_model, cat_vars)
```

### Loading mean environmental conditions from observations

This dataset includes the mean environmental conditions per month
(November and December) over the entire period of study (1981 to 2013).

``` r
mean_obs <- read_csv("../../Environmental_Data/Env_obs/All_values_month_Obs_env_vars.csv") %>% 
  mutate(month = factor(month),
         month_ordered = factor(month, ordered = T))
```

    ## Rows: 1147094 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (11): yt_ocean, xt_ocean, bottom_slope_deg, dist_shelf_km, dist_coast_km...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of categorical variables
cat_vars <- c("month", "month_ordered")

mean_obs_baked <- prep_pred(mean_obs, cat_vars)
```

## Loading layers for plotting

We will extract this layer from the `rnaturalearth` package. We will
then reproject this layer to South Polar Stereographic (`EPSG 3976`).

We will also load the MEASO layers

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

### Environmental variables matching observations

First, we will look only at the variables with no multicollinearity.
This means that sea surface temperature (`SST`) is excluded even though
this variable is available in the observational dataset.

The variable `month` will be included as an ordinal factor in our
analysis.

``` r
#Loading data
mod_match_obs <- read_csv(str_subset(file_list, "match")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month),
         month_ordered = factor(month, ordered = T))
```

    ## Rows: 37368 Columns: 13
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

#### Building GAM formula

We will run an initial GAM using all available covariates. Sea ice
related variables will be fitted individually to each `month` as we
expect sea ice to vary in extent between November and December because
sea ice should be retreating during this time.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)
```

#### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations.

``` r
#Getting training data
mod_match_obs_split <- prep_data(mod_match_obs, cat_vars)

#Selecting model data
model_data <- mod_match_obs_split$baked_train %>%
  select(all_of(covars) | "presence")
```

#### Modelling

Background data (`presence` == `0`) will be down-weighted because they
do not truly represent absences. The down-weighting applied represents
the proportion of presence to background points. In this way, the sum of
the weighted background points will be the same as the total of
presences.

We will run the GAM using all covariates available and we will
downweight background points so they have the same effect as all
presences.

**Other option: Infinitely weighted - 1e3^(1-y) -
`model_data %>% mutate(weight = 1e3^(1-presence)) %>% pull(weight)`**

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
    ##     s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + 
    ##     s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -0.8133     0.1492  -5.451    5e-08 ***
    ## month12       0.3245     0.1740   1.865   0.0621 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                           edf Ref.df Chi.sq  p-value    
    ## s(year)                             0.0003323      9  0.000   0.6675    
    ## s(bottom_slope_deg)                 0.0001508      9  0.000   0.8322    
    ## s(dist_shelf_km)                    0.6792296      9  2.196   0.0611 .  
    ## s(dist_coast_km)                    0.0011885      9  0.001   0.3811    
    ## s(depth_m)                          1.8942476      9 18.315 5.58e-06 ***
    ## s(SIC)                              3.9369584      9 30.247  < 2e-16 ***
    ## s(SIC):month_ordered12              0.4496818      9  0.560   0.1493    
    ## s(lt_pack_ice)                      4.2689087      9 61.826  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12      4.4632891      9 41.325  < 2e-16 ***
    ## s(dist_ice_edge_km)                 4.2677893      9 22.422 1.49e-06 ***
    ## s(dist_ice_edge_km):month_ordered12 4.6017066      9 35.097  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0691   Deviance explained = 7.12%
    ## -REML = 1854.9  Scale est. = 1         n = 28026

Removing variables with effects reduced to near zero (`edf` \< 0.001)
and that are non-significant: `year`, `bottom_slope_deg` (seafloor
slope), `dist_coast_km` (distance to coast in km). After evaluating the
diagnostic plots from the model above, the `dist_shelf_km` (distance to
continental shelf in km) will also be removed because its effect is
negligible.

``` r
#Defining new formula without non-significant covariates
significant_only <-  presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)

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
    ## presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC) + s(SIC, 
    ##     by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + 
    ##     s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -0.8132     0.1491  -5.453 4.97e-08 ***
    ## month12       0.3244     0.1739   1.865   0.0622 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                        edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)                    0.6792      9  2.196   0.0611 .  
    ## s(depth_m)                          1.8943      9 18.316 5.59e-06 ***
    ## s(SIC)                              3.9374      9 30.293  < 2e-16 ***
    ## s(SIC):month_ordered12              0.4482      9  0.557   0.1496    
    ## s(lt_pack_ice)                      4.2689      9 61.824  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12      4.4632      9 41.324  < 2e-16 ***
    ## s(dist_ice_edge_km)                 4.2691      9 22.435 1.46e-06 ***
    ## s(dist_ice_edge_km):month_ordered12 4.6023      9 35.110  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0691   Deviance explained = 7.12%
    ## -REML = 1854.9  Scale est. = 1         n = 28026

#### Comparing models

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
    ## 1       full_model 1842.895 0.06907069
    ## 2 significant_only 1842.886 0.06907159

These models performed almost the same, so we will keep the most
parsimonious model (i.e., `significant_only`), which has a slightly
lower AIC and slightly higher $r^{2}$. We will include `SST` to this
model and check the effect of `SST` on model performance. The smoothing
for this variable will vary by month as we have done with sea ice
related variables.

``` r
significant_only_SST <-  presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered) + s(SST_degC) + s(SST_degC, by = month_ordered)

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
    ## presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC) + s(SIC, 
    ##     by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + 
    ##     s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered) + 
    ##     s(SST_degC) + s(SST_degC, by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.9192     0.4279  -4.486 7.27e-06 ***
    ## month12       1.3767     0.4362   3.156   0.0016 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                          edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)                    0.001346      9  0.001 0.335948    
    ## s(depth_m)                          1.804959      9 12.590 0.000282 ***
    ## s(SIC)                              0.002490      9  0.001 0.635751    
    ## s(SIC):month_ordered12              4.645646      9 44.465  < 2e-16 ***
    ## s(lt_pack_ice)                      4.404155      9 68.578  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12      4.583419      9 48.521  < 2e-16 ***
    ## s(dist_ice_edge_km)                 2.360633      9  7.907 0.003508 ** 
    ## s(dist_ice_edge_km):month_ordered12 2.989646      9 10.829 0.000747 ***
    ## s(SST_degC)                         4.971576      9 21.818 3.24e-06 ***
    ## s(SST_degC):month_ordered12         3.449215      9 28.365  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0836   Deviance explained = 8.89%
    ## -REML = 1831.6  Scale est. = 1         n = 28026

The model that includes `SST` explains almost 9% of variability in the
data, which is higher than the previous models. We can see that now sea
ice concentration (`SIC`) has become insignificant, and this is likely
because of its high correlation to `SST`. We will keep the `SIC` by
month because it is highly significant. We will remove this variable,
together with distance to the shelf (`dist_shelf_km`) because its effect
is near zero (`edf` \< 0.001) and it is also non-significant.

``` r
significant_only_no_SIC <-  presence ~ month + s(depth_m) + s(SIC, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered) + s(SST_degC) + s(SST_degC, by = month_ordered)

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
    ## presence ~ month + s(depth_m) + s(SIC, by = month_ordered) + 
    ##     s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered) + s(SST_degC) + s(SST_degC, 
    ##     by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.9221     0.4285  -4.486 7.27e-06 ***
    ## month12       1.3714     0.4367   3.141  0.00169 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                       edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                          1.788      9 12.580 0.000278 ***
    ## s(SIC):month_ordered12              4.734      9 44.937  < 2e-16 ***
    ## s(lt_pack_ice)                      4.405      9 68.610  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12      4.582      9 48.419  < 2e-16 ***
    ## s(dist_ice_edge_km)                 2.351      9  7.845 0.003611 ** 
    ## s(dist_ice_edge_km):month_ordered12 2.982      9 10.803 0.000748 ***
    ## s(SST_degC)                         4.975      9 21.869 3.23e-06 ***
    ## s(SST_degC):month_ordered12         3.451      9 28.323  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0837   Deviance explained = 8.89%
    ## -REML = 1831.7  Scale est. = 1         n = 28026

The performance of this reduced model is basically the same as the
previous one. We will calculate AIC and extract $r^{2}$ values and add
them to our summary table.

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
    ## 1 significant_only_no_SIC 1810.308 0.08371041
    ## 2    significant_only_SST 1810.402 0.08361118
    ## 3        significant_only 1842.886 0.06907159
    ## 4              full_model 1842.895 0.06907069

We can see that the best performing model based on AIC is the last one
we tested, which includes `SST`, excludes `SIC`, but keeps `SIC` by
month (`significant_only_SST` in table above). We can check if the
performance of this model differs significantly from the model that did
not include `SST` (`significant_only` in table above). To do this, we
will performance an ANOVA test.

``` r
anova(significant_only_gam, significant_only_no_SIC_gam, test = "Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: presence ~ month + s(dist_shelf_km) + s(depth_m) + s(SIC) + s(SIC, 
    ##     by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + 
    ##     s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)
    ## Model 2: presence ~ month + s(depth_m) + s(SIC, by = month_ordered) + 
    ##     s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered) + s(SST_degC) + s(SST_degC, 
    ##     by = month_ordered)
    ##   Resid. Df Resid. Dev     Df Deviance  Pr(>Chi)    
    ## 1     27988     3600.2                              
    ## 2     27982     3531.3 5.7619   68.833 5.163e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

We can see these models differ significantly, so we will use the
`significant_only_no_SIC` to predict crabeater distribution.

#### Predictions

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
  lims(x = c(0, 5200000))+
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
  lims(x = c(0, 5200000))+
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
plot_match_obs <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Mean crabeater seal distribution \n(ACCESS-OM2-01 - simplified)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
plot_match_obs <- plot_grid(title, plot_match_obs, ncol = 1, 
          rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_match_obs.png"), 
       plot = plot_match_obs, device = "png", bg = "white")
```

    ## Saving 7 x 5 in image

### All environmental variables available in ACCESS-OM2-01

We will now use the full range of covariates available in the
ACCESS-OM2-01 model. We will test if the inclusion of additional
environmental variables improve predictive ability. We will follow the
same approach as we did before, first we will test a model with all
covariates and exclude any with non-significant effects.

We tested this dataset for multicollinearity, so we will get the names
of variables with low VIF plus SST as this was found to be an important
covariate before.

``` r
#Loading data
full_mod <- read_csv(str_subset(file_list, "model_env_pres")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month),
         month_ordered = factor(month, ordered = T))
```

    ## Rows: 37368 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## dbl (22): year, month, xt_ocean, yt_ocean, presence, bottom_slope_deg, dist_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#List of variables to be extracted from main dataset
mod_vars <- read_csv("../../Environmental_Data/ACCESS-OM2-01/Obs_BG_5x_Indian_weaning_LowVIF.csv") %>% 
  select(!c(sector, zone, season_year:decade)) %>% 
  names() %>% 
  append(c("SST_degC", "month_ordered"))
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

#### Building GAM formula

We will run an initial GAM using all available covariates. Non-static
variables (i.e., variables that vary with time) will be fitted
individually to each `month` as we expect them to vary between November
and December as we move from spring to summer.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + s(bottom_temp_degC) + s(bottom_temp_degC, by = month_ordered) + s(SSS_psu) + s(SSS_psu, by = month_ordered) + s(vel_lat_surf_msec) + s(vel_lat_surf_msec, by = month_ordered) + s(vel_lat_bottom_msec) + s(vel_lat_bottom_msec, by = month_ordered) + s(vel_lon_surf_msec) + s(vel_lon_surf_msec, by = month_ordered) + s(vel_lon_bottom_msec) + s(vel_lon_bottom_msec, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)
```

#### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations.

``` r
#Getting training data
full_mod_split <- prep_data(full_mod, cat_vars)

#Selecting model data
model_data <- full_mod_split$baked_train %>%
  select(all_of(covars) | "presence")
```

#### Modelling

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
    ##     s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + 
    ##     s(bottom_temp_degC) + s(bottom_temp_degC, by = month_ordered) + 
    ##     s(SSS_psu) + s(SSS_psu, by = month_ordered) + s(vel_lat_surf_msec) + 
    ##     s(vel_lat_surf_msec, by = month_ordered) + s(vel_lat_bottom_msec) + 
    ##     s(vel_lat_bottom_msec, by = month_ordered) + s(vel_lon_surf_msec) + 
    ##     s(vel_lon_surf_msec, by = month_ordered) + s(vel_lon_bottom_msec) + 
    ##     s(vel_lon_bottom_msec, by = month_ordered) + s(lt_pack_ice) + 
    ##     s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept) -0.53327    0.17613  -3.028  0.00246 **
    ## month12     -0.08614    0.20135  -0.428  0.66881   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                              edf Ref.df Chi.sq  p-value    
    ## s(year)                                0.0002998      9  0.000 0.555282    
    ## s(bottom_slope_deg)                    0.0001709      9  0.000 0.853350    
    ## s(dist_shelf_km)                       0.8924647      9  8.495 0.001430 ** 
    ## s(dist_coast_km)                       0.8206653      9  4.278 0.017782 *  
    ## s(depth_m)                             0.9316954      9 13.389 8.45e-05 ***
    ## s(SIC)                                 2.8867299      9 23.813 2.72e-06 ***
    ## s(SIC):month_ordered12                 0.0002795      9  0.000 0.501249    
    ## s(bottom_temp_degC)                    5.5951040      9 55.640  < 2e-16 ***
    ## s(bottom_temp_degC):month_ordered12    6.5146620      9 36.023  < 2e-16 ***
    ## s(SSS_psu)                             3.1855825      9 24.250 1.72e-07 ***
    ## s(SSS_psu):month_ordered12             3.2378244      9 10.696 0.000883 ***
    ## s(vel_lat_surf_msec)                   0.0001996      9  0.000 0.482792    
    ## s(vel_lat_surf_msec):month_ordered12   0.7693990      9  3.314 0.032146 *  
    ## s(vel_lat_bottom_msec)                 0.4867581      9  0.741 0.191371    
    ## s(vel_lat_bottom_msec):month_ordered12 0.0001682      9  0.000 0.937063    
    ## s(vel_lon_surf_msec)                   1.9425535      9  6.614 0.014929 *  
    ## s(vel_lon_surf_msec):month_ordered12   0.0008855      9  0.001 0.254727    
    ## s(vel_lon_bottom_msec)                 0.0013089      9  0.001 0.459148    
    ## s(vel_lon_bottom_msec):month_ordered12 0.3832751      9  0.466 0.257706    
    ## s(lt_pack_ice)                         4.9292826      9 97.035  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12         4.7911611      9 67.106  < 2e-16 ***
    ## s(dist_ice_edge_km)                    3.6727886      9 18.855 1.84e-05 ***
    ## s(dist_ice_edge_km):month_ordered12    4.0911788      9 22.366 4.01e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.113   Deviance explained = 11.7%
    ## -REML = 1812.4  Scale est. = 1         n = 28026

This model explains more of the variability in our data than the model
that uses variables for which we also have observations. We will remove
the variables with a low contribution (`edf` $\sim 0.001$). This
includes `year`, bottom slope (`bottom_slope_deg`), zonal water velocity
at the surface (`vel_lat_surf_msec`, but we will keep monthly fits),
meridional water velocity at the bottom (`vel_lon_bottom_msec`), and the
monthly fits for meridional water velocity at the
surface(`vel_lon_surf_msec`), bottom zonal velocity
(`vel_lat_bottom_msec`) and sea ice concentration (`SIC`).

We will defined this simplified model and test it below.

``` r
# Simplified model
simpler_model <- presence ~ s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC) + s(bottom_temp_degC) + s(bottom_temp_degC, by = month_ordered) + s(SSS_psu) + s(SSS_psu, by = month_ordered) + s(vel_lat_surf_msec, by = month_ordered) + s(vel_lat_bottom_msec) + s(vel_lon_surf_msec) + s(vel_lon_bottom_msec, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)

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
    ## presence ~ s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + 
    ##     s(SIC) + s(bottom_temp_degC) + s(bottom_temp_degC, by = month_ordered) + 
    ##     s(SSS_psu) + s(SSS_psu, by = month_ordered) + s(vel_lat_surf_msec, 
    ##     by = month_ordered) + s(vel_lat_bottom_msec) + s(vel_lon_surf_msec) + 
    ##     s(vel_lon_bottom_msec, by = month_ordered) + s(lt_pack_ice) + 
    ##     s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.60265    0.04785  -12.59   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                           edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)                       0.8948      9  8.820 0.001220 ** 
    ## s(dist_coast_km)                       0.8148      9  4.259 0.017773 *  
    ## s(depth_m)                             0.9341      9 14.023 5.68e-05 ***
    ## s(SIC)                                 2.9551      9 24.012 1.82e-06 ***
    ## s(bottom_temp_degC)                    5.4848      9 57.236  < 2e-16 ***
    ## s(bottom_temp_degC):month_ordered12    6.5867      9 37.607  < 2e-16 ***
    ## s(SSS_psu)                             3.1416      9 27.697  < 2e-16 ***
    ## s(SSS_psu):month_ordered12             3.5537      9 12.233 0.000567 ***
    ## s(vel_lat_surf_msec):month_ordered12   0.7569      9  3.098 0.036603 *  
    ## s(vel_lat_bottom_msec)                 0.5440      9  0.881 0.173941    
    ## s(vel_lon_surf_msec)                   1.9797      9  6.954 0.012373 *  
    ## s(vel_lon_bottom_msec):month_ordered12 0.4995      9  0.645 0.240264    
    ## s(lt_pack_ice)                         4.9322      9 96.440  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12         4.8062      9 69.418  < 2e-16 ***
    ## s(dist_ice_edge_km)                    3.7425      9 17.683 2.85e-05 ***
    ## s(dist_ice_edge_km):month_ordered12    4.1345      9 22.556 3.31e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.114   Deviance explained = 11.8%
    ## -REML = 1810.3  Scale est. = 1         n = 28026

We can see that this simplified model performs slightly better than the
original model including all environmental variables available.

We will now explore the influence of `SST` as it was identified as an
influential variable in the previous dataset.

``` r
# Simplified model with SST
simpler_model_SST <- presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC) + s(bottom_temp_degC) + s(bottom_temp_degC, by = month_ordered) + s(SSS_psu) + s(SSS_psu, by = month_ordered) + s(vel_lat_surf_msec, by = month_ordered) + s(vel_lat_bottom_msec) + s(vel_lon_surf_msec) + s(vel_lon_bottom_msec, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered) + s(SST_degC) + s(SST_degC, by = month_ordered)

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
    ## presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + 
    ##     s(SIC) + s(bottom_temp_degC) + s(bottom_temp_degC, by = month_ordered) + 
    ##     s(SSS_psu) + s(SSS_psu, by = month_ordered) + s(vel_lat_surf_msec, 
    ##     by = month_ordered) + s(vel_lat_bottom_msec) + s(vel_lon_surf_msec) + 
    ##     s(vel_lon_bottom_msec, by = month_ordered) + s(lt_pack_ice) + 
    ##     s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered) + s(SST_degC) + s(SST_degC, 
    ##     by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.0480     0.3054  -3.432 0.000599 ***
    ## month12       0.4266     0.3192   1.337 0.181343    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                              edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)                       7.388e-01      9  3.060 0.034010 *  
    ## s(dist_coast_km)                       6.701e-01      9  1.965 0.071773 .  
    ## s(depth_m)                             9.481e-01      9 18.203 3.28e-06 ***
    ## s(SIC)                                 3.821e+00      9 29.850  < 2e-16 ***
    ## s(bottom_temp_degC)                    4.268e+00      9 26.008 7.74e-06 ***
    ## s(bottom_temp_degC):month_ordered12    2.572e-05      9  0.000 0.647321    
    ## s(SSS_psu)                             3.675e+00      9 27.799  < 2e-16 ***
    ## s(SSS_psu):month_ordered12             2.334e+00      9 15.790 1.31e-05 ***
    ## s(vel_lat_surf_msec):month_ordered12   8.284e-01      9  4.796 0.013075 *  
    ## s(vel_lat_bottom_msec)                 3.286e-01      9  0.434 0.230833    
    ## s(vel_lon_surf_msec)                   2.142e+00      9  8.781 0.004946 ** 
    ## s(vel_lon_bottom_msec):month_ordered12 1.587e+00      9  3.005 0.117013    
    ## s(lt_pack_ice)                         4.785e+00      9 91.466  < 2e-16 ***
    ## s(lt_pack_ice):month_ordered12         4.674e+00      9 58.917  < 2e-16 ***
    ## s(dist_ice_edge_km)                    2.481e+00      9 11.184 0.000294 ***
    ## s(dist_ice_edge_km):month_ordered12    3.379e+00      9 12.282 0.000300 ***
    ## s(SST_degC)                            5.822e-01      9  1.232 0.092953 .  
    ## s(SST_degC):month_ordered12            4.428e+00      9 35.804  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.108   Deviance explained = 11.3%
    ## -REML =   1807  Scale est. = 1         n = 28026

This model is slightly worse than when `SST` was not included. This is
likely because `SST` was highly correlated to `SIC`. Therefore, we will
not use `SST` to predict crabeater seal distribution.

#### Comparing models

Once again, we will calculate AIC for each model and obtain the
coefficient of determination ($r^{2}$) to assess model performance.

``` r
full_sum_gam <- data.frame(model = c("full_model", "simpler_model", "simpler_model_SST"),
                      AIC = c(AIC(full_model_gam), AIC(simpler_model_gam), AIC(simpler_model_SST_gam)),
                      rsq = c(summary(full_model_gam)$r.sq, summary(simpler_model_gam)$r.sq, summary(simpler_model_SST_gam)$r.sq))

#Arranging models from lowest to highest AIC
full_sum_gam %>% 
  arrange(AIC)
```

    ##               model      AIC       rsq
    ## 1     simpler_model 1781.605 0.1142566
    ## 2 simpler_model_SST 1784.066 0.1075902
    ## 3        full_model 1786.276 0.1126668

Based on AIC and $r^{2}$, we can see that the simplified model is the
best performing of the three models tested here.

#### Predictions

We will use the best performing model to predict crabeater seal
distribution using mean monthly environmental conditions obtained from
ACCESS-OM2-01.

``` r
#Prediction monthly crabeater seal distribution
mean_pred_mod <- mean_model_baked %>% 
  mutate(pred = as.vector(predict(simpler_model_gam, 
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
  lims(x = c(0, 5200000))+
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
  lims(x = c(0, 5200000))+
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
plot_mod <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Mean crabeater seal distribution \n(ACCESS-OM2-01 - all variables)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
plot_mod <- plot_grid(title, plot_mod, ncol = 1, 
          rel_heights = c(0.05, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_mod.png"), 
       plot = plot_mod, device = "png", bg = "white")
```

    ## Saving 7 x 5 in image

### Environmental variables from observations

Finally, we will create a model using all environmental data obtained
from observations.

``` r
#Loading data
obs_env_data <- read_csv(str_subset(file_list, "/obs")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month),
         month_ordered = factor(month, ordered = T))
```

    ## Rows: 37296 Columns: 13
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

#### Building GAM formula

We will run an initial GAM using all available covariates. Sea ice
related variables will be fitted individually to each `month` as we
expect sea ice to vary in extent between November and December because
sea ice should be retreating during this time.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + s(SST_degC) + s(SST_degC, by = month_ordered) + s(lt_pack_ice) + s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)
```

#### Splitting data into testing and training

The `prep_data` function in the `useful_functions` script will be used
to split our data and to apply all necessary transformations.

``` r
#Getting training data
obs_data_split <- prep_data(obs_env_data, cat_vars)

#Selecting model data
model_data <- obs_data_split$baked_train %>%
  select(all_of(covars) | "presence")
```

#### Modelling

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
    ##     s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + 
    ##     s(SST_degC) + s(SST_degC, by = month_ordered) + s(lt_pack_ice) + 
    ##     s(lt_pack_ice, by = month_ordered) + s(dist_ice_edge_km) + 
    ##     s(dist_ice_edge_km, by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -2.0957     0.3634  -5.767 8.05e-09 ***
    ## month12       1.5077     0.3683   4.094 4.24e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                           edf Ref.df Chi.sq  p-value    
    ## s(year)                             3.062e-05      9  0.000  0.86584    
    ## s(bottom_slope_deg)                 2.216e-05      9  0.000  0.64848    
    ## s(dist_shelf_km)                    7.598e-01      9  3.252  0.03189 *  
    ## s(dist_coast_km)                    8.985e-01      9  8.964  0.00128 ** 
    ## s(depth_m)                          2.122e+00      9 11.002  0.00112 ** 
    ## s(SIC)                              6.359e+00      9 45.422  < 2e-16 ***
    ## s(SIC):month_ordered12              9.738e-01      9 33.730  < 2e-16 ***
    ## s(SST_degC)                         4.698e+00      9 16.996 3.65e-05 ***
    ## s(SST_degC):month_ordered12         3.170e+00      9 10.505  0.00128 ** 
    ## s(lt_pack_ice)                      3.960e-01      9  0.670  0.18741    
    ## s(lt_pack_ice):month_ordered12      3.709e-05      9  0.000  0.49001    
    ## s(dist_ice_edge_km)                 3.528e+00      9 90.073  < 2e-16 ***
    ## s(dist_ice_edge_km):month_ordered12 3.933e+00      9 48.173  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.111   Deviance explained = 12.3%
    ## -REML = 1764.3  Scale est. = 1         n = 27972

This model outperforms the GAMs trained with environmental data obtained
from ACCESS-OM2-01. As we have done before, we will remove any variables
that do not contribute much to the model: `year`, bottom slope
(`bottom_slope_deg`) and long-term pack ice stability (`lt_pack_ice`),
which is non-significant. We will run the model once more to check if we
can improve it by simplifying it.

``` r
# Simplified model
simple_obs_model <- presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC) + s(SIC, by = month_ordered) + s(SST_degC) + s(SST_degC, by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, by = month_ordered)

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
    ## presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + 
    ##     s(SIC) + s(SIC, by = month_ordered) + s(SST_degC) + s(SST_degC, 
    ##     by = month_ordered) + s(dist_ice_edge_km) + s(dist_ice_edge_km, 
    ##     by = month_ordered)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -2.0849     0.3629  -5.744 9.23e-09 ***
    ## month12       1.4922     0.3677   4.058 4.95e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                        edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)                    0.7730      9  3.498 0.027875 *  
    ## s(dist_coast_km)                    0.9007      9  9.189 0.001143 ** 
    ## s(depth_m)                          2.3084      9 12.843 0.000552 ***
    ## s(SIC)                              6.3427      9 45.627  < 2e-16 ***
    ## s(SIC):month_ordered12              0.9755      9 33.946  < 2e-16 ***
    ## s(SST_degC)                         4.7249      9 17.136 3.45e-05 ***
    ## s(SST_degC):month_ordered12         3.1512      9 10.453 0.001314 ** 
    ## s(dist_ice_edge_km)                 3.5107      9 89.742  < 2e-16 ***
    ## s(dist_ice_edge_km):month_ordered12 3.9783      9 48.205  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.111   Deviance explained = 12.3%
    ## -REML = 1764.4  Scale est. = 1         n = 27972

#### Comparing models

We will compare these models using AIC and $r^{2}$. We will summarise
everything in the table below.

``` r
obs_sum_gam <- data.frame(model = c("full_model", "simple_obs_model"),
                      AIC = c(AIC(full_model_gam), AIC(simple_obs_model_gam)),
                      rsq = c(summary(full_model_gam)$r.sq, summary(simple_obs_model_gam)$r.sq))

#Rearraging data based on AIC
arrange(obs_sum_gam, AIC)
```

    ##              model      AIC       rsq
    ## 1 simple_obs_model 1722.663 0.1106735
    ## 2       full_model 1723.460 0.1109321

The simpler model has a lower AIC, and its $r^{2}$ value is quite
similar to that from the more complex model. Since we are aiming to keep
the most parsimonious model, we will keep the simplified version to
predict crabeater seal distributions.

#### Predictions

We will use the best performing model to predict crabeater seal
distribution using mean monthly environmental conditions obtained from
ACCESS-OM2-01.

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
saveRDS(mean_pred_match_obs_ras,
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
  lims(x = c(0, 5200000))+
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
  lims(x = c(0, 5200000))+
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
  draw_label("Mean crabeater seal distribution (Observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
plot_obs <- plot_grid(title, plot_obs, ncol = 1, 
                      rel_heights = c(0.1, 1))

#Saving graph
ggsave(file.path(out_folder, "map_mean_pred_obs.png"), 
       plot = plot_obs, device = "png", bg = "white")
```

    ## Saving 7 x 5 in image

### Differences

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
  lims(x = c(0, 5200000))+
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
  lims(x = c(0, 5200000))+
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

#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_obs_mod <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Differences in crabeater seal distribution (Full ACCESS-OM2-01\n vs Observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
plot_obs_mod <- plot_grid(title, plot_obs_mod, ncol = 1, 
                      rel_heights = c(0.075, 1))

#Saving graph
ggsave(file.path(out_folder, "map_comp_pred_obs_mod.png"), 
       plot = plot_obs_mod, device = "png", bg = "white")
```

    ## Saving 7 x 5 in image

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
  lims(x = c(0, 5200000))+
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
  lims(x = c(0, 5200000))+
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

#Remove legend from December plot
dec_plot <- dec_plot + theme(legend.position = 'none')

#Plotting together
plot_obs_mod_lim <- plot_grid(nov_plot, dec_plot, legend, ncol = 3, nrow = 1,
                            rel_widths = c(1, 1, 0.3))

#Add title
title <- ggdraw()+
  draw_label("Differences in crabeater seal distribution \n(Simplified ACCESS-OM2-01 vs Observations)",
             fontface = "bold", hjust = 0.5)+
  theme(plot.margin = margin(0, 0, 0, 0))

#Putting everything together
plot_obs_mod_lim <- plot_grid(title, plot_obs_mod_lim, ncol = 1, 
                      rel_heights = c(0.075, 1))

#Saving graph
ggsave(file.path(out_folder, "map_comp_pred_obs_mod_lim.png"), 
       plot = plot_obs_mod_lim, device = "png", bg = "white")
```

    ## Saving 7 x 5 in image
