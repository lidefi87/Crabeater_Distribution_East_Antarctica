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
  mutate(month = factor(month)) %>% 
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

    ## Rows: 1147094 Columns: 11
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

#### Building GAM formula

We will run an initial GAM using all available covariates. Sea ice
related variables will be fitted individually to each `month` as we
expect sea ice to vary in extent between November and December because
sea ice should be retreating during this time.

``` r
# Most complex model
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)
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
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept) -0.41338    0.19539  -2.116   0.0344 *
    ## month12      0.02187    0.21223   0.103   0.9179  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(year)                     0.1943607      9  0.240 0.258690    
    ## s(bottom_slope_deg)         0.0003329      9  0.000 0.661442    
    ## s(dist_shelf_km)            0.0003342      9  0.000 0.551597    
    ## s(dist_coast_km)            0.0005639      9  0.000 0.407636    
    ## s(depth_m)                  1.5255669      9  8.605 0.002028 ** 
    ## s(SIC):month11              1.5407435      9  3.882 0.045425 *  
    ## s(SIC):month12              4.7777467      9 47.089  < 2e-16 ***
    ## s(lt_pack_ice):month11      2.5367298      9 15.352 5.45e-05 ***
    ## s(lt_pack_ice):month12      0.0003215      9  0.000 0.536723    
    ## s(dist_ice_edge_km):month11 1.8075696      9  9.978 0.000829 ***
    ## s(dist_ice_edge_km):month12 3.4099997      9 26.579 3.14e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0324   Deviance explained = 2.94%
    ## -REML = 1870.5  Scale est. = 1         n = 24276

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
    ## (Intercept) -0.40639    0.20149  -2.017   0.0437 *
    ## month12      0.01464    0.21758   0.067   0.9464  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  1.5139002      9  8.504 0.002286 ** 
    ## s(SIC):month11              1.7580525      9  5.358 0.025011 *  
    ## s(SIC):month12              4.7510137      9 47.051  < 2e-16 ***
    ## s(lt_pack_ice):month11      2.4274214      9 15.186 8.00e-05 ***
    ## s(lt_pack_ice):month12      0.0001512      9  0.000 0.494854    
    ## s(dist_ice_edge_km):month11 0.9083244      9  9.537 0.000695 ***
    ## s(dist_ice_edge_km):month12 3.4151266      9 26.580 3.33e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0316   Deviance explained =  2.9%
    ## -REML = 1870.3  Scale est. = 1         n = 24276

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
    ## 1       full_model 1872.768 0.03244056
    ## 2 significant_only 1870.465 0.03163787

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
    ## (Intercept)  -0.9796     0.4798  -2.042   0.0412 *
    ## month12       0.5474     0.4883   1.121   0.2623  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  0.8788914      9  6.903 0.002426 ** 
    ## s(SIC):month11              0.0001622      9  0.000 0.710360    
    ## s(SIC):month12              5.5737991      9 36.482  < 2e-16 ***
    ## s(lt_pack_ice):month11      2.7168525      9 17.686 2.26e-05 ***
    ## s(lt_pack_ice):month12      1.1767952      9  2.987 0.055273 .  
    ## s(dist_ice_edge_km):month11 1.9161047      9 11.526 0.000426 ***
    ## s(dist_ice_edge_km):month12 2.5014750      9 14.384 0.000211 ***
    ## s(SST_degC):month11         2.3221183      9  6.765 0.024393 *  
    ## s(SST_degC):month12         5.2331425      9 43.638  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0488   Deviance explained = 4.76%
    ## -REML = 1850.4  Scale est. = 1         n = 24276

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
    ## (Intercept)  -0.9810     0.4764  -2.059   0.0395 *
    ## month12       0.5337     0.4795   1.113   0.2656  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq  p-value    
    ## s(depth_m)                  0.9109974      9 10.016 0.000582 ***
    ## s(lt_pack_ice):month11      2.7116209      9 17.786 2.11e-05 ***
    ## s(lt_pack_ice):month12      0.0002319      9  0.000 0.943751    
    ## s(dist_ice_edge_km):month11 1.9257020      9 11.650 0.000407 ***
    ## s(dist_ice_edge_km):month12 3.4201110      9 24.091 1.09e-05 ***
    ## s(SST_degC):month11         2.3063267      9  6.649 0.025742 *  
    ## s(SST_degC):month12         4.9434692      9 50.288  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0361   Deviance explained =  3.6%
    ## -REML = 1860.1  Scale est. = 1         n = 24276

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
    ## 1    significant_only_SST 1844.186 0.04878829
    ## 2 significant_only_no_SIC 1853.806 0.03613065
    ## 3        significant_only 1870.465 0.03163787
    ## 4              full_model 1872.768 0.03244056

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
    ## 1     24241     3604.5                               
    ## 2     24250     3648.1 -8.9112  -43.687 1.491e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

We can see these models differ significantly, but since
multicollinearity has a negative effect on GAM predicitive ability, we
will use the `significant_only_no_SIC` to predict crabeater
distribution.

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
  mutate(month = factor(month))
```

    ## Rows: 32368 Columns: 22
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
  append(c("SST_degC"))
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
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(vel_lon_bottom_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)
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
    ##     s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, 
    ##     by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, 
    ##     by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, 
    ##     by = month) + s(vel_lon_bottom_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)  -0.4442     0.2519  -1.764   0.0778 .
    ## month12      -0.1052     0.2594  -0.406   0.6849  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(year)                        0.0003151      9  0.000 0.517675    
    ## s(bottom_slope_deg)            0.0002029      9  0.000 0.611350    
    ## s(dist_shelf_km)               0.8930031      9  8.642 0.000609 ***
    ## s(dist_coast_km)               2.1288806      9  8.376 0.006257 ** 
    ## s(depth_m)                     0.8281525      9  4.855 0.009059 ** 
    ## s(SIC):month11                 1.8447611      9  6.314 0.013363 *  
    ## s(SIC):month12                 2.1815213      9  6.663 0.022670 *  
    ## s(bottom_temp_degC):month11    2.6643447      9 15.852 0.000121 ***
    ## s(bottom_temp_degC):month12    0.0002193      9  0.000 0.869086    
    ## s(SSS_psu):month11             1.0666118      9  2.665 0.053239 .  
    ## s(SSS_psu):month12             4.1752225      9 24.974 8.13e-06 ***
    ## s(vel_lat_surf_msec):month11   4.4742517      9 21.695 6.76e-05 ***
    ## s(vel_lat_surf_msec):month12   0.8319600      9  4.940 0.014208 *  
    ## s(vel_lat_bottom_msec):month11 0.8717355      9  6.757 0.004854 ** 
    ## s(vel_lat_bottom_msec):month12 0.0002119      9  0.000 0.970733    
    ## s(vel_lon_surf_msec):month11   1.7854495      9 10.721 0.000720 ***
    ## s(vel_lon_surf_msec):month12   4.2652443      9 17.489 0.000631 ***
    ## s(vel_lon_bottom_msec):month11 0.0003354      9  0.000 0.602947    
    ## s(vel_lon_bottom_msec):month12 0.0003596      9  0.000 0.460549    
    ## s(lt_pack_ice):month11         2.8157975      9 21.466 1.56e-06 ***
    ## s(lt_pack_ice):month12         0.8492468      9  3.665 0.028235 *  
    ## s(dist_ice_edge_km):month11    0.9407129      9 15.086 1.91e-05 ***
    ## s(dist_ice_edge_km):month12    3.0392275      9 16.579 0.000152 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0678   Deviance explained =  6.3%
    ## -REML = 1844.1  Scale est. = 1         n = 24276

This model explains more of the variability in our data than the model
that uses variables for which we also have observations. We will remove
the variables with a low contribution (`edf` $\sim 0.001$). This
includes `year`, bottom slope (`bottom_slope_deg`), and bottom
meridional water velocity (`vel_lon_bottom_msec`).

We will defined this simplified model and test it below.

``` r
# Simplified model
simpler_model <- presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

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
    ## presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + 
    ##     s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, 
    ##     by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, 
    ##     by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)  -0.4446     0.2524  -1.762   0.0781 .
    ## month12      -0.1048     0.2599  -0.403   0.6869  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)               8.928e-01      9  8.633 0.000612 ***
    ## s(dist_coast_km)               2.125e+00      9  8.330 0.006404 ** 
    ## s(depth_m)                     8.279e-01      9  4.850 0.009090 ** 
    ## s(SIC):month11                 1.854e+00      9  6.398 0.012803 *  
    ## s(SIC):month12                 2.183e+00      9  6.673 0.022578 *  
    ## s(bottom_temp_degC):month11    2.668e+00      9 15.892 0.000119 ***
    ## s(bottom_temp_degC):month12    1.500e-05      9  0.000 0.867558    
    ## s(SSS_psu):month11             1.076e+00      9  2.703 0.052233 .  
    ## s(SSS_psu):month12             4.175e+00      9 24.961 8.17e-06 ***
    ## s(vel_lat_surf_msec):month11   4.474e+00      9 21.675 6.86e-05 ***
    ## s(vel_lat_surf_msec):month12   8.319e-01      9  4.940 0.014215 *  
    ## s(vel_lat_bottom_msec):month11 8.717e-01      9  6.760 0.004851 ** 
    ## s(vel_lat_bottom_msec):month12 1.754e-05      9  0.000 0.976517    
    ## s(vel_lon_surf_msec):month11   1.717e+00      9 10.623 0.000716 ***
    ## s(vel_lon_surf_msec):month12   4.265e+00      9 17.491 0.000630 ***
    ## s(lt_pack_ice):month11         2.815e+00      9 21.447 1.62e-06 ***
    ## s(lt_pack_ice):month12         8.515e-01      9  3.672 0.028184 *  
    ## s(dist_ice_edge_km):month11    9.397e-01      9 15.091 1.90e-05 ***
    ## s(dist_ice_edge_km):month12    3.039e+00      9 16.578 0.000152 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0678   Deviance explained = 6.29%
    ## -REML = 1844.1  Scale est. = 1         n = 24276

We can see that this simplified model performs almost the same as the
model including all environmental variables available.

We will now explore the influence of `SST` as it was identified as an
influential variable in the previous dataset.

``` r
# Simplified model with SST
simpler_model_SST <- presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, by = month)

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
    ##     s(SIC, by = month) + s(bottom_temp_degC, by = month) + s(SSS_psu, 
    ##     by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, 
    ##     by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -0.6039     0.2276  -2.653  0.00797 **
    ## month12       0.1291     0.2461   0.525  0.59972   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)               0.4036910      9  0.752 0.138809    
    ## s(dist_coast_km)               0.9322265      9  1.960 0.102063    
    ## s(depth_m)                     0.8838414      9  7.549 0.000797 ***
    ## s(SIC):month11                 1.8807990      9  6.602 0.011828 *  
    ## s(SIC):month12                 5.7025217      9 39.603  < 2e-16 ***
    ## s(bottom_temp_degC):month11    2.6664893      9 17.162 9.88e-05 ***
    ## s(bottom_temp_degC):month12    0.0002078      9  0.000 0.662363    
    ## s(SSS_psu):month11             0.0526950      9  0.054 0.279930    
    ## s(SSS_psu):month12             0.0008130      9  0.000 0.566080    
    ## s(vel_lat_surf_msec):month11   4.5178675      9 23.485 2.64e-05 ***
    ## s(vel_lat_surf_msec):month12   1.2578481      9  4.484 0.024441 *  
    ## s(vel_lat_bottom_msec):month11 0.8654703      9  6.409 0.006038 ** 
    ## s(vel_lat_bottom_msec):month12 0.0001628      9  0.000 0.877933    
    ## s(vel_lon_surf_msec):month11   1.2295749      9  9.673 0.000865 ***
    ## s(vel_lon_surf_msec):month12   3.6620571      9 15.557 0.000800 ***
    ## s(lt_pack_ice):month11         2.8140657      9 19.772 7.25e-06 ***
    ## s(lt_pack_ice):month12         0.0009246      9  0.001 0.377936    
    ## s(dist_ice_edge_km):month11    0.9320740      9 13.086 7.15e-05 ***
    ## s(dist_ice_edge_km):month12    2.5976511      9 16.901 0.000114 ***
    ## s(SST_degC):month11            0.0002582      9  0.000 0.511333    
    ## s(SST_degC):month12            5.2930472      9 47.248  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0744   Deviance explained = 7.25%
    ## -REML = 1829.1  Scale est. = 1         n = 24276

This model is slightly worse than when `SST` was not included. This is
likely because `SST` was highly correlated to `SIC`. Therefore, we will
not use `SST` to predict crabeater seal distribution.

Now that `SST` has been included, a few variables have become
non-significant (e.g., distance to the continental shelf,
`dist_shelf_km`, and distance to the coast, `dist_coast_km`). We can
also see that its predictive ability increased. Since we know `SIC` is
hight correlated to `SST`, we will check the effect that removing `SIC`
has on predictive ability.

``` r
# Simplified model with SST and no SIC
simpler_model_SST_noSIC <- presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, by = month)

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
    ## presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + 
    ##     s(bottom_temp_degC, by = month) + s(SSS_psu, by = month) + 
    ##     s(vel_lat_surf_msec, by = month) + s(vel_lat_bottom_msec, 
    ##     by = month) + s(vel_lon_surf_msec, by = month) + s(lt_pack_ice, 
    ##     by = month) + s(dist_ice_edge_km, by = month) + s(SST_degC, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)  -1.0613     0.5218  -2.034   0.0419 *
    ## month12       0.5220     0.5257   0.993   0.3207  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                      edf Ref.df Chi.sq  p-value    
    ## s(dist_shelf_km)               8.430e-01      9  5.587 0.003930 ** 
    ## s(dist_coast_km)               1.889e+00      9  6.816 0.011457 *  
    ## s(depth_m)                     7.740e-01      9  3.456 0.021055 *  
    ## s(bottom_temp_degC):month11    2.585e+00      9 15.091 0.000193 ***
    ## s(bottom_temp_degC):month12    1.082e-04      9  0.000 0.740281    
    ## s(SSS_psu):month11             6.356e-01      9  1.097 0.138190    
    ## s(SSS_psu):month12             3.615e+00      9 21.418 2.15e-05 ***
    ## s(vel_lat_surf_msec):month11   4.455e+00      9 21.626 6.81e-05 ***
    ## s(vel_lat_surf_msec):month12   8.286e-01      9  4.781 0.015527 *  
    ## s(vel_lat_bottom_msec):month11 8.647e-01      9  6.375 0.006087 ** 
    ## s(vel_lat_bottom_msec):month12 4.508e-05      9  0.000 0.907245    
    ## s(vel_lon_surf_msec):month11   1.746e+00      9 10.467 0.000802 ***
    ## s(vel_lon_surf_msec):month12   4.339e+00      9 15.878 0.001445 ** 
    ## s(lt_pack_ice):month11         2.828e+00      9 21.296 5.91e-06 ***
    ## s(lt_pack_ice):month12         3.504e-01      9  0.536 0.199137    
    ## s(dist_ice_edge_km):month11    9.388e-01      9 14.639 2.67e-05 ***
    ## s(dist_ice_edge_km):month12    2.768e+00      9 16.056 0.000137 ***
    ## s(SST_degC):month11            2.397e+00      9  7.316 0.019281 *  
    ## s(SST_degC):month12            4.753e+00      9 29.023  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0716   Deviance explained = 6.93%
    ## -REML = 1835.4  Scale est. = 1         n = 24276

#### Comparing models

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
    ## 1       simpler_model_SST 1820.146 0.07437048
    ## 2 simpler_model_SST_noSIC 1828.149 0.07163717
    ## 3           simpler_model 1842.591 0.06778572
    ## 4              full_model 1842.686 0.06781712

Based on AIC and $r^{2}$, we can see that the models including `SST` are
the best performing. The exclusion of `SIC` results in a slightly lower
predictive ability, but we will use this latter model for predictions as
we know GAMs are affected by multicollinearity.

#### Predictions

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

### Environmental variables from observations

Finally, we will create a model using all environmental data obtained
from observations.

``` r
#Loading data
obs_env_data <- read_csv(str_subset(file_list, "/obs")) %>% 
  #Setting month as factor and ordered factor
  mutate(month = factor(month))
```

    ## Rows: 32315 Columns: 13
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
full_model <- presence ~ month + s(year) + s(bottom_slope_deg) + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(SST_degC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)
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
    ##     s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(SST_degC, 
    ##     by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.1049     0.3404  -3.246  0.00117 **
    ## month12       0.6185     0.3452   1.792  0.07318 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq p-value    
    ## s(year)                     0.3606070      9  0.520 0.21600    
    ## s(bottom_slope_deg)         0.0001093      9  0.000 0.84793    
    ## s(dist_shelf_km)            0.0001364      9  0.000 0.93537    
    ## s(dist_coast_km)            0.1305543      9  0.149 0.28364    
    ## s(depth_m)                  0.7228757      9  2.549 0.05405 .  
    ## s(SIC):month11              0.3723062      9  0.582 0.20669    
    ## s(SIC):month12              5.9197207      9 77.286 < 2e-16 ***
    ## s(SST_degC):month11         2.3732842      9  7.974 0.01033 *  
    ## s(SST_degC):month12         4.5287493      9 14.901 0.00210 ** 
    ## s(lt_pack_ice):month11      0.0001880      9  0.000 0.62617    
    ## s(lt_pack_ice):month12      2.0693935      9  9.266 0.00433 ** 
    ## s(dist_ice_edge_km):month11 0.0003432      9  0.000 0.60804    
    ## s(dist_ice_edge_km):month12 0.7596504      9  3.081 0.03827 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.054   Deviance explained = 5.85%
    ## -REML = 1816.4  Scale est. = 1         n = 24236

This model outperforms the GAMs trained with environmental data obtained
from ACCESS-OM2-01. As we have done before, we will remove any variables
that do not contribute much to the model: bottom slope
(`bottom_slope_deg`) and distance to the continental shelf
(`dist_shelf_km`).

``` r
# Simplified model
simple_obs_model <- presence ~ month + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(SST_degC, by = month)+ s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

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
    ## presence ~ month + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + 
    ##     s(SST_degC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -1.1165     0.3387  -3.296 0.000981 ***
    ## month12       0.6330     0.3433   1.844 0.065249 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq p-value    
    ## s(dist_coast_km)            0.1106666      9  0.125 0.28729    
    ## s(depth_m)                  0.7270594      9  2.597 0.05256 .  
    ## s(SIC):month11              0.3395906      9  0.507 0.21848    
    ## s(SIC):month12              5.9169414      9 76.996 < 2e-16 ***
    ## s(SST_degC):month11         2.3692546      9  8.016 0.01027 *  
    ## s(SST_degC):month12         4.5043515      9 14.717 0.00229 ** 
    ## s(lt_pack_ice):month11      0.0001673      9  0.000 0.60940    
    ## s(lt_pack_ice):month12      2.1160232      9  9.376 0.00448 ** 
    ## s(dist_ice_edge_km):month11 0.0001487      9  0.000 0.63027    
    ## s(dist_ice_edge_km):month12 0.7563236      9  3.049 0.03908 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0536   Deviance explained = 5.83%
    ## -REML = 1816.4  Scale est. = 1         n = 24236

This simplified model explains almost the same variance as the full
model.

Although multicollinearity was not identified in the observational
dataset, we will assess the effect of keeping either `SST` or `SIC` on
the predictive ability of the GAM.

``` r
# Simplified model no SST
simple_obs_model_noSST <- presence ~ month + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)

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
    ## presence ~ month + s(dist_coast_km) + s(depth_m) + s(SIC, by = month) + 
    ##     s(lt_pack_ice, by = month) + s(dist_ice_edge_km, by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.57130    0.07975  -7.164 7.86e-13 ***
    ## month12      0.11922    0.09443   1.263    0.207    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df  Chi.sq p-value    
    ## s(dist_coast_km)            0.7984065      9   1.384 0.15798    
    ## s(depth_m)                  0.7967682      9   3.733 0.02639 *  
    ## s(SIC):month11              1.2563555      9   4.331 0.02925 *  
    ## s(SIC):month12              6.6374081      9 108.541 < 2e-16 ***
    ## s(lt_pack_ice):month11      0.0006402      9   0.000 0.77816    
    ## s(lt_pack_ice):month12      2.0098060      9   7.956 0.00865 ** 
    ## s(dist_ice_edge_km):month11 0.0011088      9   0.000 0.53481    
    ## s(dist_ice_edge_km):month12 1.0552114      9   4.036 0.02456 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0467   Deviance explained = 5.21%
    ## -REML = 1819.6  Scale est. = 1         n = 24236

Removing `SST` from the model results in a decline in the model
predictive ability of about 0.5%. We will now check results if `SIC` is
removed.

``` r
# Simplified model no SIC
simple_obs_model_noSIC <- presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + s(SST_degC, by = month) + s(SST_degC, by = month) + s(dist_ice_edge_km, by = month)

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
    ## presence ~ month + s(dist_shelf_km) + s(dist_coast_km) + s(depth_m) + 
    ##     s(SST_degC, by = month) + s(SST_degC, by = month) + s(dist_ice_edge_km, 
    ##     by = month)
    ## 
    ## Parametric coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.1032     0.3368  -3.276  0.00105 **
    ## month12       0.6469     0.3484   1.857  0.06336 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                                   edf Ref.df Chi.sq p-value    
    ## s(dist_shelf_km)            0.0009432      9  0.001 0.41694    
    ## s(dist_coast_km)            0.3713814      9  0.601 0.19840    
    ## s(depth_m)                  1.3642176      9  6.237 0.00885 ** 
    ## s(SST_degC):month11         2.3499126      9  8.457 0.01096 *  
    ## s(SST_degC):month12         5.9977528      9 38.345 < 2e-16 ***
    ## s(dist_ice_edge_km):month11 0.0017107      9  0.001 0.57759    
    ## s(dist_ice_edge_km):month12 3.6430688      9 11.279 0.01017 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.0319   Deviance explained = 3.39%
    ## -REML = 1858.1  Scale est. = 1         n = 24236

The effect of removing `SIC` is much larger on the model predictive
ability. Finally, we will apply the same model fitted to the ACCESS
data.

We will not apply the same model fitted to the ACCESS data because we
did not use `SIC`.

#### Comparing models

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
    ## 1       simple_obs_model 1801.604 0.05361745
    ## 2             full_model 1802.363 0.05402606
    ## 3 simple_obs_model_noSST 1806.082 0.04667752
    ## 4 simple_obs_model_noSIC 1846.489 0.03186003

The simple model has the highest AIC, and its $r^{2}$, so we will use
this to predict crabeater distribution.

#### Predictions

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

### Differences

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
