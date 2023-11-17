Merging data
================
Denisse Fierro Arcos
2023-09-22

- <a href="#merging-background-and-presence-data"
  id="toc-merging-background-and-presence-data">Merging background and
  presence data</a>
  - <a href="#loading-libraries" id="toc-loading-libraries">Loading
    libraries</a>
  - <a href="#loading-data" id="toc-loading-data">Loading data</a>
    - <a href="#environmental-data-from-access-om2-01-model"
      id="toc-environmental-data-from-access-om2-01-model">Environmental data
      from ACCESS-OM2-01 model</a>
    - <a href="#environmental-data-from-observations"
      id="toc-environmental-data-from-observations">Environmental data from
      observations</a>
    - <a href="#matching-environmental-variables-in-model-with-observations"
      id="toc-matching-environmental-variables-in-model-with-observations">Matching
      environmental variables in model with observations</a>

# Merging background and presence data

Now that we have identified the best number of background points and had
performed an exploratory analysis of all datasets, we will merge
together background points and presence data that we will use in our
models.

## Loading libraries

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.3     ✔ purrr   1.0.2
    ## ✔ tibble  3.2.1     ✔ dplyr   1.1.2
    ## ✔ tidyr   1.3.0     ✔ stringr 1.5.0
    ## ✔ readr   2.1.3     ✔ forcats 0.5.2
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

## Loading data

### Environmental data from ACCESS-OM2-01 model

``` r
#Base folder containing data
base_folder <- "../../Environmental_Data"

#Loading crabeaters presence data with environmental data 
crabeaters <- file.path(base_folder, "ACCESS-OM2-01/unique_crabeater_obs_all_env.csv") %>% 
  read_csv() %>% 
  #Selecting observations for the Indian sector during the weaning period
  filter(str_detect(sector, "Indian") & life_stage == "weaning")
```

    ## Rows: 3240 Columns: 32
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr   (6): date, source, sector, zone, season_year, life_stage
    ## dbl  (25): latitude, longitude, year, yt_ocean, xt_ocean, month, decade, pre...
    ## dttm  (1): event_date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Loading background data
mod_bg_20 <- file.path(base_folder, "ACCESS-OM2-01/unique_background_20x_obs_all_env.csv") %>% 
  read_csv()
```

    ## Rows: 30671 Columns: 30
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (5): date, sector, zone, season_year, life_stage
    ## dbl (25): year, longitude, latitude, xt_ocean, yt_ocean, month, decade, pres...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Joining data together
crab_mod_all <- crabeaters %>% 
  #Ensure both datasets have the same columns
  select(all_of(names(mod_bg_20))) %>% 
  #Bind them together
  bind_rows(mod_bg_20)

#All environmental variables from ACCESS-OM2-01
#Define path out
path_out <- file.path(base_folder, "model_env_pres_bg_20x_Indian_weaning.csv")
#Save data
crab_mod_all %>% 
  #Keeping predictive variables
  select(c(year, month, xt_ocean, yt_ocean, presence:dist_ice_edge_km)) %>% 
  #Removing rows with NA values
  drop_na() %>% 
  #Saving to disk
  write_csv(path_out)
```

### Environmental data from observations

``` r
#Loading crabeaters presence data with environmental data 
crabeaters_obs <- file.path(base_folder, "Env_obs/unique_crabeater_obs_all_env.csv") %>% 
  read_csv() %>% 
  #Selecting observations for the Indian sector during the weaning period
  filter(str_detect(sector, "Indian") & life_stage == "weaning")
```

    ## Rows: 3240 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr   (6): date, source, sector, zone, season_year, life_stage
    ## dbl  (16): latitude, longitude, year, yt_ocean, xt_ocean, month, decade, pre...
    ## dttm  (1): event_date
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Loading background data
obs_bg_20 <- file.path(base_folder, "Env_obs/unique_background_20x_obs_all_env.csv") %>% 
  read_csv()
```

    ## Rows: 30671 Columns: 21
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (5): date, sector, zone, season_year, life_stage
    ## dbl (16): year, longitude, latitude, xt_ocean, yt_ocean, month, decade, pres...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Joining data together
crab_obs_all <- crabeaters_obs %>% 
  #Ensure both datasets have the same columns
  select(all_of(names(obs_bg_20))) %>% 
  #Bind them together
  bind_rows(obs_bg_20) %>% 
  #Ensuring areas north of sea ice edge are given a value of 0 for sea ice metrics
  mutate(SIC = case_when(is.na(SIC) & !is.na(SST_degC) ~ 0,
                         T ~ SIC),
         lt_pack_ice = case_when(is.na(lt_pack_ice) & !is.na(SST_degC) ~ 0,
                         T ~ lt_pack_ice))
  
  
#All environmental variables from observations
#Define path out
path_out <- file.path(base_folder, "obs_env_pres_bg_20x_Indian_weaning.csv")
#Save data
crab_obs_all %>% 
  select(c(year, month, xt_ocean, yt_ocean, presence:dist_ice_edge_km)) %>% 
  drop_na() %>% 
  write_csv(path_out)

#Getting names of variables available in observations
obs_var <- crab_obs_all %>% 
  select(c(year, month, xt_ocean, yt_ocean, presence:dist_ice_edge_km)) %>% 
  names()
```

### Matching environmental variables in model with observations

``` r
#Define path out
path_out <- file.path(base_folder, "mod-match-obs_env_pres_bg_20x_Indian_weaning.csv")
#Save data
crab_mod_all %>% 
  select(all_of(obs_var)) %>% 
  drop_na() %>% 
  write_csv(path_out)
```
