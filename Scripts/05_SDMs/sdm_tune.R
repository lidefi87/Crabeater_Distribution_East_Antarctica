library(SDMtune)
library(tidyverse)
test <- read.csv("Environmental_Data/ACCESS-OM2-01/Obs_BG_5x_Indian_weaning_LowVIF_match-obs.csv") %>% 
  drop_na()

pres_coord <- test %>% 
  filter(presence == 1) %>% 
  select(xt_ocean, yt_ocean) %>% 
  rename("x" = "xt_ocean", "y" = "yt_ocean")
back_coord <- test %>% 
  filter(presence == 0) %>% 
  select(xt_ocean, yt_ocean) %>% 
  rename("x" = "xt_ocean", "y" = "yt_ocean")
pred <- test %>% 
  select(bottom_slope_deg:dist_ice_edge_km)
  # select(year, month, decade, bottom_slope_deg:dist_ice_edge_km) %>% 
  # mutate_at(vars("year", "month", "decade"), as.factor)

data_test <- SWD(species = "crabeaters",
    coords = rbind(pres_coord, back_coord),
    data = pred,
    pa = test$presence)

datasets <- trainValTest(data_test, test = 0.25, only_presence = T)

default_model <- train(method = "Maxent", data = datasets[[1]])

slotNames(default_model@model)

default_model <- train(method = "Maxent", data = data_test,
                       fc = "lh", reg = 0.5, iter = 500)

pred_1 <- predict(default_model, data = data_test, type = "cloglog")

all_vals <- read_csv("Environmental_Data/ACCESS-OM2-01/All_values_ACCESS-OM2-01_env_vars.csv")

pred_2 <- predict(default_model, data = all_vals, type = "cloglog")

predictions <- all_vals %>% 
  drop_na() %>% 
  select(xt_ocean, yt_ocean) %>% 
  mutate(pred = as.vector(pred_2))

library(stars)
library(cmocean)
ras_pred <- predictions %>%
  st_as_stars() %>% 
  st_set_crs(4326)

ras_pred_trans <- ras_pred %>% 
  st_transform(crs = 3976)

ggplot() +
  geom_stars(data = ras_pred_trans) +
  scale_fill_cmocean(name = "ice", direction = -1, 
                     guide = guide_colorbar(barwidth = 1, barheight = 10, 
                                            ticks = FALSE, nbin = 1000, 
                                            frame.colour = "black"), 
                     limits = c(0, 1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Crabeater seal distribution during weaning season (ACCESS-OM2-01)",
       x = "Longitude",
       y = "Latitude",
       fill = "Probability")

ths <- thresholds(default_model, 
                  type = "cloglog")
ths

auc(default_model)

plotROC(default_model, test =  datasets[[2]])


