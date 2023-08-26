library(maxnet)
library(tidyverse)
library(stars)
library(starsExtra)
library(cmocean)

crab_ind <- read_csv("Environmental_Data/ACCESS-OM2-01/Obs_BG_5x_Indian_weaning_LowVIF_match-obs.csv")%>% 
  #Removing any rows with NA values
  drop_na()

#Selecting variables (response and predictors) to be included in the model
model_data <- crab_ind %>% 
  select(bottom_slope_deg:dist_ice_edge_km) %>% 
  scale() %>% 
  as.data.frame() %>% 
  bind_cols(crab_ind %>% select(presence))

# model_data <- crab_ind %>% 
#   select(presence:dist_ice_edge_km, month) %>% 
#   mutate(month = factor(month))

max_model <- maxnet(model_data$presence, 
                    model_data %>% select(!presence))

mean_model <- read_csv("Environmental_Data/ACCESS-OM2-01/All_values_ACCESS-OM2-01_env_vars.csv")

mean_model <- mean_model %>% 
  select(!c(xt_ocean, yt_ocean)) %>%
  scale() %>% 
  as.data.frame() %>% 
  cbind(mean_model %>% select(c(xt_ocean, yt_ocean))) %>% 
  drop_na()

predicted <- predict(max_model, mean_model, 
                     clamp = F, type = "cloglog")

pred <- mean_model %>% 
  select(xt_ocean, yt_ocean) %>% 
  mutate(pred = as.vector(predicted))

ras_pred <- pred %>%
  st_as_stars()


crab_loc <- crab_ind %>%
  select(xt_ocean, yt_ocean, presence)

ggplot() +
  geom_stars(data = ras_pred) +
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
       fill = "Probability")#+
  # geom_point(data = crab_loc, aes(x = xt_ocean, y = yt_ocean, 
  #                                 color = as.factor(presence)), alpha = 0.3) 
  #scale_x_continuous(breaks = seq(40, 70, 10), limits = c(42, 70))+
  





