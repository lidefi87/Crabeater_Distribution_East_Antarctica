library(maxnet)
library(tidyverse)
library(stars)
library(starsExtra)
library(cmocean)

crab_ind <- read_csv("Environmental_Data/ACCESS-OM2-01/Obs_BG_5x_Indian_weaning_LowVIF.csv")%>% 
  #Removing any rows with NA values
  drop_na()

#Selecting variables (response and predictors) to be included in the model
model_data <- crab_ind %>% 
  select(presence:dist_ice_edge_km)

max_model <- maxnet(model_data$presence, model_data%>% select(!presence) %>% scale())

mean_model <- read_csv("Cleaned_Data/Env_obs/All_values_ACCESS-OM2-01_env_vars.csv") %>%
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
  labs(title = "Crabeater seal distribution during weaning season (Observations)",
       x = "Longitude",
       y = "Latitude",
       fill = "Probability")+
  coord_map("ortho", orientation = c(-90, 0, 0))
  # geom_point(data = crab_loc, aes(x = xt_ocean, y = yt_ocean, 
  #                                 color = as.factor(presence)), alpha = 0.3) 
  #scale_x_continuous(breaks = seq(40, 70, 10), limits = c(42, 70))+
  





