library(raster)
library(sf)
library(tidyverse)
library(spatialEco)
library(dismo)

#Create sample raster
mod_grid <- raster(ymn = -80, ymx = -45, xmn = 38, xmx = 160, resolution = 0.1, 
                   vals = 1, crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")

ind_crab <- crabeaters %>% 
  filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = F)

yext <- (max(ind_crab$latitude) - min(ind_crab$latitude))/5
xext <- (max(ind_crab$longitude) - min(ind_crab$longitude))/5
max_ext <- max(c(yext, xext))

tgb_kde <- sp.kde(x = ind_crab, bw = round(max_ext, 2),
                  ref = terra::rast(mod_grid), 
                  standardize = TRUE, scale.factor = 10000)

# generate 50k random samples from the KDE raster file
kde_samples <- randomPoints(raster(tgb_kde), n = 10000, prob = TRUE)
kde_samples <- kde_samples %>% 
  as.data.frame() %>% 
  rename("longitude" = "x", "latitude" = "y") 

set.seed(42)
back_samples <- kde_samples %>% 
  slice_sample(n = nrow(ind_crab), replace = F)
  
back_samples <- ind_crab %>% 
  st_drop_geometry() %>% 
  dplyr::select(event_date, date, year, month, season_yea, life_stage, decade) %>% 
  cbind(back_samples) %>% 
  mutate(number_ind = 0, basis_reco = "BACKGROUND")

back_samples


ind_crab %>% 
  ggplot()+
  geom_sf(color = "red")+
  geom_point(inherit.aes = F, data = back_samples, aes(x = longitude, y = latitude),
             color = "green", alpha = 0.5)+
  facet_grid(basis_reco~.)

back_samples %>% 
  write_csv("Cleaned_Data/Background_Points_Indian-Sectors_weaning.csv")
