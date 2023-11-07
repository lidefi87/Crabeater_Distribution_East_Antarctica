# Loading libraries
library(terra)
library(sf)
library(tidyverse)
library(measoshapes)
source("Scripts/01_Bio_data/04a_randomPoints.R")

#Loading MEASO regions and assign a more current version of CRS (as per sf warning)
indian <- measo_regions05_ll
#Maintaining original WGS84
st_crs(indian) <- 4326

#Getting extent for Central Indian (CI) and East Indian (EI) MEASO sectors to
#create mask
indian_ext <- indian %>%
  filter((str_starts(name, "CI") | str_starts(name, "EI")) & !str_ends(name, "T")) %>% 
  st_bbox() %>% 
  ext()

#Creating ocean grid using Indian MEASO sectors
grid_ind <- rast(indian_ext, resolution = 0.1, vals = 1, crs = "WGS84")

#Loading original crabeater observations
ind_crab <- read_csv("Biological_Data/unique_crabeater_obs_grid.csv") %>% 
  filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = F)

#Creating background points
#November
bg_pts_nov <- ind_crab %>% 
  filter(month == 11) %>% 
  bg_pts(grid_ind)

#December
bg_pts_dec <- ind_crab %>% 
  filter(month == 12) %>% 
  bg_pts(grid_ind)

#Stacking background points
bg_pts_all <- bind_rows(bg_pts_nov, bg_pts_dec)

#Checking results
ggplot()+
  geom_point(data = bg_pts_all, aes(x = longitude, y = latitude),
             color = "orange", alpha = 0.5)+
  geom_point(data = ind_crab, aes(xt_ocean, yt_ocean))+
  facet_wrap(~month)

bg_pts_all %>% 
  write_csv("Biological_Data/BG_points/Background_20xPoints_Indian-Sectors_weaning.csv")
