# Loading libraries
library(raster)
library(sf)
library(stars)
library(tidyverse)
library(spatialEco)
library(dismo)
library(measoshapes)

#Getting the Central Indian (CI) and East Indian (EI) MEASO sectors to apply as mask
indian <- measo_regions05_ll %>%
  filter((str_starts(name, "CI") | str_starts(name, "EI")) & !str_ends(name, "T")) %>% 
  st_union()

#Creating ocean grid
grid_ocean <- raster(xmn = 0, xmx = 180, ymn = -80, ymx = -45, resolution = 0.1,
               vals = 1, crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")

#Subsetting grid to Indian sectors
grid_ind <- crop(grid_ocean, extent(st_bbox(indian)))
#Checking results
plot(grid_ind)

#Loading original crabeater observations
ind_crab <- read_csv("Cleaned_Data/All_sources_clean_data_MEASO.csv") %>% 
  filter(basis_record == "HUMAN_OBSERVATION") %>% 
  filter(str_detect(sector, "Indian") & life_stage == "weaning") %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = F)

#Getting maximum extent of crabeater observations to calculate kernel for KDE
yext <- (max(ind_crab$latitude) - min(ind_crab$latitude))/5
xext <- (max(ind_crab$longitude) - min(ind_crab$longitude))/5
max_ext <- max(c(yext, xext))

#Calculating KDE
tgb_kde <- sp.kde(x = ind_crab, bw = round(max_ext, 2),
                  ref = terra::rast(grid_ind), 
                  standardize = TRUE, scale.factor = 10000)
#Checking results
plot(tgb_kde)


#Generating random samples from the KDE raster file
set.seed(42)

#Getting total number of samples from gridded data
uniq_obs_grid <- read_csv("Cleaned_Data/unique_crabeater_obs_all_env.csv") %>% 
  filter(str_detect(sector, "Indian") & life_stage == "weaning")

#Getting 20 times the number of unique observations
source("Scripts/randomPoints.R")
kde_samples <- randomPoints2(mask = raster(tgb_kde), n = nrow(uniq_obs_grid)*20, prob = T, 
                             tryf = 2, replace = T)

#Extracting background data points
back_samples <- kde_samples %>% 
  as.data.frame() %>% 
  rename("longitude" = "x", "latitude" = "y") 

#Matching background points to observation dates  
#Creating a vector with observation indices selected at random 20 times to match background
ind <- vector()
#Repeating process 20 times because we have 20 times as many backgrounds than observations
for (i in 1:20) {
  rand <- sample(1:nrow(uniq_obs_grid), size = nrow(uniq_obs_grid), replace = F)
  ind <- append(ind, rand)
}

#Applying indices to data frame with unique obs
back_samples <- uniq_obs_grid[ind,] %>% 
  select(date:month, season_year:presence) %>% 
  mutate(presence  = 0) %>% 
  cbind(back_samples)

uniq_obs_grid %>% 
  ggplot(aes(xt_ocean, yt_ocean))+
  geom_point(color = "red")+
  geom_point(inherit.aes = F, data = back_samples, aes(x = longitude, y = latitude),
             color = "green", alpha = 0.5)

back_samples %>% 
  write_csv("Cleaned_Data/Background_20xPoints_Indian-Sectors_weaning.csv")
