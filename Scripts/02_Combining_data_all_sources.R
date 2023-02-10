###############################################################################
# Joining crabeater seal data obtained from varios sources
# Author: Denisse Fierro Arcos
# Date: 2023-02-09
# 
# Crabeater seal data downloaded from various sources

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)

# Loading data ------------------------------------------------------------
emage_data <- read_csv("Cleaned_Data/EMAGE-II_seal_census_clean_data.csv") %>% 
  #Select columns of interest
  select(date_time, latitude, longitude, seals_number)

fil_data <- read_csv("Cleaned_Data/FIL_2014_Filchner_Outflow_Trough_seal_census_cleaned.csv") %>% 
  #Select columns of interest
  select(dt_trans_start_gmt, time_obs, lat_obs, lon_obs, number_ind) %>% 
  #Fix date time from time of transect start to time of observation
  mutate(dt_trans_start_gmt = as_datetime(paste(as_date(dt_trans_start_gmt),
                                                time_obs, " "))) %>% 
  #Remove time column
  select(-time_obs)

gbif_data <- read_csv("Cleaned_Data/GBIF_rgbif_cleaned.csv") %>% 
  #Select only observations from humans and material citation (data collected for a publication)
  filter(basisOfRecord != "MACHINE_OBSERVATION") %>% 
  #Select columns of interest
  select(eventDate, decimalLatitude, decimalLongitude, individualCount)

obis_data <- read_csv("Cleaned_Data/OBIS_cleaned.csv") %>% 
  #Select columns of interest
  select(eventdate, decimallatitude, decimallongitude, individualcount)

scar_data <- read_csv("Cleaned_Data/SCAR-APIS_cleaned.csv") %>% 
  #Select columns of interest
  select(eventDate, decimalLatitude, decimalLongitude, individualCount)

scar_bio <- read_csv("Cleaned_Data/SCAR-Biology_cleaned.csv") %>% 
  filter(basisOfRecord == "HUMAN_OBSERVATION") %>% 
  #Select columns of interest
  select(eventDate, decimalLatitude, decimalLongitude, individualCount)

#Renaming datasets
new_names <- c("date_time", "latitude", "longitude", "number_individuals")
names(emage_data) <- new_names
names(fil_data) <- new_names
names(gbif_data) <- new_names
names(obis_data) <- new_names
names(scar_data) <- new_names
names(scar_bio) <- new_names

#Bind them together
merged_data <- bind_rows(emage_data, fil_data, gbif_data, obis_data, scar_data, scar_bio) %>% 
  distinct(date_time, latitude, longitude, .keep_all = T) %>% 
  mutate(year = year(date_time), 
         month = month(date_time))
#Classifying data per decade
decades <- seq(1970, 2020, 10)
merged_data$decade<- decades[findInterval(merged_data$year, decades)]

#Removing datasets
rm(emage_data, fil_data, gbif_data, obis_data, scar_data, scar_bio)

#Plotting data
world <- rnaturalearth::ne_countries(returnclass = 'sf')
merged_data %>% 
  ggplot(aes(longitude, latitude))+
  geom_point(aes(colour = as.factor(decade)))+
  geom_sf(inherit.aes = F, data = world)+
  lims(y = c(-90, -45))+
  facet_wrap(~decade)

table(merged_data$decade)
table(merged_data$month)

  