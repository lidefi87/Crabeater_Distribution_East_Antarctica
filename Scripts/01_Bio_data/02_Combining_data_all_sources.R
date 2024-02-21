###############################################################################
# Joining crabeater seal data obtained from varios sources
# Author: Denisse Fierro Arcos
# Date: 2023-02-09
# 
# Crabeater seal data downloaded from various sources

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)

# Setting new column names ------------------------------------------------
new_names <- c("event_date", "latitude", "longitude", "number_individuals", "basis_record", "source")

# Loading data ------------------------------------------------------------
#ANT/EMAGE censuses
emage_data <- read_csv("Biological_Data/Cleaned_Data/EMAGE-II_seal_census_clean_data.csv") %>% 
  #Select columns of interest
  select(date_time, latitude, longitude, seals_number) %>% 
  mutate(basisOfRecord = "HUMAN_OBSERVATION",
         source = "ANT-EMAGE")
names(emage_data) <- new_names

#Belgica expedition
belgica <- read_csv("Biological_Data/Cleaned_Data/Belgica121_cleaned.csv") %>% 
  select(eventDate, decimalLatitude, decimalLongitude, individualCount, basisOfRecord) %>% 
  mutate(source = "Belgica_121")
names(belgica) <- new_names

#Bornemann ARGOS tracking data for 13 seals
bornemann <- read_csv("Biological_Data/Cleaned_Data/Bornemann_ARGOS.csv") %>% 
  select(date_time, latitude, longitude, l_carcinophaga_number, basisOfRecord) %>% 
  mutate(source = "Bornemann_Argos")
names(bornemann) <- new_names

#SCAR seabirds
scar_seabirds <- read_csv("Biological_Data/Cleaned_Data/Cleaned_ASAC_2208_seabirds.csv") %>% 
  select(observation_date, latitude, longitude, species_count, basisOfRecord) %>% 
  mutate(source = "SCAR_seabirds")
names(scar_seabirds) <- new_names

#Filchner data
fil_data <- read_csv("Biological_Data/Cleaned_Data/FIL_2014_Filchner_Outflow_Trough_seal_census_cleaned.csv") %>% 
  #Select columns of interest
  select(dt_trans_start_gmt, time_obs, lat_obs, lon_obs, number_ind) %>% 
  #Fix date time from time of transect start to time of observation
  mutate(dt_trans_start_gmt = as_datetime(paste(as_date(dt_trans_start_gmt),
                                                time_obs, " "))) %>% 
  #Remove time column
  select(-time_obs) %>% 
  mutate(basisOfRecord = "HUMAN_OBSERVATION",
         source = "Filchner")
names(fil_data) <- new_names

#GBIF data
gbif_data <- read.csv("Biological_Data/Cleaned_Data/GBIF_rgbif_cleaned.csv")  %>% 
  unite("source", institutionCode, collectionCode, remove = T) %>% 
  #Select columns of interest
  select(eventDate, decimalLatitude, decimalLongitude, individualCount, basisOfRecord, source) %>% 
  mutate(eventDate = as_date(eventDate),
         source = paste0("GBIF_", source))
names(gbif_data) <- new_names

#MEOP data
meop_data <- read_csv("Biological_Data/Cleaned_Data/MEOP_cleaned.csv") %>% 
  #Adding individual count column
  mutate(number_ind = 1) %>% 
  #Select columns of interest
  select(JULD, LATITUDE, LONGITUDE, number_ind, basisOfRecord) %>% 
  mutate(source = "MEOP")
names(meop_data) <- new_names

#OBIS data
obis_data <- read.csv("Biological_Data/Cleaned_Data/OBIS_cleaned.csv") %>% 
  unite("source", dataset_id, institutioncode, remove = T) %>% 
  #Select columns of interest
  select(eventdate, decimallatitude, decimallongitude, individualcount, basisofrecord, source) %>% 
  mutate(eventdate = as_date(eventdate),
         basisofrecord = case_when(basisofrecord == "Humanobservation" ~ "HUMAN_OBSERVATION"),
         source = paste0("OBIS_", source))
names(obis_data) <- new_names

#SCAR APIS data
scar_data <- read_csv("Biological_Data/Cleaned_Data/SCAR-APIS_cleaned.csv") %>% 
  #Select columns of interest
  select(eventDate, decimalLatitude, decimalLongitude, individualCount, basisOfRecord) %>% 
  mutate(source = "SCAR_APIS")
names(scar_data) <- new_names

#SCAR Biology
scar_bio <- read_csv("Biological_Data/Cleaned_Data/SCAR-Biology_cleaned.csv") %>% 
  #Select columns of interest
  select(eventDate, decimalLatitude, decimalLongitude, individualCount, basisOfRecord) %>% 
  mutate(source = "SCAR_Biology")
names(scar_bio) <- new_names


# Merging all datasets ----------------------------------------------------
merged_data <- bind_rows(belgica, bornemann, emage_data, fil_data, gbif_data,
                         meop_data, obis_data, scar_data, scar_bio, scar_seabirds) %>% 
  distinct(event_date, latitude, longitude, .keep_all = T) %>% 
  mutate(year = year(event_date), 
         month = month(event_date),
         season_year = case_when(month == 12 | month <= 2 ~ "summer",
                                 month >= 3 | month <= 5 ~ "autumn",
                                 month >= 6 | month <= 8 ~ "winter",
                                 month >= 9 | month <= 11 ~ "spring"),
         life_stage = case_when(month == 9 | month == 10 ~ "breeding",
                                month == 11 | month == 12 ~ "weaning",
                                month == 1 | month == 2 ~ "moulting",
                                T ~ "in-between"))

#Classifying data per decade
decades <- seq(1970, 2020, 10)
merged_data$decade<- decades[findInterval(merged_data$year, decades)]

#Removing single datasets
rm(belgica, bornemann, emage_data, fil_data, gbif_data, meop_data, obis_data,
   scar_data, scar_bio, scar_seabirds)

#Plotting data
world <- rnaturalearth::ne_countries(returnclass = 'sf')

merged_data %>% 
  ggplot(aes(longitude, latitude))+
  geom_point(aes(colour = source))+
  geom_sf(inherit.aes = F, data = world)+
  lims(y = c(-90, -45))+
  facet_wrap(~decade)

table(merged_data$decade)
hist(merged_data$decade)

table(merged_data$month)
hist(merged_data$month)

#Saving merged dataset
merged_data %>% 
  write_csv("Biological_Data/Cleaned_Data/All_sources_clean_data.csv")
