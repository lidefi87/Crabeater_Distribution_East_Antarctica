###############################################################################
# Cleaning crabeater seal data obtained from SCAR Biodiversity
# Author: Denisse Fierro Arcos
# Date: 2023-02-10
# 
# Crabeater seal data downloaded from SCAR Biodiversity

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(CoordinateCleaner)

# Loading data downloaded from SCAR ---------------------------------------
scar_bio <- read_csv("Data/SCAR_biodiversity/bb903809-5652-4afd-a326-b64064dbdbf1.csv") %>%
  #Selecting only observations for crabeater seals
  filter(str_detect(scientificName, "Lobodon carcinophaga")) %>% 
  #Remove empty columns
  janitor::remove_empty("cols") %>% 
  #Removing any observations for which occurrence status is not present
  filter(occurrenceStatus == "PRESENT") %>%
  #Removing any observations north of 45S as it is beyond the Southern Ocean boundaries
  filter(decimalLatitude <= -45) %>%
  #Adding "individual count" column and setting to one as observations refer to
  #location of a particular individual over time
  mutate(individualCount = 1) %>% 
  #Removing observations for which there is no date
  drop_na(eventDate) %>% 
  #Removing duplicate observations - based on lat/lon coordinates and date of observation
  distinct(eventDate, decimalLatitude, decimalLongitude, .keep_all = T) %>% 
  #Creating one column with the name of the species
  unite(species, genus, specificEpithet, sep = " ")

#Second filter - Coordinate cleaner
scar_bio_CC <- scar_bio %>% 
  #Removing records with invalid coordinates and potential outliers
  clean_coordinates(lon = "decimalLongitude", lat = "decimalLatitude") %>% 
  #Checking potential issues with coordinate conversions and rounding
  cd_ddmm(lon = "decimalLongitude", lat = "decimalLatitude", ds = "datasetKey") %>% 
  #Removing duplicated records - based on coordinates and date of observation
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", additions = c("year", "month", "day"))

#Checking results
summary(scar_bio_CC)

#Saving clean dataset 
scar_bio_CC %>% 
  write_csv("Cleaned_Data/SCAR-Biology_cleaned.csv")
