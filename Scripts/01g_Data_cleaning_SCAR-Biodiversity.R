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
scar_bio <- read_csv("Data/SCAR_biodiversity/07b4986f-a39b-4443-b87e-10dd8bce5203.csv") %>% 
  #Selecting only observations for crabeater seals
  filter(str_detect(scientificName, "Lobodon carcinophaga")) 

scar_bio2 <- read_csv("Data/SCAR_biodiversity/bb903809-5652-4afd-a326-b64064dbdbf1.csv") %>% 
  filter(str_detect(scientificName, "Lobodon carcinophaga"))

scar_bio <- scar_bio %>% 
  bind_rows(scar_bio2) %>% 
  #Removing any observations for which occurrence status is not present
  filter(occurrenceStatus == "PRESENT") %>% 
  #Since there is no information in the "individual count" column,
  #we assume there was at least one animal present
  mutate(individualCount = 1) %>% 
  #Remove empty columns
  janitor::remove_empty("cols") %>% 
  #Removing observations for which there is no date
  drop_na(eventDate) %>% 
  #Removing any observations north of 45S as it is beyond the Southern Ocean boundaries
  filter(decimalLatitude <= -45) %>%
  #Removing duplicate observations - based on lat/lon coordinates and date of observation
  distinct(eventDate, decimalLatitude, decimalLongitude, .keep_all = T) %>% 
  #Creating one column with the name of the species
  unite(species, genus, specificEpithet, sep = " ")

#Remove second dataset because it has been merged already
rm(scar_bio2)

#Second filter - Coordinate cleaner
scar_bio_CC <- scar_bio %>% 
  #Removing records with invalid coordinates and potential outliers
  clean_coordinates(lon = "decimalLongitude", lat = "decimalLatitude")
#Checking results
summary(scar_bio_CC)
#No problematic data points detected

#Final check for coordinate issues
scar_bio_CC <- scar_bio %>% 
  #Checking potential issues with coordinate conversions and rounding
  cd_ddmm(lon = "decimalLongitude", lat = "decimalLatitude", ds = "datasetKey") %>% 
  #Removing duplicated records - based on coordinates and date of observation
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", additions = c("eventDate"))

#Saving clean dataset 
scar_bio_CC %>% 
  write_csv("Cleaned_Data/SCAR-Biology_cleaned.csv")
