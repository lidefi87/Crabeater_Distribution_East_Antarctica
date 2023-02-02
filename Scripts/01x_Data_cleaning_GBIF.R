###############################################################################
# Cleaning crabeater seal data obtained from GBIF using rgbif
# Author: Denisse Fierro Arcos
# Date: 2023-01-27
# 
# Crabeater seal data downloaded from GBIF using rgif package.


# Loading libraries -------------------------------------------------------
library(rgbif)
library(tidyverse)


# Searching data in GBIF --------------------------------------------------
#Accessing metadata for crabeater seals
crabeater <- name_backbone("Lobodon carcinophaga")

#Querying GBIF database
query <- occ_download(pred("taxonKey", crabeater$usageKey),+
                      #Keeping only records that are georeferenced
                      pred("hasCoordinate", TRUE))
#Starting download
occ_download_wait(query)

#Downloading all data to local disk - zip file is downloaded
crabeater_query <- occ_download_get(key = query, path = "Data/GBIF/rgbif_download/") %>% 
  occ_download_import(crabeaters_gbif, path = "Data/GBIF/rgbif_download/", na.strings = c("", NA))


# Cleaning data -----------------------------------------------------------
crabeater <- crabeater_query %>% 
  #Removing absent records
  filter(occurrenceStatus == "PRESENT") %>% 
  #Removing any observations north of 45S as it is beyond the Southern Ocean boundaries
  filter(decimalLatitude <= -45) %>% 
  #Removing any fossil or preserved specimen records
  filter(!basisOfRecord %in% c("PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN")) %>% 
  #Removing any records with reported georeference issues 
  filter(hasGeospatialIssues == F) %>% 
  #Removing known default values for coordinate uncertainty
  filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 9999)) %>% 
  #Removing any records of dead or mummified animals
  filter(!str_detect(str_to_lower(occurrenceRemarks), "deceased|dead|died|mumm") | is.na(occurrenceRemarks)) %>% 
  #Removing any records identified as "approximate coordinates"
  filter(!str_detect(str_to_lower(occurrenceRemarks), 'coord.* approx') | is.na(occurrenceRemarks)) %>% 
  #Reclassifying bio-logging locations as 'MACHINE_OBSERVATIONS'
  mutate(basisOfRecord = case_when(str_detect(samplingProtocol, 'bio-log') & basisOfRecord != 'MACHINE_OBSERVATION' ~ 'MACHINE_OBSERVATION',
                                   T ~ basisOfRecord)) %>%
  #Changing observations with zeroes (0) in individual counts to 1. One record has life stage included and other records
  #belong to AAD ANARE database which only contains sightings and no absences (https://data.aad.gov.au/metadata/records/DB_Historic_WoV)
  mutate(individualCount = case_when(individualCount == 0 ~ 1, 
                                     T ~ individualCount)) %>% 
  #Removing observations from collection "ANTXXIII" because coordinates are all integers
  #The error associated to the reported location can be as large as 120 km
  filter(!str_detect(collectionCode, "ANTXXIII") | is.na(collectionCode)) %>% 
  #Removing observations with uncertainty about species
  filter(!str_detect(str_to_lower(occurrenceRemarks), 'seal.*\\?') | is.na(occurrenceRemarks)) %>% 
  #Removing GBIF record #1563413650 because there is no date associated to observation
  filter(gbifID != 1563413650) %>% 
  #Removing observation with "PRESUMED_NEGATED_LATITUDE" as issue because it is too far inland
  filter(!str_detect(issue, "PRESUMED_NEGATED_LATITUDE") | is.na(issue)) %>% 
  #Adding year 2008 to Belgian expedition as that is the year the expedition took place
  mutate(year = case_when(collectionCode == "Belgian Antarctic Research Expedition 2008/09" ~ 2008, T ~ year)) %>%
  
  #Removing any empty columns
  janitor::remove_empty(which = "cols")

crabeater %>% write_csv("Cleaned_Data/GBIF_rgbif_cleaned.csv")

# Dividing data between 

