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
  filter(!coordinateUncertaintyInMeters %in% c(301,3036,999,9999)) %>% 
  #Removing any empty columns
  janitor::remove_empty(which = "cols")

crabeater %>% write_csv("Cleaned_Data/GBIF_rgbif_cleaned.csv")
