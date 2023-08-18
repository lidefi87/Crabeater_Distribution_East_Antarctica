###############################################################################
# Cleaning crabeater seal data obtained from GBIF
# Author: Denisse Fierro Arcos
# Date: 2023-01-27
# 
# Crabeater seal data downloaded from GBIF using rgbif package.
# Dataset citation: GBIF Occurrence Download https://doi.org/10.15468/dl.6ea6vv.
# Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2023-01-27.


# Loading libraries -------------------------------------------------------
library(rgbif)
library(tidyverse)
library(lubridate)
library(CoordinateCleaner)

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
crabeater_query <- occ_download_get(key = query, path = "Original_Data/GBIF/rgbif_download/") %>% 
  occ_download_import(crabeaters_gbif, path = "Original_Data/GBIF/rgbif_download/", na.strings = c("", NA))

#Getting citation details
gbif_citation(occ_download_meta("0260596-220831081235567"))$download

# Cleaning data -----------------------------------------------------------
#Loading data from disk
crabeater <- read.delim("Original_Data/GBIF/occurrence.txt", na.strings = c(" ", "", NA)) %>% 
  #Removing absent records
  filter(occurrenceStatus == "PRESENT") %>% 
  #Removing any observations north of 45S as it is beyond the Southern Ocean boundaries
  filter(decimalLatitude <= -45) %>% 
  #Removing any fossil or preserved specimen records
  filter(!basisOfRecord %in% c("PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN")) %>% 
  #Removing any records with reported georeference issues 
  filter(hasGeospatialIssues == 'false') %>% 
  #Removing known default values for coordinate uncertainty
  filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 9999)) %>% 
  #Removing records from unknown sources
  drop_na(publisher) %>% 
  #Removing any records of dead or mummified animals
  filter(!str_detect(str_to_lower(occurrenceRemarks), "deceased|dead|died|mumm") | is.na(occurrenceRemarks)) %>% 
  #Removing any records identified as "approximate coordinates"
  filter(!str_detect(str_to_lower(occurrenceRemarks), "coord.* approx") | is.na(occurrenceRemarks)) %>% 
  #Removing any records when latitude and longitude coordinates are exactly the same as this usually indicates
  #data entry errors
  filter(decimalLatitude != decimalLongitude) %>% 
  #Reclassifying bio-logging locations as 'MACHINE_OBSERVATIONS' 
  mutate(basisOfRecord = case_when(str_detect(samplingProtocol, "bio-log") & basisOfRecord != "MACHINE_OBSERVATION" ~ "MACHINE_OBSERVATION",
                                   T ~ basisOfRecord)) %>%
  #Reclassifying observations from tagged seal program by SCAR - AntOBIS as 'MACHINE_OBSERVATIONS'
  #More info: https://obis.org/dataset/089ef71a-b32c-4e35-8c74-79cde58dbfff
  mutate(basisOfRecord = case_when(str_detect(identifier, "tagged") & basisOfRecord != "MACHINE_OBSERVATION" ~ "MACHINE_OBSERVATION",
                                   T ~ basisOfRecord)) %>% 
  #Reclassifying observations from Bornemann as 'MACHINE_OBSERVATIONS' because this data was
  #obtained from tagged seals (https://doi.org/10.1594/PANGAEA.819980)
  mutate(basisOfRecord = case_when(str_detect(datasetID, "PANGAEA.819980") & basisOfRecord != "MACHINE_OBSERVATION" ~ "MACHINE_OBSERVATION",
                                   T ~ basisOfRecord)) %>% 
  #Removing duplicate entry from iNaturalist
  filter(!gbifID == 3355056489) %>% 
  #Removing entries from grey literature
  filter(!str_detect(publisher, "Ukrainian")) %>% 
  #Removing observations where individual number reported did not match status (i.e, presence reported as 0 individuals)
  filter(!str_detect(issue, "COUNT_CONFLICTS") | is.na(issue)) %>% 
  #Removing observations from collection "ANTXXIII" because coordinates are all integers
  #The error associated to the reported location can be as large as 120 km
  filter(!str_detect(collectionCode, "ANTXXIII") | is.na(collectionCode)) %>% 
  #Removing observations with uncertainty about species - indicated by question mark
  filter(!str_detect(str_to_lower(occurrenceRemarks), 'seal.*\\?|crabeater.*\\?') | is.na(occurrenceRemarks)) %>% 
  #Removing observation with "PRESUMED_NEGATED_LATITUDE" as issue because it is too far inland
  filter(!str_detect(issue, "PRESUMED_NEGATED_LATITUDE") | is.na(issue)) %>% 
  #Setting date column to date format
  mutate(eventDate = as_date(eventDate)) %>%
  #Belgian expedition took place between 2008-12-07 and 2009-01-14. Records have a date range, rather than
  #date of observation. Updating all dates to 2008-12-26 as it represents half point of expedition
  mutate(eventDate = case_when(str_detect(collectionCode, "Belgian") ~ as_date("2008-12-26"),
                               T ~ eventDate),
         #Ensuring year, month and date columns are filled in for all observations with an event date
         year = case_when(is.na(year) ~ year(eventDate), 
                          T ~ year),
         month = case_when(is.na(month) ~ month(eventDate), 
                           T ~ month),
         day = case_when(is.na(day) ~ day(eventDate), 
                T ~ day)) %>%
  #Removing any records with no date
  drop_na(eventDate) %>% 
  #Removing records prior to 1968 (beginning of modelled environmental data)
  filter(year >= 1968) %>% 
  #Removing duplicate observations - based on lat/lon coordinates and date of observation
  distinct(eventDate, decimalLatitude, decimalLongitude, .keep_all = T) %>% 
  #Removing observations with low coordinate precision. We chose to remove observations with precision over
  #10 km because this is the nominal horizontal resolution of our environmental data
  filter(coordinateUncertaintyInMeters <= 10000 | is.na(coordinateUncertaintyInMeters)) %>% 
  #Removing datasets downloaded directly from original source
  #Belgica 121 expedition
  filter(!datasetKey %in% c("b635be2e-76ea-4600-8f83-549601653c0a", 
                            #SCAR APIS
                            "82ec8f46-f762-11e1-a439-00145eb45e9a",
                            #SCAR RAATD
                            "b86fe411-8e62-4cd0-aab2-914d75401598",
                            #Seabirds of the Southern and South Indian Ocean
                            "82dd797a-f762-11e1-a439-00145eb45e9a",
                            #ARGOS tracking from Bornemann
                            "8ea94e2a-118c-4a53-a957-5348d48d61e1")) %>% 
  #Given that all records included in this dataset are based on observations,
  #the "individualCount" column must have at least a value of 1. It cannot contain NAs
  mutate(individualCount = case_when(is.na(individualCount) ~ 1,
                                     T ~ individualCount)) %>% 
  #Removing any empty columns
  janitor::remove_empty(which = "cols")

#Second filter - Coordinate cleaner
crabeater_CC <- crabeater %>% 
  #Removing records with invalid coordinates and potential outliers
  clean_coordinates(lon = "decimalLongitude", lat = "decimalLatitude",
                    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "zeros")) %>% 
  #Checking potential issues with coordinate conversions and rounding
  cd_ddmm(lon = "decimalLongitude", lat = "decimalLatitude", ds = "datasetKey") %>% 
  #Removing duplicated records - based on coordinates and date of observation
  cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", additions = c("year", "month", "day"))

#No problematic data points detected
summary(crabeater_CC)

#Saving clean dataset 
crabeater_CC %>% 
  write_csv("Biological_Data/Cleaned_Data/GBIF_rgbif_cleaned.csv")

