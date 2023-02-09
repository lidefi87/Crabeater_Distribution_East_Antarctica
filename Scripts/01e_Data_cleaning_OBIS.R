###############################################################################
# Cleaning crabeater seal data obtained from OBIS
# Author: Denisse Fierro Arcos
# Date: 2023-02-06
# 
# Crabeater seal data downloaded from OBIS.

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)
library(CoordinateCleaner)

# Loading data downloaded from OBIS ---------------------------------------
obis_data <- read_csv("Data/OBIS/Occurrence.csv") %>% 
  #Removing date start/mid/end because it is based on eventDate column
  select(!date_start:date_end) %>% 
  #Removing any empty columns
  janitor::remove_empty("cols") %>% 
  #Removing observations for which there is no date
  drop_na(eventdate) %>% 
  #Removing observations for preserved specimens
  filter(basisofrecord != "PreservedSpecimen") %>% 
  #Ensuring basis of record categories are recorded in the same way
  mutate(basisofrecord = str_to_title(basisofrecord)) %>% 
  #Removing any observations north of 45S as it is beyond the Southern Ocean boundaries
  filter(decimallatitude <= -45) %>% 
  #Retaining only presence records
  filter(occurrencestatus == "present") %>% 
  #Removing known default values for coordinate uncertainty
  filter(!coordinateuncertaintyinmeters %in% c(301, 3036, 999, 9999)) %>% 
  #Removing observations with low coordinate precision. We chose to remove observations with precision over
  #10 km because this is the nominal horizontal resolution of our environmental data
  filter(coordinateuncertaintyinmeters <= 10000 | is.na(coordinateuncertaintyinmeters)) %>% 
  #Removing records prior to 1968 (beginning of modelled environmental data)
  filter(date_year >= 1968) %>% 
  #Standardising date formats
  mutate(dateidentified = case_when(is.na(dateidentified) ~ parse_datetime(eventdate), 
                                    T ~ dateidentified),
         #Ensuring year, month and date columns are filled in for all observations with an event date
         year = case_when(is.na(year) ~ year(eventdate), 
                          T ~ year),
         month = case_when(is.na(month) ~ month(eventdate), 
                           T ~ month),
         day = case_when(is.na(day) ~ day(eventdate), 
                         T ~ day)) %>%
  #Removing observation for which date could not be parsed because a date range was provided
  #instead of a single date
  drop_na(dateidentified) %>% 
  #Removing observations with no information in the "individual counts" column
  drop_na(individualcount) %>% 
  #Removing duplicate observations - based on lat/lon coordinates and date of observation
  distinct(eventdate, decimallatitude, decimallongitude, .keep_all = T) %>% 
  #Removing columns that are not needed
  select(!c(originalscientificname:maximumdepthinmeters, taxonrank,
            superdomain:infraorder, datasetid, verbatimeventdate)) %>% 
  #Removing any empty columns
  janitor::remove_empty(which = "cols")

#Second filter - Coordinate cleaner
obis_CC <- obis_data %>% 
  #Removing records with invalid coordinates and potential outliers
  clean_coordinates(lon = "decimallongitude", lat = "decimallatitude")
#Checking results
summary(obis_CC)
plot(obis_CC)
#Flagged observations are not removed because crabeaters can occur in the WAP

#Final check for coordinate issues
obis_CC <- obis_data %>% 
  #Checking potential issues with coordinate conversions and rounding
  cd_ddmm(lon = "decimallongitude", lat = "decimallatitude", ds = "dataset_id") %>% 
  #Removing duplicated records - based on coordinates and date of observation
  cc_dupl(lon = "decimallongitude", lat = "decimallatitude", additions = c("dateidentified"))

#Saving clean dataset 
obis_CC %>% 
  write_csv("Cleaned_Data/OBIS_cleaned.csv")

