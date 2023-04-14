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
obis_data <- read.csv("Data/OBIS/Occurrence.csv", na.strings = c("", " ", NA)) %>% 
  #Removing date start/mid/end because it is based on eventDate column and
  #removing repeated columns with scientific name
  select(!c(date_start:originalscientificname, marine:species)) %>% 
  #Removing any empty columns
  janitor::remove_empty("cols") %>% 
  #Removing observations for which there is no date
  drop_na(eventdate) %>% 
  #Removing observations for preserved specimens
  filter(basisofrecord != "PreservedSpecimen") %>% 
  #Removing observations from Iziko SA Museum incorrectly labelled as "human observation"
  filter(dataset_id != "2e471257-d3af-40f7-ba72-fccdb460fc2b") %>% 
  #Removing dataset from Terra Nova Expedition in 1910-1913
  filter(dataset_id != "159421b6-a764-4af1-9cc3-bc652dbeb12e") %>% 
  #Removing dataset from US Antarctic expeditions 1939-1941
  filter(dataset_id != "2888968e-dc59-4ae4-873b-236b0c7a3e5a") %>% 
  #Removing datasets downloaded directly from original source
  #Belgica 121 expedition
  filter(!dataset_id %in% c("7eb2ccad-2063-4ea8-b4e3-164dbc852fce", 
                            #ARGOS tracking
                            "8a1109fe-bf82-42b2-b296-a46c8059b82b",
                            #SCAR RAATD
                            "48cb8624-a221-47ed-9a6d-b99b0bb394e0",
                            #Marine birds and mammals of the SO (BELARE)
                            "7350daa3-7c45-4f0f-b897-32b2a2b9c241",
                            #Seabirds of the Southern and South Indian Ocean
                            "832a055e-49ee-40bd-b138-c29cc07606d7")) %>% 
  #Removing any records dropped by OBIS quality control
  filter(dropped == 'f') %>% 
  #Ensuring basis of record categories are recorded in the same way
  mutate(basisofrecord = str_to_title(basisofrecord)) %>% 
  #Removing any observations north of 45S as it is beyond the Southern Ocean boundaries
  filter(decimallatitude <= -45) %>% 
  #Removing known default values for coordinate uncertainty
  filter(!coordinateuncertaintyinmeters %in% c(301, 3036, 999, 9999)) %>% 
  #Removing observations with low coordinate precision. We chose to remove observations with precision over
  #10 km because this is the nominal horizontal resolution of our environmental data
  filter(coordinateuncertaintyinmeters <= 10000 | is.na(coordinateuncertaintyinmeters)) %>% 
  #Correcting date format
  #Ensuring year, month and date columns are filled in for all observations with an event date
  mutate(eventdate = parse_date_time(eventdate, c("%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S")), 
         year = case_when(is.na(year) ~ year(eventdate), 
                          T ~ year),
         month = case_when(is.na(month) ~ month(eventdate), 
                           T ~ month),
         day = case_when(is.na(day) ~ day(eventdate), 
                         T ~ day)) %>%
  #Removing records prior to 1968 (beginning of modelled environmental data)
  filter(year >= 1968) %>% 
  #Removing observation for which date could not be parsed because a date range was provided
  #instead of a single date
  drop_na(eventdate) %>% 
  #Removing observations with no information in the "individual counts" column
  drop_na(individualcount) %>% 
  #Removing duplicate observations - based on lat/lon coordinates and date of observation
  distinct(eventdate, decimallatitude, decimallongitude, .keep_all = T) %>% 
  #Adding species column before applying coordinate cleaner workflow
  mutate(species = "Lobodon carcinophagus")

#Second filter - Coordinate cleaner
obis_CC <- obis_data %>% 
  #Removing records with invalid coordinates and potential outliers
  clean_coordinates(lon = "decimallongitude", lat = "decimallatitude") %>% 
  #Checking potential issues with coordinate conversions and rounding
  cd_ddmm(lon = "decimallongitude", lat = "decimallatitude", ds = "dataset_id") %>% 
  #Removing duplicated records - based on coordinates and date of observation
  cc_dupl(lon = "decimallongitude", lat = "decimallatitude", additions = c("year", "month", "day"))

#Checking results
summary(obis_CC)

#Removing points with equal latitude and longitude coordinates (highlighted in .equ column)
obis_CC <- obis_CC %>% 
  filter(.equ == T)

#Saving clean dataset 
obis_CC %>% 
  write_csv("Cleaned_Data/OBIS_cleaned.csv")

