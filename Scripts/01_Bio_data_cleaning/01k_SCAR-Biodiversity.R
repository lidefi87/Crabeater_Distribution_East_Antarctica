###############################################################################
# Cleaning crabeater seal data obtained from RAATD tracking data
# Author: Denisse Fierro Arcos
# Date: 2023-02-10
# 
# Crabeater seal data downloaded from SCAR Biodiversity
# Dataset citation: Ropert-Coudert Y, Van de Putte A P, Bornemann H, Charrassin 
# J, Costa D P, Danis B, Hückstädt L A, Jonsen I D, Lea M, Reisinger R R, 
# Thompson D, Torres L G, Trathan P N, Wotherspoon S, Ainley D G, Alderman R, 
# Andrews-Goff V, Arthur B, Ballard G, Bengtson J, Bester M N, Boehme L, Bost C, 
# Boveng P, Cleeland J, Constantine R, Crawford R J M, Dalla Rosa L, de Bruyn P 
# N, Delord K, Descamps S, Double M, Emmerson L, Fedak M, Friedlander A, Gales 
# N, Goebel M, Goetz K T, Guinet C, Goldsworthy S D, Harcourt R, Hinke J, 
# Jerosch K, Kato A, Kerry K R, Kirkwood R, Kooyma G L, Kovacs K M, Lawton K, 
# Lowther A D, Lydersen C, Lyver P O, Makhado A B, Márquez M E I, McDonald B, 
# McMahon C, Muelbert M, Nachtsheim D, Nicholls K W, Nordøy E S, Olmastroni S, 
# Phillips R A, Pistorius P, Plötz J, Pütz K, Ratcliffe N, Ryan P G, Santos M, 
# Schytte Blix A, Southwell C, Staniland I, Takahashi A, Tarroux A, Trivelpiece 
#W, Wakefield E, Weimerskirch H, Wienecke B, Xavier J C, Raymond B, Hindell M A 
# (2020): The Retrospective Analysis of Antarctic Tracking (Standardised) Data 
# from the Scientific Committee on Antarctic Research. v1.3. SCAR - AntOBIS. 
# Dataset/Metadata. 
# https://ipt.biodiversity.aq/resource?r=raatd_scar_trackingdata&v=1.3

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(CoordinateCleaner)

# Loading data downloaded from SCAR ---------------------------------------
scar_bio <- read_csv("Original_Data/SCAR_biodiversity/bb903809-5652-4afd-a326-b64064dbdbf1.csv") %>%
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
  write_csv("Biological_Data/Cleaned_Data/SCAR-Biology_cleaned.csv")
