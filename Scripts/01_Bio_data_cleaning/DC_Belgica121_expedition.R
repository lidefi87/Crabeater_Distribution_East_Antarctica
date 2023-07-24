###############################################################################
# Cleaning crabeater seal data obtained from Belgica 121 expedition
# Author: Denisse Fierro Arcos
# Date: 2023-04-13
# 
# Data downloaded from GBIF.
# Dataset citation: Danis, B., Christiansen, H., Guillaumot, C., Heindler, F., 
# Houston, R., Jossart, Q., Lucas, K., Moreau, C., Pasotti, F., Robert. H., 
# Wallis, B., Saucède, T., [2021]. The Belgica 121 expedition to the 
# Western Antarctic Peninsula: a high resolution biodiversity census. 
# https://doi.org/10.15468/56bv6z.


# Loading libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)
library(CoordinateCleaner)

# Loading data ------------------------------------------------------------
belgica <- read_delim("Data/Belgica121Expedition/occurrence.txt") %>% 
  #Selecting observations for crabeater seals
  filter(genus == "Lobodon") %>% 
  #Removing any records with reported georeference issues 
  filter(hasGeospatialIssues == F) %>% 
  #Removing observations with low coordinate precision. We chose to remove observations with precision over
  #10 km because this is the nominal horizontal resolution of our environmental data
  filter(coordinateUncertaintyInMeters <= 10000 | is.na(coordinateUncertaintyInMeters)) %>% 
  #Removing empty columns
  janitor::remove_empty("cols") %>% 
  #Keeping columns with useful information
  select(c(gbifID, basisOfRecord, occurrenceID, individualCount, 
           eventDate:samplingProtocol, decimalLatitude:decimalLongitude, datasetKey))

#Saving clean dataset 
belgica %>% 
  write_csv("Cleaned_Data/Belgica121_cleaned.csv")

