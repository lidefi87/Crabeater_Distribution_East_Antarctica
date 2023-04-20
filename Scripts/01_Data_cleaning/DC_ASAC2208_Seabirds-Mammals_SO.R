###############################################################################
# Cleaning crabeater seal data obtained from the Seabirds of the Southern and
# South Indian Ocean
# Author: Denisse Fierro Arcos
# Date: 2023-04-14
# 
# Full dataset downloaded from the Australian Antarctic Data Centre
# Dataset citation: Woehler, E. and Raymond, B. (1999) Distribution and 
# abundance of seabirds in the Southern Indian Ocean, 1978/1979+, Ver. 1, 
# Australian Antarctic Data Centre. doi:10.4225/15/5643E8C0743C2.

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(lubridate)
library(CoordinateCleaner)


# Loading data ------------------------------------------------------------
asac <- read_csv("Data/ASAC_2208_seabirds/ASAC_2208_seabirds.csv") %>% 
  #Keeping observations for crabeater seals only
  filter(genus == "Lobodon") %>% 
  #Remove any observations without date or coordinates
  drop_na(observation_date, longitude, latitude) %>% 
  #Removing empty columns
  janitor::remove_empty("cols") %>% 
  #Removing columns not needed
  select(!c(observation_date_hour:observation_date_time_zone, bird_age)) %>% 
  #Correcting number of individuals observed based on notes column
  mutate(species_count = case_when(is.na(species_count) ~ as.integer(str_extract(notes, "[0-9] ")), 
                                   T ~ species_count)) %>% 
  mutate(species_count = case_when(is.na(species_count) & str_detect(str_to_lower(notes), "crabeater") ~ 1, 
                                   T ~ species_count)) %>% 
  #Removing observations of zero individuals as the dataset description does not
  #state absences were recorded
  filter(species_count > 0) %>% 
  #Ensuring observations are within SO boundaries
  filter(latitude <= -45) %>% 
  #Removing any records when latitude and longitude coordinates are exactly the same as this usually indicates
  #data entry errors
  filter(latitude != longitude) %>% 
  #Adding basis of record column
  mutate(basisOfRecord = "HUMAN_OBSERVATION")


# Saving data -------------------------------------------------------------
asac %>% 
  write_csv("Cleaned_Data/Cleaned_ASAC_2208_seabirds.csv")
