###############################################################################
# Cleaning crabeater seal data obtained from Colin Southwell - AAD
# Author: Denisse Fierro Arcos
# Date: 2023-02-07
# 
# Crabeater seal data downloaded from Colin Southwell - AAD

# Loading libraries -------------------------------------------------------
library(tidyverse)
library(CoordinateCleaner)

# Loading data downloaded from SCAR ---------------------------------------
southwell_Feb85 <- read_delim("Data/Southwell_AAD/SpaceTimeDist_Ship_Effort_85_99_8Sep03_03.txt") %>% 
  
southwell_Sep85 <- read_delim("Data/Southwell_AAD/SpaceTimeDist_Ship_Sightings_85_99_8Sep03_05.txt") %>% 
  #Extract crabeater seals only
