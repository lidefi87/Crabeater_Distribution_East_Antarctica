###############################################################################
# Cleaning crabeater seal data from tracking data - 13 tagged seals
# Author: Denisse Fierro Arcos
# Date: 2023-04-20
# 
# Open source data obtained from PANGAEA.
# Dataset citation: Bornemann, Horst; Plötz, Joachim (2013): Occurrence of 13 
# crabeater seals (Lobodon carcinophaga) from ARGOS satellite track observations
# during campaign DRE1998. Alfred Wegener Institute, Helmholtz Centre for Polar 
# and Marine Research, Bremerhaven, PANGAEA, https://doi.org/10.1594/PANGAEA.819980


# Loading libraries -------------------------------------------------------
library(tidyverse)


# Load data ---------------------------------------------------------------
lob_argos <- read_delim("Original_Data/Bornemann_2013/DRE1998_ARGOS_sat_tracks.tab", 
                        skip = 40) %>% 
  janitor::clean_names() %>% 
  #Adding column to classify data as obtained by satellite tags
  mutate(basisOfRecord = "MACHINE_OBSERVATION")


# Saving dataset ----------------------------------------------------------
write_csv(lob_argos, "Biological_Data/Cleaned_Data/Bornemann_ARGOS.csv")


