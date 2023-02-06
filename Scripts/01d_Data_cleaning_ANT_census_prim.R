###############################################################################
# Cleaning crabeater seal data obtained from different sources
# Author: Denisse Fierro Arcos
# Date: 2023-02-06
# 
# Open source data obtained from PANGEA.

# Loading libraries -------------------------------------------------------
library(tidyverse)

# Loading individual datasets ---------------------------------------------
# EMAGE-II seal census (Gurarie_2016 folder) ------------------------------
emage_census <- read_delim("Data/Gurarie_2016/datasets/EMAGE-II_seal_census_raw-data.tab", 
                         #Skipping the first 35 rows containing metadata
                         delim = "\t", skip = 35) %>% 
  #Shortening column names
  janitor::clean_names() 

#Keeping only columns with crabeater seal abundance and where no animals were detected
crab_emage <- emage_census %>% 
  filter(comment == "Lobodon carcinophaga" | comment == "Lobodon carcinophaga, subadult") %>% 
  bind_rows(emage_census %>% filter(seals_number == 0))


# Saving clean data -------------------------------------------------------
crab_emage %>% 
  write_csv("Cleaned_Data/EMAGE-II_seal_census_clean_data.csv")
