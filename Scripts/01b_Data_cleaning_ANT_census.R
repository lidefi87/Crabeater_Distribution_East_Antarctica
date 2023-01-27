###############################################################################
# Cleaning crabeater seal data obtained from different sources
# Author: Denisse Fierro Arcos
# Date: 2023-01-27
# 
# Open source data obtained from PANGEA and SCAR. Colin Southwell from AAD
# shared data collected by him and his team in 1985
# Each data set is cleaned separately before they are combined into a single
# dataset.


# Loading libraries -------------------------------------------------------
library(tidyverse)

# Loading individual datasets ---------------------------------------------
# ANT-XV_3 seal census (Bester_2015 folder) -------------------------------
ant_census <- read_delim("Data/Bester_2015/datasets/ANT-XV_3_seal_census_raw_data.tab", 
                             #Skipping the first 82 rows containing metadata
                             delim = "\t", skip = 82) %>% 
  #Shortening column names
  janitor::clean_names() %>% 
  #Keeping only columns with crabeater seal abundance and supporting data
  select_at(vars(id_census_leg_frame_each_day_is:side_observation_from_left_or_righ, contains("carcino"), ice_cov_percent:comment)) 

#Simplify column names
col_names_simple <- c("census_code", "lat_seg_start", "lon_seg_start", 
                      "lat_seg_end", "lon_seg_end", "dt_seg_start_gmt", 
                      "dt_seg_end_gmt", "time_seg_start_gmt", "time_seg_end_gmt", 
                      "seg_length_km", "altitude_m", "speed_km_h", "bin_width",
                      "cum_seg_length_km", "side_helicopter", "number_ind_tot",
                      "number_ind_bin_1", "number_ind_bin_2", "number_ind_bin_3", 
                      "number_ind_bin_4", "number_ind_bin_5", "number_ind_bin_6", 
                      "sic_per", "ice_floe_per", "brash_ice_per", "cake_ice_per",
                      "small_floe_10-100m_per", "med_floe_100-500m_per",
                      "large_floe_500m_per", "air_temp_trans_C", "wind_speed_trans_ms",
                      "visibility", "glare", "contrast", "heli_type", "comment")
                       
#Renaming columns
names(ant_census) <- col_names_simple

#Cleaning data
ant_census <- ant_census %>% 
  #Removing census 1 because it was a test flight
  filter(!str_detect(census_code, "^1_")) %>% 
  #Ensuring counts with NA recorded are changed to zeroes (0)
  mutate_at(vars(contains("ind_bin")), ~replace_na(., 0))
  
  
  
  
  #Ensuring all unidentified animals are spelled the same way
  mutate(species_lobodon_carcinophaga_or_lepto = str_to_sentence(species_lobodon_carcinophaga_or_lepto),
         #Changing NA values under individuals to zero (0)
         ind_no_number_size_group = replace_na(ind_no_number_size_group, 0)) 