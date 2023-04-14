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
  #Removing rows where no bin data is available
  filter(comment != "No bin specific data available" | is.na(comment)) %>% 
  #Breaking down census code into day, transect and section to complete coordinates
  separate_wider_delim(census_code, "_", names = c("census_day", "census_transect", "census_segment"), cols_remove = F) %>% 
  #Group by census day and transect
  group_by(census_day, census_transect) %>%
  #Fill start and end coordinates within grouped data
  fill(lat_seg_start:lon_seg_end, .direction = "downup") %>% 
  #Removing rows with no coordinate information
  drop_na(lat_seg_start:lon_seg_end) %>% 
  #Removing rows with no date information
  drop_na(dt_seg_start_gmt, dt_seg_end_gmt) %>% 
  #Removing rows identified as containing erroneous coordinates
  filter(!comment == "Erraneous Lat Lon data" | is.na(comment)) %>% 
  #Removing observations from last bin starting at 345 m from the helicopter up to the horizon
  #because reliability of species ID diminishes with distance from helicopter
  select(!number_ind_bin_6) %>% 
  #Ensuring total number of individuals per day/transect/segment are correct
  rowwise() %>% 
  mutate(number_ind_tot = sum(c_across(number_ind_bin_1:number_ind_bin_5)))

#Calculate mid point for all transects
ant_census <- geosphere::midPoint(p1 = ant_census[, c("lon_seg_start", "lat_seg_start")],
                    p2 = ant_census[, c("lon_seg_end", "lat_seg_end")]) %>% 
  as_tibble() %>% 
  #Rename column outputs
  rename(lon_seg_mid = lon, lat_seg_mid = lat) %>% 
  #Attaching to clean dataset
  bind_cols(ant_census, .) %>% 
  #Removing empty columns
  janitor::remove_empty("cols") %>% 
  #Adding basis of record column
  mutate(basisOfRecord = "HUMAN_OBSERVATION")


# Saving clean data -------------------------------------------------------
ant_census %>% 
  write_csv("Cleaned_Data/ANT-XV_3_seal_census_clean_data.csv")
  