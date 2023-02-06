###############################################################################
# Cleaning crabeater seal data obtained from different sources
# Author: Denisse Fierro Arcos
# Date: 2023-01-25
# 
# Open source data obtained from PANGEA and SCAR. Colin Southwell from AAD
# shared data collected by him and his team in 1985
# Each data set is cleaned separately before they are combined into a single
# dataset.


# Loading libraries -------------------------------------------------------
library(tidyverse)

# Loading individual datasets ---------------------------------------------
# Bester_2015 -----
#Simplify column names
col_names_simple <- c("transect_code", "dt_trans_start_gmt", "dt_trans_end_gmt", 
                      "lat_trans_start", "lon_trans_start", "lat_trans_end", 
                      "lon_trans_end", "area_name", "transect_id",
                      "transect_length_km", "altitude_m", "duration_trans_hms", 
                      "segment_id", "segment_label", "lat_seg_start", "lon_seg_start",
                      "lat_seg_end", "lon_seg_end", "segment_length_km",
                      "dt_seg_start_gmt", "dt_seg_end_gmt", "duration_seg_hms",
                      "observer", "side_helicopter", "sampling_method", "species",
                      "number_ind", "distance_interval", "distance_helip_m", 
                      "comment_inconsistencies", "comment_individual", "comment_dup",
                      "time_obs", "lat_obs", "lon_obs", "floe_class", "segment_sic",
                      "segment_si_surface", "air_temp_trans_C", "wind_speed_trans_ms",
                      "vis_trans_m", "comment")


# Filchner Ice Shelf Outflow dataset --------------------------------------
bester_outflow <- read_delim("Data/Bester_2015/datasets/FIL_2014_Filchner_Outflow_seal_census.tab", 
                             #Skipping the first 76 rows containing metadata
                             delim = "\t", skip = 76) %>% 
  #Shortening column names
  janitor::clean_names() %>% 
  #Ensuring all unidentified animals are spelled the same way
  mutate(species_lobodon_carcinophaga_or_lepto = str_to_sentence(species_lobodon_carcinophaga_or_lepto),
         #Changing NA values under individuals to zero (0)
         ind_no_number_size_group = replace_na(ind_no_number_size_group, 0)) 

bester_trough <- read_delim("Data/Bester_2015/datasets/FIL_2014_Filchner_Trough_seal_census.tab", 
                            #Skipping the first 57 rows containing metadata
                            delim = "\t", skip = 57) %>% 
  #Shortening column names
  janitor::clean_names() %>% 
  #Changing NA values under individuals to zero (0)
  mutate(ind_no_number_size_group = replace_na(ind_no_number_size_group, 0)) 
  
bester_all <- bester_outflow %>% 
  bind_rows(bester_trough)

#Renaming columns to simple column names
names(bester_all) <- col_names_simple 

#Extracting list of unique transects
transects_unique <- bester_all %>% 
  distinct(transect_code, segment_label, area_name)

# Extracting crabeater data -----------------------------------------------
#Extracting ALL crabeater observations in dataset
lobodon_all <- bester_all %>% 
  filter(str_detect(species, "Lobodon"))

# Inconsistencies ---------------------------------------------------------
#Checking observations with inconsistencies reported in species ID
bester_all_incon <- bester_all %>% 
  filter(str_detect(comment_inconsistencies, ".*spp.*"))

#Getting unique individual ID for discrepancies
ind_id <- bester_all_incon %>% 
  distinct(comment_individual) %>% 
  pull()

#Splitting between conflicting species and crabeaters/unrecognised
lob_disc <- data.frame()
lob_unr <- data.frame()
for(id in ind_id){
  comp <- bester_all_incon %>% filter(comment_individual == id)
  lob <- sum(str_detect(comp$species, "Lobodon"))
  un <- sum(str_detect(comp$species, "Unidentified"))
  if(lob == 0){
    print(paste(id, "no crabeaters"))
    next}
  if(lob == 1 & un == 1){
    comp <- comp %>% 
      #Ensuring the lowest number of individuals is kept
      mutate(number_ind = min(number_ind)) %>% 
      #Keep only crabeaters
      filter(species != "Unidentified")
    lob_unr <- rbind(lob_unr, comp)
    print(paste(id, "crabeater"))}
  else if(lob == 1 & un == 0){
    print(paste(id, "discrepancy"))
    lob_disc <- rbind(lob_disc, comp)}}
rm(comp, id, lob, un)

#Duplicate detections, but no inconsistencies in ID
lob_dup <- lobodon_all %>% 
  filter(!is.na(comment_dup)) %>%
  filter(is.na(comment_inconsistencies))

#Duplicate detections, but different bins
lob_bin_dif <- lobodon_all %>% 
  filter(str_detect(comment_inconsistencies, ".*bin.*") & !str_detect(comment_inconsistencies, ".*spp.*")) %>% 
  group_by(comment_individual) %>% 
  #Keeping minimum number of individuals identified per observation
  mutate(number_ind = case_when(str_detect(comment_inconsistencies, ".*size.*") ~ min(number_ind),
                                T ~ number_ind)) %>% 
  filter(!is.na(comment_dup))

#Duplicate detections, but differences in number of individuals
lob_num_dif <- lobodon_all %>% 
  filter(comment_inconsistencies == "size") %>% 
  group_by(comment_individual) %>% 
  #Keeping minimum number of individuals identified per observation
  mutate(number_ind = min(number_ind)) %>% 
  filter(!is.na(comment_dup))
  
# Merging clean data ------------------------------------------------------
#Crabeaters - no duplicates or inconsistencies
lob_confirmed <- lobodon_all %>% 
  filter(is.na(comment_dup)) %>%
  filter(is.na(comment_inconsistencies))

#Adding duplicates where mismatch ID was Unidentified
lob_confirmed <- lob_confirmed %>% 
  full_join(lob_unr, keep = F) %>% 
  full_join(lob_dup, keep = F) %>% 
  full_join(lob_bin_dif, keep = F) %>% 
  full_join(lob_num_dif, keep = F) %>% 
  distinct(comment_individual, .keep_all = T)

#Final check - Ensuring no individuals with conflicted species recorded were kept
sum(lob_disc$comment_individual %in% lob_confirmed$comment_individual)


# Final data cleaning -----------------------------------------------------
lob_confirmed <- lob_confirmed %>% 
  #Ensuring there are no observations without coordinates
  drop_na(lat_obs, lon_obs) %>% 
  #Removing observations from last bin starting at 345 m from the helicopter up to the horizon
  #because reliability of species ID diminishes with distance from helicopter
  filter(!str_detect(distance_interval, "Bin 6.*horizon.*"))

#Saving clean dataset
write_csv(lob_confirmed, "Cleaned_Data/FIL_2014_Filchner_Outflow_Trough_seal_census_cleaned.csv")


lob_confirmed %>% 
  group_by(transect_code, transect_length_km) %>% 
  summarise(tot_ind = sum(number_ind)) %>% 
  mutate(trans_width_km = ((345.73*2)/1000),
         trans_area_km = transect_length_km*trans_width_km,
         density_ind_km = tot_ind / trans_area_km)
