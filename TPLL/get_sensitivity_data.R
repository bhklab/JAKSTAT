library(dplyr)
library(stringr)
source(file='functions_drug_screen.R') # functions used to process drug sensitivity data

#######################################
##### PROCESS DRUG SCREENING DATA #####
#######################################

##### 1. Parse the drug screen data into raw.sensitivity array #####

# load('./Data/drug_screen_correction/p1332_response_data_corrected.Rdata')
# rownames(dose.df) <- str_replace(rownames(dose.df), "CD3_3", "p1332")
# save("dose.df", "viability.df", file='./Data/drug_screen_correction/p1332_response_data_corrected.Rdata')

raw.sensitivity <- get_raw_sensitivity(
  './Data/drug_screen/', 
  c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5', 'Dose6', 'Dose7'), 
  './Data/p_number_drug_sensitivity.csv',
  './Data/drug_screen_correction/'
)

##### 2. Create sensitivity_info #####
sensitivity_info <- get_sensitivity_info(raw.sensitivity)

##### 3. Create sensitivity_profile #####
sensitivity_profile <- get_sensitivity_profile(raw.sensitivity)

saveRDS(raw.sensitivity, './Data/curated_data/raw_sensitivity.rds')
saveRDS(sensitivity_info, './Data/curated_data/sensitivity_info.rds')
saveRDS(sensitivity_profile, './Data/curated_data/sensitivity_profile.rds')