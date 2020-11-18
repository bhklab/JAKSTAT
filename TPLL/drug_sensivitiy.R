library(dplyr)
library(stringr)

source(file='functions_drug_screen.R')

##### 1. Parse the drug screen data into raw.sensitivity array #####
raw.sensitivity <- get_raw_sensitivity('./Data/drug_screen/', c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5', 'Dose6', 'Dose7'))

##### 2. Create sensitivity_info #####
sensitivity_info <- get_sensitivity_info(raw.sensitivity)

##### Create sensitivity_profile #####
sensitivity_profile <- data.frame(matrix(data=NA, ncol=7, nrow=0))
colnames(sensitivity_profile) <- c('aac_recomputed', 'ic50_recomputed', 'ic50_published', 'meanviability_published', 'HS', 'E_inf', 'EC50')

