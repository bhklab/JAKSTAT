library(PharmacoGx)
library(dplyr)
library(stringr)

source(file='functions_drug_screen.R')

##### 1. Parse the drug screen data into raw.sensitivity array #####
raw.sensitivity <- get_raw_sensitivity('./Data/drug_screen/', c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5', 'Dose6', 'Dose7'))

##### 2. Create sensitivity_info #####
sensitivity_info <- get_sensitivity_info(raw.sensitivity)

##### 3. Create sensitivity_profile #####
sensitivity_profile <- get_sensitivity_profile(raw.sensitivity)