library(stringr)
library(dplyr)
library(readxl)

get_sample_raw_sensitivity <- function(sample_name, corrected_name, datatype, path){
  sample_data <- read_excel(path, sheet=datatype)
  sample_data <- data.frame(sample_data[, c('DRUG_NAME', 'Max.Conc.tested', 'Min.Conc.tested', 'D1', 'D2', 'D3', 'D4', 'D5')])
  row_names <- paste(corrected_name, sample_data$DRUG_NAME, sep='_')
  doses <- c('D1', 'D2', 'D3', 'D4', 'D5')
  
  dose_df <- sample_data
  num <- 0
  for(d in doses){
    dose_df[, d] <- (dose_df$Min.Conc.tested * (10^num))/1000
    num <- num + 1
  }
  dose_df <- dose_df[, doses]
  rownames(dose_df) <- row_names
  
  viability_df <- sample_data
  viability_df <- viability_df[, doses]
  rownames(viability_df) <- row_names
  return(list('dose'=dose_df, 'viability'=viability_df))
}

root_dir <- './Data/additional_drug_screen_andersson_et_al/'

summary_path <- paste0(root_dir, 'DSS_PLL.xlsx')
summary_sheet <- read_excel(summary_path)
samples <- colnames(summary_sheet)[2:length(colnames(summary_sheet))]

clinical_sample <- readRDS('./Data/curated_data/clinical_sample_data.rds')
p_number_drug_sensitivity <- read.csv('./Data/p_number_drug_sensitivity.csv')

clinical_sample <- clinical_sample[, c('p_number', 'cellid')]
p_number_drug_sensitivity <- p_number_drug_sensitivity[, c('p_number', 'tp_number')]
pnumbers <- c(p_number_drug_sensitivity$p_number, clinical_sample$p_number)
tp_numbers <- c(p_number_drug_sensitivity$tp_number, clinical_sample$cellid)
sample_map <- data.frame(matrix(data=NA, ncol=0, nrow=length(pnumbers)))
sample_map$p_number <- pnumbers
sample_map$tp_number <- tp_numbers
sample_map <- sample_map[!duplicated(sample_map[1:2]), ]

corrected_names <- lapply(samples, function(name){
  print(name)
  if(str_detect(name, '^TP(\\d+)')){
    return(paste0('PSLRU', name))
  }
  if(str_detect(name, '^FM.')){
    return(str_replace_all(name, '\\.', ''))
  }
  if(str_detect(name, '^p(\\d+)')){
    pnum <- sample_map[sample_map$p_number == name, ]
    if(length(pnum$tp_number) > 0){
      return(pnum$tp_number[1])
    }else{
      return(name)
    }
  }
  if(str_detect(name, '^TPLL-')){
    return(str_replace_all(name, '-', ''))
  }
  return(name)
})

samples_df <- data.frame(matrix(data=NA, ncol=0, nrow=length(samples)))
samples_df$samples <- samples
samples_df$corrected <- unlist(corrected_names)

doses <- c('D1', 'D2', 'D3', 'D4', 'D5')
dose_df <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
viability_df <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
colnames(dose_df) <- doses 
colnames(viability_df) <- doses

datatype <- 'EC50'

files <- list.files(paste0(root_dir, "processed_data"))
sample_replace <- c("P0436", "P0619", "P0750")
for(sample in samples){
  # print(sample)
  filename <- ''
  if(sample %in% sample_replace){
    filename <- grep(str_replace(sample, '0', 'O'), files, ignore.case=TRUE, value=T)[1]
  }else{
    filename <- grep(sample, files, ignore.case=TRUE, value=T)[1]
  }
  if(!is.na(filename)){
    path <- paste0(root_dir, 'processed_data/', filename)
    corrected_name <- samples_df[samples_df$samples == sample, ]$corrected[1]
    sample_raw_sensitivity <- get_sample_raw_sensitivity(sample, corrected_name, datatype, path)
    dose_df <- rbind(dose_df, sample_raw_sensitivity[['dose']])
    viability_df <- rbind(viability_df, sample_raw_sensitivity[['viability']])
  }else{
    print(paste('sample not found:', sample))
  }
}

colnames(dose_df) <- c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5') 
colnames(viability_df) <- c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5') 

# Merge dose and viability data frames into one dataframe
raw.sensitivity <- cbind(dose_df, viability_df)

saveRDS(raw.sensitivity, './Data/curated_data/new_raw_sensitivity.rds')

# # order data by row names
# raw.sensitivity <- raw.sensitivity[order(rownames(raw.sensitivity)),]
# 
# conc_tested <- 5
# 
# # Make the merged dataframe into an array
# raw.sensitivity <- array(
#   c(
#     as.matrix(raw.sensitivity[ ,1:conc_tested]), 
#     as.matrix(raw.sensitivity[ ,(conc_tested + 1):(2 * conc_tested)])
#   ), 
#   c(nrow(raw.sensitivity), conc_tested , 2), 
#   dimnames=list(
#     rownames(raw.sensitivity), 
#     colnames(raw.sensitivity[ ,1:conc_tested]), 
#     c("Dose", "Viability")
#   )
# )

# ##### 2. Create sensitivity_info #####
# sensitivity_info <- get_sensitivity_info(raw.sensitivity, 'Dose1', 'Dose5')
# 
# ##### 3. Create sensitivity_profile #####
# sensitivity_profile <- get_sensitivity_profile(raw.sensitivity)
# 
# dose_rownames <- rownames(dose_df)
# viability_rownames <- rownames(viability_df)
# info_rownames <- rownames(sensitivity_info)
# profile_rownames <- rownames(sensitivity_profile)
# 
# setdiff(info_rownames, profile_rownames)
# 
# saveRDS(raw.sensitivity, './Data/curated_data/new_raw_sensitivity.rds')
# saveRDS(sensitivity_info, './Data/curated_data/new_sensitivity_info.rds')
# saveRDS(sensitivity_profile, './Data/curated_data/new_sensitivity_profile.rds')

# [1] "sample not found: FM.1158"                                                                                                                                                       
# [1] "sample not found: FM.1241"
# [1] "sample not found: FM.615"
# [1] "sample not found: FM.616"
