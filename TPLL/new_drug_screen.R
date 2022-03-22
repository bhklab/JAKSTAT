library(stringr)
library(dplyr)
library(readxl)

get_sample_raw_sensitivity <- function(sample_name, datatype, path){
  sample_data <- read_excel(path, sheet=datatype)
  sample_data <- data.frame(sample_data[, c('DRUG_NAME', 'Max.Conc.tested', 'Min.Conc.tested', 'D1', 'D2', 'D3', 'D4', 'D5')])
  row_names <- paste(sample_name, sample_data$DRUG_NAME, sep='_')
  doses <- c('D1', 'D2', 'D3', 'D4', 'D5')
  
  dose_df <- sample_data
  num <- 0
  for(d in doses){
    dose_df[, d] <- dose_df$Min.Conc.tested * (10^num)
    num <- num + 1
  }
  dose_df <- dose_df[, doses]
  rownames(dose_df) <- row_names
  
  viability_df <- sample_data
  viability_df <- viability_df[, doses]
  rownames(viability_df) <- row_names
  return(list('dose'=dose_df, 'viability'=viability_df))
}

summary_path <- './Data/drug_screen_andersson_et_al/DSS_PLL.xlsx'
summary_sheet <- read_excel(summary_path)
samples <- colnames(summary_sheet)[2:length(colnames(summary_sheet))]

doses <- c('D1', 'D2', 'D3', 'D4', 'D5')
ic50_dose <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
ic50 <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
colnames(ic50_dose) <- doses 
colnames(ic50) <- doses

ec50_dose <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
ec50 <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
colnames(ec50_dose) <- doses 
colnames(ec50) <- doses

datatype <- 'EC50'

files <- list.files("./Data/drug_screen_andersson_et_al/processed_data")
sample_replace <- c("P0436", "P0619", "P0750")
for(sample in samples){
  print(sample)
  filename <- ''
  if(sample %in% sample_replace){
    filename <- grep(str_replace(sample, '0', 'O'), files, ignore.case=TRUE, value=T)[1]
  }else{
    filename <- grep(sample, files, ignore.case=TRUE, value=T)[1]
  }
  if(!is.na(filename)){
    path <- paste0('./Data/drug_screen_andersson_et_al/processed_data/', filename)
    sample_raw_sensitivity <- get_sample_raw_sensitivity(sample, datatype, path)
    # ic50_dose <- rbind(ic50_dose, sample_raw_sensitivity[['dose']])
    # ic50 <- rbind(ic50, sample_raw_sensitivity[['viability']])
    ec50_dose <- rbind(ec50_dose, sample_raw_sensitivity[['dose']])
    ec50 <- rbind(ec50, sample_raw_sensitivity[['viability']])
  }
}

raw.sensitivity <- array(
  c(
    as.matrix(dose_df), 
    as.matrix(viability_df)
  ), 
  c(nrow(dose_df), 5, 2), 
  dimnames=list(
    rownames(dose_df), 
    doses, 
    c("Dose", "Viability")
  )
)