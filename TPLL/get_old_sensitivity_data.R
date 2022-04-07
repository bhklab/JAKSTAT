library(dplyr)
library(stringr)

# replaces drug sensitivity sample id with a TP number based on id annotation file
replace_sample_id <- function(sensitivity_data, annotation_file_path, is_sensitivity_info = FALSE){
  annotation <- read.csv(annotation_file_path)
  annotation <- data.frame(annotation)
  annotation <- annotation[annotation$helsinki_drug_screen_id != "", ]
  rownames(annotation) <- annotation$helsinki_drug_screen_id
  
  data_rownames <- rownames(sensitivity_data)
  
  new_cellid <- vector()
  new_rownames <- vector()
  
  for(rowname in data_rownames){
    id_substr <- sub("\\_.*", "", rowname)
    if(id_substr %in% rownames(annotation)){
      new_name <- gsub(id_substr, annotation[id_substr, 'tp_number'], rowname)
      new_rownames <- c(new_rownames, new_name)
      if(is_sensitivity_info){
        new_cellid <- c(new_cellid, annotation[id_substr, 'tp_number'])
      }
    }else{
      new_rownames <- c(new_rownames, rowname)
      if(is_sensitivity_info){
        new_cellid <- c(new_cellid, sensitivity_data[rowname, 'cellid'])
      }
    }
  }
  
  rownames(sensitivity_data) <- new_rownames
  if(is_sensitivity_info){
    sensitivity_data$cellid <- new_cellid
  }
  
  return(sensitivity_data)
}

root_dir <- './Data/drug_screen/'
colNames <- c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5', 'Dose6', 'Dose7')
p_number_file <- './Data/p_number_drug_sensitivity.csv' 
corrected_data_dir <- './Data/drug_screen_correction/'

# read in all the drug screen data, and merge all dose and viability data
files <- list.files(path=root_dir, pattern = "*.Rdata$", recursive=TRUE)
corrected_files <- list.files(path=corrected_data_dir, pattern = "*.Rdata$")
corrected <- str_split(corrected_files, "_")
corrected <- unlist(lapply(corrected, function(item){return(item[1])}))

conc_tested <- as.numeric(length(colNames)) # number of concentrations tested for each experiment (same as number of columns)

dose_merged <- data.frame(matrix(data=NA, ncol=conc_tested, nrow=0))
colnames(dose_merged) <- colNames

viability_merged <- data.frame(matrix(data=NA, ncol=conc_tested, nrow=0))
colnames(viability_merged) <- colNames

for(file in files){
  load(paste0(root_dir, file))
  sample_name <- unlist(str_split(rownames(dose.df)[1], "_"))[1]
  
  if(sample_name %in% corrected){
    print(paste(sample_name, "correction"))
    rm(dose.df, viability.df)
    load(paste0(corrected_data_dir, corrected_files[str_detect(corrected, sample_name)]))
  }    
  
  dose_merged <- rbind(dose_merged, dose.df)
  viability_merged <- rbind(viability_merged, viability.df)
  rm(dose.df, viability.df)
}

# Remove Dose2 and Dose 6 columns to make total number of columns to be 5 to merge with the new data.
dose_merged <- dose_merged[, !names(dose_merged) %in% c('Dose2', 'Dose6')]
colnames(dose_merged) <-  c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5')
viability_merged <- viability_merged[, !names(viability_merged) %in% c('Dose2', 'Dose6')]
colnames(viability_merged) <-  c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5')

# Merge dose and viability data frames into one dataframe
raw.sensitivity <- cbind(dose_merged, viability_merged)

# replace p number as a sample id with TP number
raw.sensitivity <- replace_sample_id(raw.sensitivity, p_number_file)

saveRDS(raw.sensitivity, './Data/curated_data/old_raw_sensitivity.rds')
