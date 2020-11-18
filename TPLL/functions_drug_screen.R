### functions used for processing TPLL drug screen data ###

get_raw_sensitivity <- function(root_dir, colNames){
  # read in all the drug screen data, and merge all dose and viability data
  files <- list.files(path=root_dir, pattern = "*.Rdata$", recursive=TRUE)
  
  conc_tested <- as.numeric(length(colNames)) # number of concentrations tested for each experiment (same as number of columns)
  
  dose_merged <- data.frame(matrix(data=NA, ncol=conc_tested, nrow=0))
  colnames(dose_merged) <- colNames
  
  viability_merged <- data.frame(matrix(data=NA, ncol=conc_tested, nrow=0))
  colnames(viability_merged) <- colNames
  
  for(file in files){
    load(paste0(root_dir, file))
    dose_merged <- rbind(dose_merged, dose.df)
    viability_merged <- rbind(viability_merged, viability.df)
    rm(dose.df, viability.df)
  }
  
  # Merge dose and viability data frames into one dataframe
  raw.sensitivity <- cbind(dose_merged, viability_merged)
  
  # Make the merged dataframe into an array
  raw.sensitivity <- array(
    c(
      as.matrix(raw.sensitivity[ ,1:conc_tested]), 
      as.matrix(raw.sensitivity[ ,(conc_tested + 1):(2 * conc_tested)])
    ), 
    c(nrow(raw.sensitivity), conc_tested , 2), 
    dimnames=list(
      rownames(raw.sensitivity), 
      colnames(raw.sensitivity[ ,1:conc_tested]), 
      c("Dose", "Viability")
    )
  )
  return(raw.sensitivity)
}

get_sensitivity_info <- function(raw.sensitivity) {
  sensitivity_info <- data.frame(matrix(data=NA, ncol=5, nrow=0))
  colnames(sensitivity_info) <- c('cellid', 'drugid', 'chosen.min.range', 'chosen.max.range', 'rm.by.conc.range')
  doseDF <- raw.sensitivity[,,'Dose']
  for(row in rownames(doseDF)){
    split <- strsplit(row, '_')
    cell <- if(stringr::str_detect(split[[1]][1], '^CD*')) paste(split[[1]][1], split[[1]][2], sep='_') else split[[1]][1]
    drug <- if(stringr::str_detect(split[[1]][1], '^CD*')) split[[1]][3] else split[[1]][2]
    sensitivity_info[row, ] <- c(cell, drug, doseDF[row, 'Dose1'], doseDF[row, 'Dose7'], FALSE)
  }
  return(sensitivity_info)
}