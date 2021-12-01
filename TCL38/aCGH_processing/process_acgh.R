library(stringr)
library(dplyr)

source("./functions.R")

process_info_df <- data.frame(matrix(nrow=0, ncol=0))

files <- list.files("./data/original")

for(file in files){
  lines <- readLines(paste0("./data/original/", file), warn = FALSE)
  data_start <- str_which(lines, "AberrationNo")
  
  data_lines <- lines[data_start:(length(lines))] 
  
  # get cell name and number of calls - the last line of the read
  sample <- data_lines[length(data_lines)]
  cellname <- str_extract(sample, ".+(?==\\d+)")
  cellname <- str_trim(cellname, "both")
  cellname <- str_replace_all(cellname, "\\s", "_")
  num_calls <- as.numeric(str_extract(sample, "(?<==)\\d+"))
  
  # extract sample processing info, and populate the process_info_df
  process_info <- get_process_info(cellname, file, lines[1: (data_start - 1)])
  if(dim(process_info_df)[1] == 0){
    process_info_df <- data.frame(matrix(nrow=0, ncol=length(names(process_info))))
    colnames(process_info_df) <- names(process_info)
    process_info_df <- rbind(process_info_df, process_info)
  }else{
    missing_in_df <- setdiff(names(process_info), colnames(process_info_df))
    missing_in_new_row <- setdiff(colnames(process_info_df), names(process_info))
    
    if(length(missing_in_df) > 0){
      process_info_df[missing_in_df] <- NA
    }
    if(length(missing_in_new_row)){
      missing <- c(rep(NA, length(missing_in_new_row)))
      names(missing) <- missing_in_new_row
      process_info <- c(process_info, missing)
    }
    process_info_df <- rbind(process_info_df, process_info)
  }
  
  # read everything starting at the header line where the data begins.
  # extract header
  header <- str_split(data_lines[1], pattern = "\t")
  header <- unlist(header)
  header <- str_replace_all(header, "\\s", "_")
  header <- str_remove_all(header, "[^\\w_]")
  
  # extract data
  data <- data_lines[3:(num_calls + 2)]
  data_df <- data.frame(matrix(nrow=length(data), ncol=length(header)))
  colnames(data_df) <- header
  
  data_list <- lapply(data, function(line){
    tmp <- str_split(line, pattern="\t")
    tmp <- unlist(tmp)
    names(tmp) <- header
    return(tmp)
  })
  
  col_values <- lapply(header, function(col_name){
    col_item <- lapply(data_list, function(item){
      return(item[col_name])
    })
    col_item <- unlist(col_item)
    return(col_item)
  })
  names(col_values) <- header
  
  for(col_name in header){
    data_df[[col_name]] <- col_values[[col_name]]
  }
  write.csv(data_df, file=paste0("./data/formatted/", cellname, ".csv"), row.names=FALSE)
}
