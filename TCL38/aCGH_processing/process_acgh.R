library(stringr)
library(dplyr)

get_process_info <- function(cellname, filename, lines){
  process_lines <- lines[1: (data_start - 1)]
  process_lines <- process_lines[!str_detect(process_lines, "^\\t+$")]
  process_lines <- str_remove_all(process_lines, "\\t")
  process_info <- lapply(process_lines, function(line){
    unlist(str_split(line, ": "))
  })
  info_data <- lapply(process_info, function(info){
    return(str_trim(info[2], "both"))
  })
  info_col <- unlist(lapply(process_info, function(info){
    return(info[1])
  }))
  info_col <- str_remove_all(info_col, "[^\\s\\w_]")
  info_col <- str_to_lower(info_col)
  info_col <- str_replace_all(info_col, "\\s", "_")
  names(info_data) <- info_col
  info_data <- c(list(cell=cellname, file=filename), info_data)
  info_data <- info_data[names(info_data) != ""]
  return(info_data)
}

process_info_df <- data.frame(matrix(nrow=0, ncol=0))

# Data obtained from /repository/vmuvienna_neubauer/procdata/20211110/Tcell tumors_aCGH_Aberr Reports
files <- list.files("./data")

# combined = data.frame(matrix(nrow=0, ncol=0))
df_list <- list()
for(file in files){
  lines <- readLines(paste0("./data/", file), warn = FALSE)
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
  # write.csv(data_df, file=paste0("./data/formatted/", cellname, ".csv"), row.names=FALSE)
  data_df$Gene_Names <- str_remove_all(data_df$Gene_Names, '"')
  data_df$Hs_hg19_CNV_20120403 <- str_remove_all(data_df$Hs_hg19_CNV_20120403, '"')
  data_df$Hs_hg19_miRNA_20120403 <- str_remove_all(data_df$Hs_hg19_miRNA_20120403, '"')
  data_df$Mm_mm8_miRNA_20090908 <- str_remove_all(data_df$Mm_mm8_miRNA_20090908, '"')
  df_list[[cellname]] <- data_df
}

# output data
# The following data are used to create MultiAssayExperiment object.
saveRDS(df_list, "../create_multiassay/data/acgh_assay_data.rds")
write.csv(process_info_df, "../create_multiassay/data/acgh_samples.csv")
