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