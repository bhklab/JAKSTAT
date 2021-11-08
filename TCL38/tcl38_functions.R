library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
library(rlist)

#'Returns a dataframe that includes all the cellline-drug combinations tested in each PSet given as a list.
#'@param pset.list `list` is a list of PSets, presumably subsetted by a list of cell lines (TCL38 cell lines in this case).
#'@return a dataframe containing all combinations of cellid and drugid, and whether a specific cellline-drug combination is present in each PSet
get_cell_drug_combinations <- function (pset.list) {
  cells <- c()
  drugs <- c()
  for(dataset.name in names(pset.list)){
    cells <- c(cells, pset.list[[dataset.name]]@cell$cellid)
    drugs <- c(drugs, pset.list[[dataset.name]]@drug$drugid)
  }
  cells <- unique(cells)
  drugs <- unique(drugs)
  cell_drug_combinations <- expand.grid(cellid=cells, drugid=drugs)
  for(dataset.name in names(pset.list)){
    info <- pset.list[[dataset.name]]@sensitivity[["info"]]
    cell_drug_combinations[dataset.name] <- ifelse(
      is.na(match(
        paste0(cell_drug_combinations$cellid, cell_drug_combinations$drugid),
        paste0(info$cellid, info$drugid),
      )),
      NA, "Yes"
    )
  }
  return(cell_drug_combinations)
}

#'Subsets drug sensitivity data in each PSet in a given list of PSets by cellid and drugid, then merges all the sensitivity data.
#'@param pset.list `list` is a list of PSets, presumably subsetted by a list of cell lines (TCL38 cell lines in this case).
#'@param cellid `character` is a cellline id if interest. Obtained from the resulting dataframe of the get_cell_drug_combinations function.
#'@param drugid `character` is a drug id of interest. Obtained from the resulting dataframe of the get_cell_drug_combinations function.
#'@return A list that contains the following dataframes: 
#'info (Merged drug sensitivity experiment info for the given combination of cellid and drug from all the PSets in pset.list), 
#'profiles (Merged drug sensitivity experiment profile for the given combination of cellid and drug from all the PSets in pset.list),
#'raw (Merged raw sensitivity data for the given combination of cellid and drug from all the PSets in pset.list)
subset_sensitivity_data <- function(pset.list, cellid, drugid){
  sens_list <- list()
  for(dataset.name in names(pset.list)){
    info <- pset.list[[dataset.name]]@sensitivity[["info"]]
    info <- info[info$cellid == cellid & info$drugid == drugid, ]
    profiles <- pset.list[[dataset.name]]@sensitivity[["profiles"]]
    profiles <- profiles[rownames(profiles) %in% rownames(info), ]
    raw <- data.frame(pset.list[[dataset.name]]@sensitivity[["raw"]])
    raw <- raw[rownames(raw) %in% rownames(info), ]
    sens_list[[dataset.name]] <- list(info=info, profiles=profiles, raw=raw)
  }
  
  sensitivity_info_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
  sensitivity_profile_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
  
  for(dataset.name in names(sens_list)){
    sensitivity_info_total <- merge_sensitivity_data(
      sens_list[[dataset.name]][["info"]],
      dataset.name,
      sensitivity_info_total
    )
    sensitivity_profile_total <- merge_sensitivity_data(
      sens_list[[dataset.name]][["profiles"]],
      dataset.name,
      sensitivity_profile_total
    )
  }
  
  # Standardize the dose/viability col names
  dose_colnames_total <- c()
  for(dataset.name in names(sens_list)){
    tmp <- sens_list[[dataset.name]][["raw"]]
    dose_colnames_total <- c(dose_colnames_total, colnames(tmp))
  }
  dose_colnames_total <- format_dose_col(dose_colnames_total)
  dose_colnames_total <- unique(dose_colnames_total)
  
  doses <- dose_colnames_total[grepl(".Dose", dose_colnames_total)]
  doses <- str_sort(doses, numeric=TRUE)
  viabilities <- dose_colnames_total[grepl(".Viability", dose_colnames_total)]
  viabilities <- str_sort(viabilities, numeric=TRUE)
  dose_colnames_total <- c(doses, viabilities)
  
  # Merge raw sensitivity data
  sensitivity_raw_total <- data.frame(matrix(data=NA, ncol=length(dose_colnames_total), nrow=0))
  colnames(sensitivity_raw_total) <- dose_colnames_total
  for(dataset.name in names(sens_list)){
    tmp <- sens_list[[dataset.name]][["raw"]]
    rownames(tmp) <- paste0(rownames(tmp), ".", dataset.name)
    colnames_tmp <- colnames(tmp)
    colnames(tmp) <- format_dose_col(colnames_tmp)
    tmp[setdiff(colnames(sensitivity_raw_total), colnames(tmp))] <- NA
    sensitivity_raw_total <- dplyr::bind_rows(sensitivity_raw_total, tmp)
  }
  
  subset_sensitivity <- list(info=sensitivity_info_total, profiles=sensitivity_profile_total, raw=sensitivity_raw_total)
  
  return(subset_sensitivity)
}


#'A helper function used to merge sensitivity info and sensitivity profiles data
merge_sensitivity_data <- function(sensitivity_data, dataset.name, sensitivity_total){
  total <- sensitivity_total
  tmp <- sensitivity_data
  rownames(tmp) <- paste0(rownames(tmp), ".", dataset.name)
  tmp["dataset"] <- dataset.name
  if(nrow(total) == 0){
    total <- tmp
  }else{
    col_total <- c(colnames(total), colnames(tmp))
    col_total <- unique(col_total)
    total[setdiff(col_total, colnames(total))] <- NA
    tmp[setdiff(col_total, colnames(tmp))] <- NA
    total <- rbind(total, tmp)
  }
  return(total)
}

#'A helper function to format dose/viability colums so that the columns are standardized for merging raw sensitivity data from different PSets.
format_dose_col <- function(dose_colnames){
  format <- function(colname) {
    dose_num <- regmatches(colname, gregexpr("[[:digit:]]+", colname))
    dose_type <- sub('.*\\.', '', colname)
    return(paste0("dose", dose_num, ".", dose_type))
  }
  return(unlist(lapply(dose_colnames, format)))
}