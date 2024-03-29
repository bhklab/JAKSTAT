library(dplyr) 
library(tidyr)
library(stringr)

get_slot_data <- function(pset_subsets, slot_type, isCuration=FALSE){
  slot_colnames_list <- list()
  slot_data_list <- list()
  for(dataset.name in names(pset_subsets)){
    slot_data <- if(isCuration) pset_subsets[[dataset.name]]@curation[[slot_type]] else attr(pset_subsets[[dataset.name]], slot_type)
    if(isCuration & slot_type == "tissue"){
      slot_data$cellid <- rownames(slot_data)
    }
    slot_colnames_list[[dataset.name]] <- colnames(slot_data)
    slot_data_list[[dataset.name]] <- slot_data
  }
  common_colnames <- Reduce(intersect, slot_colnames_list)
  common_data <- lapply(slot_data_list, function(x) x %>% select(all_of(common_colnames)))
  slot_data_merged <- dplyr::bind_rows(common_data)
  slot_data_merged <- dplyr::distinct(slot_data_merged)
  
  for(dataset.name in names(pset_subsets)){
    tmp <- slot_data_list[[dataset.name]]
    tmp <- tmp[names(tmp)[!names(tmp) %in% common_colnames]]
    if(ncol(tmp) > 0){
      colnames(tmp) <- paste0(colnames(tmp), ".", dataset.name)
      slot_data_merged <- merge(slot_data_merged, tmp, by=0, all=TRUE)
      slot_data_merged$Row.names = NULL
      id_col_name <- if(isCuration) paste0("unique.", slot_type, "id") else paste0(slot_type, "id")
      id_col_name <- paste0(slot_type, "id")
      if(isCuration){
        id_col_name <- if(slot_type == "tissue") "cellid" else paste0("unique.", slot_type, "id")
      }
      rownames(slot_data_merged) <- slot_data_merged[[id_col_name]]
    }
  }
  if(isCuration & slot_type == "tissue"){
    slot_data_merged$cellid = NULL
  }
  return(slot_data_merged)
}

get_pset_name <- function(dataset, available.psets){
  if(str_detect(dataset, "GDSC")){
    gdsc_ver <- regmatches(dataset, gregexpr("[[:digit:]]+", dataset))
    available_versions <- available.psets$version[available.psets$`Dataset Name` == 'GDSC']
    version <- grep(paste0("v", gdsc_ver), available_versions, value=TRUE)
    return(paste0("GDSC_", version))
  }else{
    return(available.psets$`PSet Name`[available.psets$`Dataset Name` == dataset])
  }
}

merge_molecular_profile <- function(pset_subsets, mol_data_type){
  coldata_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
  rowdata_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
  for(dataset.name in names(PSet_Subsets)){
    mol_prof <- PSet_Subsets[[dataset.name]]@molecularProfiles[[mol_data_type]]
    if(!is.null(mol_prof)){
      print(dataset.name)
      # parse row data (genes)
      rowdata_tmp <- data.frame(mol_prof@elementMetadata)
      if(nrow(rowdata_total) == 0){
        rowdata_total <- rowdata_tmp
      }else{
        rowdata_total <- rbind(rowdata_total, rowdata_tmp)
        rowdata_total <- dplyr::distinct(rowdata_total)
      }
      
      # parse col data
      coldata_tmp <- data.frame(mol_prof@colData)
      if(nrow(coldata_tmp) > 0){
        rownames(coldata_tmp) <- paste0(rownames(coldata_tmp), ".", dataset.name)
        coldata_tmp["dataset"] <- dataset.name
        if(nrow(coldata_total) == 0){ 
          coldata_total <- coldata_tmp # initialize coldata_total df.
        }else{
          colnames_total <- c(colnames(coldata_total), colnames(coldata_tmp))
          colnames_total <- unique(colnames_total)
          coldata_total[setdiff(colnames_total, colnames(coldata_total))] <- NA
          coldata_tmp[setdiff(colnames_total, colnames(coldata_tmp))] <- NA
          coldata_total <- rbind(coldata_total, coldata_tmp)
        }
      }
    }
  }
  rownames(rowdata_total) <- rowdata_total$rownames
  rowdata_total <- rowdata_total[order(rownames(rowdata_total)), , drop=F]
  
  assay_total <- data.frame(matrix(data=NA, ncol=0, nrow=nrow(rowdata_total)))
  rownames(assay_total) <- rownames(rowdata_total)
  for(dataset.name in names(PSet_Subsets)){
    mol_prof <- PSet_Subsets[[dataset.name]]@molecularProfiles[[mol_data_type]]
    if(!is.null(mol_prof)){
      assay_tmp <- data.frame(mol_prof@assays@data@listData[["exprs"]])
      if(nrow(assay_tmp) > 0){
        colnames(assay_tmp) <- paste0(colnames(assay_tmp), ".", dataset.name)
        empty_df <- data.frame(matrix(data=NA, ncol=ncol(assay_tmp), nrow=length(setdiff(rownames(rowdata_total), rownames(assay_tmp)))))
        colnames(empty_df) <- colnames(assay_tmp)
        rownames(empty_df) <- setdiff(rownames(rowdata_total), rownames(assay_tmp))
        assay_tmp <- rbind(assay_tmp, empty_df)
        assay_tmp <- assay_tmp[order(rownames(assay_tmp)), ]
        assay_total <- cbind(assay_total, assay_tmp)
      }
    }
  }
  
  sum_exp <- SummarizedExperiment(
    assays=list(exprs=as.matrix(assay_total)),
    rowData=rowdata_total, 
    colData=coldata_total 
  ) 
  sum_exp@metadata$annotation <- mol_data_type
  return(sum_exp) 
}

merge_sensitivity_data <- function(sensitivity_data, dataset.name, sensitivity_total, datatype){
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

format_dose_col <- function(dose_colnames){
  format <- function(colname) {
    dose_num <- regmatches(colname, gregexpr("[[:digit:]]+", colname))
    dose_type <- sub('.*\\.', '', colname)
    return(paste0("dose", dose_num, ".", dose_type))
  }
  return(unlist(lapply(dose_colnames, format)))
}