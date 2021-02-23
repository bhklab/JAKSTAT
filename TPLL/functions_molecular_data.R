### functions used for processing TPLL molecular data ###

###load packages##
# library(SummarizedExperiment)
# library(dplyr)
# library(stringr)

##### Functions #####
read_and_format_data <- function(path){
  # Read data
  data <- read.csv(path)
  
  # unify the delimiter of the samples to '.' if it is prefixed with 'TP'
  for(col in colnames(data)){
    if(length(grep('^TP\\d+', col, value=TRUE))){
      colnames(data)[colnames(data) == col] <- gsub('_', '.', col)
    }
  }
  
  # Format all TP sample names to 'TP', followed by three digits 
  for(col in colnames(data)){
    # catch patient smaples (prefixed with TP)
    if(length(grep('^TP\\d+', col, value=TRUE))){
      # remove suffixes appended after '.' or '_' which indicates time series
      id <- strsplit(col, "[.|_]")[[1]][1]
      if(!length(grep('^TP\\d{3}$', id, value=TRUE))){
        if(length(substring(id, 3)) < 3){
          colnames(data)[colnames(data) == col] <- sub('(^TP)', '\\10', col)
        }
      }
    }
  }
  return(data)
}

get_avg_gep <- function(gep){
  ensembl_id <- unique(gep[,3])
  
  # get gene expression sample names
  gep_samples <- unique(colnames(gep))
  gep_samples <- gep_samples[- c(1,2,3)] # remove non-sample columns
  
  # average out the gene expression values of duplicated ensembl id.
  gep_avg <- data.frame(matrix(data=NA, ncol=ncol(gep[-c(1,2,3)]), nrow=length(ensembl_id)))
  colnames(gep_avg) <- gep_samples 
  rownames(gep_avg) <- ensembl_id
  gep_sample_len <- length(gep_samples)
  for(id in ensembl_id){
    found <- gep %>% filter(Ensembl_ID == id)
    means <- vector('numeric', gep_sample_len)
    for(i in 1:gep_sample_len){
      means[i] <- mean(as.numeric(found %>% pull(gep_samples[i])))
    }
    gep_avg[id, ] <- means
  }
  return(gep_avg)
}

merge_mut_dataframes <- function(df1, df2){
  merged <- df1
  
  ## identify genes that only exists in df2
  difference <- rownames(df2)[!rownames(df2) %in% rownames(df1)]
  
  ## add the the difference as new rows to merged
  merged[difference, ] <- NA
  
  ## add columns in df2 to merged
  merged[, colnames(df2)] <- NA
  
  ## populate the data for new columns
  for(col in colnames(df2)){
    for(gene in rownames(df2)){
      merged[gene, col] <- df2[gene, col]
    }
  }
  return(merged)
}

add_p_numbers_to_clinical_data <- function(clinical_data, p_number_dir){
  col_names <- c(
    'p_number_gep', 
    'p_number_miRNA', 
    'p_number_mRNASeq', 
    'p_number_snp', 
    'p_number_wes_paired',
    'p_number_wes_single',
    'p_number_wes_single_followup',
    'p_number_wgs'
  )
  clinical_data[, col_names] <- NA
  clin_rows <- rownames(clinical_data)
  
  for(col_name in col_names){
    annotation <- read.csv(paste0(p_number_dir, col_name, '.csv'))
    row_names <- annotation$TP_number
    row_names <- str_replace(row_names, '^TP\\d{3}_', format_clinical_rowname)
    rownames(annotation) <- row_names
    annotation_rows <- rownames(annotation)
    for(clin_row_name in clin_rows){
      if(clin_row_name %in% annotation_rows){
        clinical_data[clin_row_name, col_name] <- if(annotation[clin_row_name, 'p_number'] != '') annotation[clin_row_name, 'p_number'] else NA
      }
    }
  }
  
  return(clinical_data)
}

assign_p_number <- function(mol_data, datatype) {
  annotation <- NULL
  if(datatype == 'wes_single'){
    # merge wes and wes_followup annotations
    annotation <- read.csv(paste0('./Data/p-numbers/p_number_', datatype, '.csv'))
    followup <- read.csv(paste0('./Data/p-numbers/p_number_', datatype, '_followup', '.csv'))
    
    # format time series TP number to "TP###.t#
    followup$TP_number <- str_replace(followup$TP_number, '^TP\\d{3}_', format_clinical_rowname)
    
    annotation <- bind_rows(annotation, followup)
  }else{
    annotation <- read.csv(paste0('./Data/p-numbers/p_number_', datatype, '.csv'))
  }
  
  # if no p-number is assigned, assign the TP number. This is the case for control samples.
  for(i in rownames(annotation)){
    if(annotation[i, 'p_number'] == ''){
      annotation[i, 'p_number'] = annotation[i, 'TP_number']
    }
  }
  
  # find any TP numbers that does not exist in the annotation dataframe, and add them
  difference <- setdiff(colnames(mol_data), annotation$TP_number)
  if(length(difference) > 0){
    for(n in difference){
      annotation <- rbind(annotation, c(n, n))
    }
  }
  
  # note
  # 1. wes_mut_single_merged: TP029 found in the dataframe, but not in the annotation.
  # 2. wes_mut_single_merged: TP079 found in the annotation, but not in the dataframe.
  # 3. snp: timeseries data found in the dataframe but not in the annotation.
  # 4. snp: TP092, TP093, TP094, TP095 and TP096 found in the annotation but not in th dataframe.
  
  names(mol_data) <- annotation$p_number[match(names(mol_data), annotation$TP_number)]

  return(mol_data)
}

get_clinical_sample_data <- function(total_samples, filepath, timeseries_filepath){
  clinical_samples <- read.csv(filepath)
  row.names(clinical_samples) <- clinical_samples$p_number
  colnames(clinical_samples)[2] <- "cellid"
  clinical_samples <- add_column(clinical_samples, time_series = NA, .after="cellid")
  
  # add time series designation
  followup <- read.csv(timeseries_filepath)
  followup$TP_number <- str_replace(followup$TP_number, '^TP\\d{3}_', format_clinical_rowname)
  rownames(followup) <- followup$p_number
  for(pnum in followup$p_number){
    clinical_samples[pnum, 'time_series'] <- str_split(followup[pnum, 'TP_number'], '\\.')[[1]][2]
  }
  
  # add missing sample rows (control samples) to clinical_samples to ensure that all the samples are curated.
  additional <- setdiff(total_samples, rownames(clinical_samples))
  
  df <- data.frame(matrix(data=NA, ncol=ncol(clinical_samples), nrow=length(additional)))
  colnames(df) <- colnames(clinical_samples)
  rownames(df) <- additional
  
  for(n in additional){
    df[n, 'p_number'] <- n
    if(length(grep('^TP\\d{3}.t\\d{1}', n)) > 0){
      split <- str_split(n, "\\.")
      df[n, 'cellid'] <- split[[1]][1]
      df[n, 'time_series'] <- split[[1]][2]
    }else{
      df[n, 'cellid'] <- n
    }
  }
  
  clinical_samples <- rbind(clinical_samples, df)
  clinical_samples <- clinical_samples[order(rownames(clinical_samples)), ]
  
  return(clinical_samples)
}

format_clinical_rowname <- function(match){
  return(str_replace(match, '_', '.'))
}

get_gene_info <- function(genes, data){
  return(genes[which(genes$gene_id %in% rownames(data)),])
}

get_gene_names <- function(genes, data){
  gene_names <- genes[which(genes$gene_name %in% rownames(data)),] 
  gene_names <- gene_names[c('gene_id', 'gene_name')]
  rownames(gene_names) <- NULL
  return(gene_names %>% group_by(gene_name) %>% summarize(Ensembl_ID=sapply(list(gene_id), paste, collapse=","))) 
}

get_summarized_experiment <- function(data, clinical_sample_data, genes, gene_column_name, annotation) {
  filtered_clinical_sample_data <- clinical_sample_data[row.names(clinical_sample_data) %in% colnames(data), ]
  
  if(gene_column_name == 'gene_id'){
    filtered_gene_info <- get_gene_info(genes, data)
  }else{
    filtered_gene_info <- get_gene_names(genes, data)
  }
  
  # remove any genes that are not in both data, and filtered_gene_info
  filtered_data <- data[rownames(data) %in% filtered_gene_info[[gene_column_name]], ]
  
  # add everything in SummarizedExperiment
  Sum_Exp <- SummarizedExperiment(
    assays=list(exprs=as.matrix(filtered_data)), #gene expression data in matrix
    rowData=filtered_gene_info, #gene features Data
    colData=filtered_clinical_sample_data) #patient/sample info
  
  Sum_Exp@metadata$annotation <- annotation
  return(Sum_Exp)
}