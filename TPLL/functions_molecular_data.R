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

order_dataframes <- function(df){
  return(df[,order(colnames(df))])
}

read_and_format_clinical_data <- function(path, all_samples, data_colnames, col_names){
  clinical_data <- read.csv(path) 
  rownames(clinical_data) <- clinical_data[,1]
  clinical_data <- clinical_data[-c(1)]
  clinical_data <- data.frame(t(clinical_data)) # transpose clinical_data patient id as rows and attirbutes as columns.
  
  # add any missing samples (patients and healthy individuals) into clinical_data and populate the columns with NA
  # this is done to ensure that all the samples in each SummarizedExperiment are included and not ommitted.
  clinical_missing <- all_samples[!all_samples %in% rownames(clinical_data)]
  clinical_data[clinical_missing,] <- NA
  
  # add cellid column, and populate it with rownames for SummarizedExeriment object to be created.
  clinical_data[, 'cellid'] <- rownames(clinical_data)
  
  # add time series columns and populate the cell with timeseries sample names if they exist
  clinical_data[,col_names] <- NA
  
  timeseries_samples <- unique(rownames(clinical_data)[grepl('^TP\\d{3}.', rownames(clinical_data))])
  for(sample in timeseries_samples){
    for(i in 1:length(data_colnames)){
      timeseries <- grep(sample, data_colnames[[i]], value=TRUE)
      if(length(timeseries)){
        clinical_data[sample, col_names[i]] <- substring(timeseries, first=1, last=5)
      }
    }
  }
  
  # make sure that the clinical_data is ordered alphabetically
  clinical_data <- clinical_data[order(rownames(clinical_data)), ]
  return(clinical_data)
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

get_filtered_clinical_data <- function(clinical_data, data){
  return(clinical_data[row.names(clinical_data) %in% colnames(data), ])
}

get_summarized_experiment <- function(data, clinical_data, genes, gene_column_name, annotation) {
  filtered_clinical_data <- get_filtered_clinical_data(clinical_data, data)
  
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
    colData=filtered_clinical_data) #patient/sample info
  
  Sum_Exp@metadata$annotation <- annotation
  return(Sum_Exp)
}