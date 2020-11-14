###load packages##
library(PharmacoGx)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

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
  return(clinical_data)
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

##################################
##### Prep: Read in all data #####
##################################

#1. Read all data to be included in the PSet and format patient names.

# load Ensembl reference
load("./Data/Ensembl.v75.annotation.RData") #for hg17 (Ensembl v75)

gep <- read_and_format_data("./Data/GEP_matrix_69_TPLL_vs_10_healthy_CD3_pan_T.csv")
snp <- read_and_format_data("./Data/sCNA_of_82_TPLL_and_6_FU_samples.csv")
wgs_mut <- read_and_format_data("./Data/WGS_mutations_3_TPLL_tg_pairs.csv")
wes_mut_paired <- read_and_format_data("./Data/WES_mutations_17_TPLL_tg_pairs.csv")
wes_mut_single <- read_and_format_data("./Data/WES_mutations_32_TPLL_tumor_single.csv")
wes_mut_single_followup <- read_and_format_data("./Data/WES_mutations_5_patients_with_FU.csv")

#2. for gene expression profiles, average out the gene expression values of duplicated ensembl id, and set ensembl gene id as rown name.
gep <- gep[!is.na(gep$Ensembl_ID),] # remove rows where Ensembl_ID is NA
gep_avg <- get_avg_gep(gep)

#3. remove unnecessary columns and assign gene names as row names for SNP, WES and WGS data.

# SNP
rownames(snp) <- snp[,1] # change row names to gene names.
snp <- snp[-c(1,2)] # remove the gene names column and mean_CN column

# WES Mutation Paired
rownames(wes_mut_paired) <- wes_mut_paired$gene
wes_mut_paired <- wes_mut_paired[-c(1)]

# WES Mutaation Single
rownames(wes_mut_single) <- wes_mut_single$gene
wes_mut_single <- wes_mut_single[-c(1)]

# WES Mutation Single Followup
rownames(wes_mut_single_followup) <- wes_mut_single_followup$gene
wes_mut_single_followup <- wes_mut_single_followup[-c(1)]

# merge wes_mut_single and wes_mut_single_followup (timeseries data for wes_mut_single)
wes_mut_single_merged <- merge_mut_dataframes(wes_mut_single, wes_mut_single_followup)

# WGS Mutation 
rownames(wgs_mut) <- wgs_mut$gene
wgs_mut <- wgs_mut[-c(1)]


#############################################
##### 1: Create curationCell data.frame #####
#############################################

# NOTE: Each "cell" represents each smaple
# Concatenate all sample names from molecular profiles + drug sensitivity data into one vector
total_samples <- unique(c(colnames(gep_avg), colnames(snp), colnames(wes_mut_paired), colnames(wes_mut_single_merged), colnames(wgs_mut)))

# create data.frame with unique sample names
curationCell <- data.frame(matrix(ncol=1, nrow=(length(total_samples))))
colnames(curationCell) <- c("unique.cellid")
curationCell$unique.cellid <- total_samples
rownames(curationCell) <- curationCell$unique.cellid


#############################################
##### 2: Create curationTissue data.frame #####
#############################################
curationTissue <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) 
curationTissue$cellid <- curationCell$unique.cellid
curationTissue$unique.tissueid <- rep('Lymphoid', length(total_samples))
curationTissue$tissue_type <- ifelse(startsWith(curationTissue$cellid, 'TP'), 'T-PLL cells', 'CD3 Pan T cells')
rownames(curationTissue) <- curationCell$unique.cellid


#####################################
##### 3: Create Cell data.frame #####
#####################################

# rownames must be unique sample names. All samples names from molecular profiles + drug sensitivity data must be included here.
# this is the reference dataframe used to dictate which smaples are included in the PSet. Any data samples that are not in this data frame will be removed.
cell <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) #ncol depends on number of sample features to include.
cell$cellid <- curationCell$unique.cellid
cell$tissueid <- curationTissue$unique.tissueid
cell$tissue_type <- curationTissue$tissue_type
rownames(cell) <- curationCell$unique.cellid


################################################################
##### 4: Create summarizedExperiments (molecular profiles) #####
################################################################

##### 4.0.1 Prepare the clinical_data data frame to be used as colData in each SummarizedExperiment object #####
clinical_data <- read_and_format_clinical_data(
  path = "./Data/clinical_data_2020_10_23.csv", 
  all_samples = rownames(cell),
  data_colnames = list(colnames(snp), colnames(wes_mut_single_merged)), 
  col_names = c('timeseries_snp', 'timeseries_wes_mut_single'))

##### 4.1 Create Gene Expression Data SummarizedExperiment #####
TPLL_GEP <- get_summarized_experiment(gep_avg, clinical_data, features_gene, 'gene_id', 'rnaseq')

##### 4.2 Create SNP Data SummarizedExperiment #####
TPLL_SNP <- get_summarized_experiment(snp, clinical_data, features_gene, 'gene_name', 'mutation')

##### 4.3 Create WGS (Whole Genome Sequencing) mutation data SummarizedExperiment #####
TPLL_WGS_MUT <- get_summarized_experiment(wgs_mut, clinical_data, features_gene, 'gene_name', 'mutation')

##### 4.4 Create WES (Whole Exome Sequencing) single-end and paired-end mutation data SummarizedExperiments#####
### Single-end Mutation Data ###
# wes_mut_single_merged[wes_mut_single_merged == 'unmutated'] <- 'wt'
# tmp <- wes_mut_single_merged
# for(col in colnames(wes_mut_single_merged)){
#   print(col)
#   for(row in rownames(wes_mut_single_merged)){
#     if(!is.na(wes_mut_single_merged[row, col]) && wes_mut_single_merged[row, col] != 'wt'){
#       wes_mut_single_merged[row, col] <- str_match(wes_mut_single_merged[row, col], ":(p.*?):")[1,2]
#     }
#   }
# }
TPLL_WES_MUT_SINGLE <- get_summarized_experiment(wes_mut_single_merged, clinical_data, features_gene, 'gene_name', 'mutation')

### Paied-end Mutation Data ###
TPLL_WES_MUT_PAIRED <- get_summarized_experiment(wes_mut_paired, clinical_data, features_gene, 'gene_name', 'mutation')

##########################
##### 4: Create PSet #####
##########################

PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=list(
    "rnaseq"=TPLL_GEP, 
    "snp"=TPLL_SNP,
    "mut_wes_single"=TPLL_WES_MUT_SINGLE, 
    "mut_wes_paired"=TPLL_WES_MUT_PAIRED,
    "mut_wgs"=TPLL_WGS_MUT),
  name="TPLL",
  cell=cell,
  drug=NULL,
  sensitivityInfo=NULL,
  sensitivityRaw=NULL,
  sensitivityProfiles=NULL,
  sensitivityN=NULL,
  curationCell=curationCell,
  curationDrug=NULL,
  curationTissue=curationTissue,
  datasetType=c("perturbation"))

PSet@annotation$version <- 1	

saveRDS(PSet, file="TPLL_PSet.RDS")

#########################
##### 5: Using PSet #####
#########################

# 1. get summary of molecular profiles:
PSet@molecularProfiles

# 2. get summary of each molecular profile:
## for GEP
summarizeMolecularProfiles(PSet, mData='rnaseq')

## for SNP, WES and WGS mutation data (change 'mData' parameter to corresponding profile name)
summarizeMolecularProfiles(PSet, mData='snp', summary.stat='and')

# 3. summary of a molecular profile with a specific sample and  specific gene
# 'cell.lines' accespts one or more sample names
# for GEP, 'features' parameter accepts Ensembl ID(s)
gep_summary <- summarizeMolecularProfiles(PSet, mData='rnaseq', cell.lines = 'TP047', features='ENSG00000223972')
gep_summary@assays@data$exprs # displays the summary 

# for SNP, WES and WGS mutation data, 'features' parameter accepts gene names
snp_summary <- summarizeMolecularProfiles(PSet, mData='snp', summary.stat='and', cell.lines='TP047', features='ZXDB')
snp_summary@assays@data$exprs

wes_single_summary <- summarizeMolecularProfiles(PSet, mData='mut_wes_single', summary.stat='or', cell.lines=c('TP014'), features=c('A1CF', 'ABCB6'))
wes_single_summary@assays@data$exprs

# to summarize the time series data,
# 1. obtain samples that are in the time series
samples <- data.frame(PSet@molecularProfiles$mut_wes_single@colData)
timeseries_samples <- samples[!is.na(samples$timeseries_wes_mut_single), c('cellid', 'timeseries_wes_mut_single')]
timeseries_cells <- rownames(timeseries_samples[timeseries_samples$timeseries_wes_mut_single == 'TP092', ])
timeseries_summary <- summarizeMolecularProfiles(PSet, mData='mut_wes_single', summary.stat='and', cell.lines=timeseries_cells, features=c('ABCB5', 'ABCB6')) 
timeseries_summary@assays@data$exprs
