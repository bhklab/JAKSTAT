# Parses the ACGH data into SummarizedExperiment as it is the requireped format for the PSet.

library(stringr)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]

# work_dir <- '/Users/minoru/Code/bhklab/JAKSTAT/TCL38_PSet-snakemake/data/processed'

# create sample dataframe
# The csv file is created with aCGH_processing/process_acgh.R
acgh_samples <- read.csv(file.path(work_dir, "acgh_samples.csv"), row.names=1) 
# The csv file is created with rnaseq_processing/process_abundance.R
rnaseq_samples <- read.csv(file.path(work_dir, "rnaseq_samples.csv"), row.names=1) 

# The RDS is created with aCGH_processing/process_acgh.R.
acgh_data <- readRDS(file.path(work_dir, "acgh_assay_data.rds"))

cell_map <- data.frame(cellid=unique(rnaseq_samples$cell), acgh_cells=NA)
cell_map$acgh_cells <- unlist(lapply(cell_map$cellid, function(cellid){
  for(cell in acgh_samples$cell){
    cellid_comp <- str_to_lower(cellid)
    cellid_comp <- str_replace_all(cellid_comp, '\\W', '')
    cell_comp <- str_to_lower(cell)
    cell_comp <- str_replace_all(cell_comp, '\\W', '')
    if(cellid_comp == cell_comp){
      return(cell)
    }
  }
  return(NA)
}))
cell_map[cell_map$cellid == 'CCRF-H-SB2', c('acgh_cells')] <- 'H-SB2'
cell_map[cell_map$cellid == 'OCI-Ly13.2', c('acgh_cells')] <- 'OCI_13'
cell_map[cell_map$cellid == 'Karpas299', c('acgh_cells')] <- 'K299'
cell_map[cell_map$cellid == 'Karpas384', c('acgh_cells')] <- 'Karpas_384'
cell_map <- cell_map[!is.na(cell_map$acgh_cells), ]
cell_map$cellid <- str_replace_all(cell_map$cellid, "\\W", "_")

acgh_samples$cell <- unlist(lapply(acgh_samples$cell, function(cell){
  return(cell_map[cell_map$acgh_cells == cell, c('cellid')])
}))
names(acgh_data) <- unlist(lapply(names(acgh_data), function(cell){
  return(cell_map[cell_map$acgh_cells == cell, c('cellid')])
}))

assay_rownames <- c()
new_colnames <- c(
  'cellid', 
  'abberation_num',
  'chr',
  'cytoband', 
  'start', 
  'stop', 
  'probes', 
  'amplification', 
  'deletion', 
  'pval', 
  'gene_names', 
  'hg19_cnv', 
  'hg19_mirna', 
  'mm8_mirna'
)
merged_df <- data.frame(matrix(ncol=length(new_colnames)))
colnames(merged_df) <- new_colnames
for(sample in names(acgh_data)){
  col_names <- colnames(acgh_data[[sample]])
  acgh_data[[sample]]$cellid <- sample
  rownames(acgh_data[[sample]]) <- paste0(sample, '_', acgh_data[[sample]]$AberrationNo)
  acgh_data[[sample]] <- acgh_data[[sample]][, c('cellid', col_names)]
  colnames(acgh_data[[sample]]) <- new_colnames
  merged_df <- rbind(merged_df, acgh_data[[sample]])
}
merged_df <- merged_df[rownames(merged_df) != '1', ]

rownames(acgh_samples) <- acgh_samples$cell

assay_df <- data.frame(matrix(nrow=length(rownames(merged_df)), ncol=length(rownames(acgh_samples))))
rownames(assay_df) <- rownames(merged_df)
colnames(assay_df) <- rownames(acgh_samples)
for(cell in rownames(acgh_samples)){
  rows <- rownames(merged_df[merged_df$cellid == cell, ])
  assay_df[rownames(assay_df) %in% rows, c(cell)] <- 'yes'
}

acgh_se <- SummarizedExperiment(
  assays=list(assay=assay_df),
  rowData=merged_df,
  colData = acgh_samples
)
acgh_se@metadata$annotation <- "acgh"

saveRDS(acgh_se, file.path(work_dir, 'acgh_se.rds'))