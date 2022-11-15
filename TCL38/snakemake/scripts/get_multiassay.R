library(stringr)
library(dplyr)
library(MultiAssayExperiment)
library(Biobase)
library(GenomicRanges)
library(RaggedExperiment)

#' To run this script, 
#' you need the following files by running aCGH_processing/process_acgh.R
#' 1. acgh_samples.csv
#' 2. acgh_assay_data.rds
#' and the following files by running rnaseq_processing/process_abundance.R
#' 1. rnaseq_samples.csv
#' 2. gene_annotation.csv
#' 3. tx_annotation.csv
#' 4. expr_list.rds

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
 
# input_dir <- '/Users/minoru/Documents/bhklab/jakstat/tcl38/data/multiassay_exp'

# create sample dataframe
# The csv file is created with aCGH_processing/process_acgh.R
acgh_samples <- read.csv(file.path(input_dir, "acgh_samples.csv"), row.names=1) 
# The csv file is created with rnaseq_processing/process_abundance.R
rnaseq_samples <- read.csv(file.path(input_dir, "rnaseq_samples.csv"), row.names=1) 

acgh_cells <- unique(acgh_samples$cell)
rnaseq_cells <- unique(rnaseq_samples$cell)

# create sample info df
samples_df <- as.data.frame(matrix(nrow=39, ncol=2, dimnames=list(NULL, c("acgh_cells", "rnaseq_cells"))))
samples_df$rnaseq_cells <- rnaseq_cells

# map cell names in aCGH and RNA-seq data, and come up with a common name as the row name.
cellname_map <- list(`CCRF-H-SB2`="H-SB2", `OCI-Ly13.2`="OCI_13", `Karpas299`="K299", `BC`=NA)
for(cellname in rnaseq_cells){
  found <- acgh_cells[str_detect(tolower(str_replace_all(acgh_cells, "[-_\\.]", "")), tolower(str_replace_all(cellname, "[-_\\.]", "")))]
  samples_df[samples_df$rnaseq_cells == cellname, ]$acgh_cells <- ifelse(length(found) > 0, found, cellname_map[[cellname]])
}
rownames(samples_df) <- str_replace_all(samples_df$rnaseq_cells, "[-\\.]", "_")

# add RNA-seq samples data
samples_df$rnaseq_sample_names <- lapply(samples_df$rnaseq_cells, function(cell){
  return(paste(rnaseq_samples[rnaseq_samples$cell == cell, ]$sample_name, collapse=", "))
})
samples_df$rnaseq_kallisto_sample_names <- lapply(samples_df$rnaseq_cells, function(cell){
  return(paste(rnaseq_samples[rnaseq_samples$cell == cell, ]$kallisto_sample_name, collapse=", "))
})
samples_df$rnaseq_healthy <- lapply(samples_df$rnaseq_cells, function(cell){
  ifelse(cell == "BC", "Yes", "")
})

# merge aCGH samples data
samples_df[colnames(acgh_samples)[colnames(acgh_samples) != "cell"]] <- NA
for(cell in acgh_cells){
  values <- as.list(acgh_samples[acgh_samples$cell == cell, -which(colnames(acgh_samples) == "cell")])
  for(name in names(values)){
    samples_df[which(samples_df$acgh_cells == cell), which(colnames(samples_df) == name)] <- values[[name]]
  }
}
colnames(samples_df)[colnames(samples_df) %in% colnames(acgh_samples)] <- paste0("acgh_", colnames(samples_df)[colnames(samples_df) %in% colnames(acgh_samples)])

# Format assay data
assay_list <- list()

# read RNA-seq data
# create ExpressionSet objects for the RNA-seq data
expr_pheno <- as(rnaseq_samples, "AnnotatedDataFrame")
expr_list <- readRDS(file.path(input_dir, 'expr_list.rds'))
for(expr_name in names(expr_list)){
  assay_list[[expr_name]] <- ExpressionSet(assayData=as.matrix(expr_list[[expr_name]]), phenoData=expr_pheno)
}

# create GenomicRanges objects for aCGH data, and parse them into a RaggedExperiment object
# The RDS is created with aCGH_processing/process_acgh.R.
acgh_data <- readRDS(file.path(input_dir, "acgh_assay_data.rds"))
genomic_ranges_list <- GRangesList()
for(sample in names(acgh_data)){
  tmp <- acgh_data[[sample]]
  genomic_ranges_list[[sample]] <- GRanges(
    seqnames=tmp$Chr,
    ranges=IRanges(start=as.numeric(tmp$Start), end=as.numeric(tmp$Stop), names=paste0(sample, "_", tmp$AberrationNo)),
    strand="*",
    abberation_num=tmp$AberrationNo,
    cytoband=tmp$Cytoband,
    probes=tmp$Probes,
    amplification=tmp$Amplification,
    deletion=tmp$Deletion,
    pval=tmp$pval,
    gene_names=tmp$Gene_Names,
    hg19_cnv=tmp$Hs_hg19_CNV_20120403,
    hg19_mirna=tmp$Hs_hg19_miRNA_20120403,
    mm8_mirna=tmp$Mm_mm8_miRNA_20090908
  )
}
acgh_coldata <- DataFrame(id=acgh_samples$cell)
assay_list$acgh <- RaggedExperiment(
  genomic_ranges_list,
  colData=acgh_coldata
)

# create sampleMap dataframe
sample_map_df <- data.frame()
cells <- c()
samples <- c()
for(sample in rownames(expr_pheno)){
  cell <- rownames(samples_df[which(samples_df$rnaseq_cells == rnaseq_samples[sample, ]$cell), ])
  cells <- c(cells, cell)
  samples <- c(samples, sample)
}

expr_data_map <- data.frame(colname=samples, primary=cells, stringsAsFactors = FALSE)
expr_data_map$assay <- NA 
for(expr_name in names(expr_list)){
  expr_data_map$assay <- expr_name
  sample_map_df <- dplyr::bind_rows(sample_map_df, expr_data_map)
}

# add aCGH data sample info to the sample map
acgh_sample_map <- as.data.frame(matrix(nrow=length(names(acgh_data)), ncol=3, dimnames=list(NULL, c("colname", "primary", "assay"))), stringsAsFactors = FALSE)
acgh_sample_map$colname <- names(acgh_data)
cells <- lapply(names(acgh_data), function(sample){
  return(rownames(samples_df[which(samples_df$acgh_cells == sample), ]))
})
acgh_sample_map$primary <- unlist(cells)
acgh_sample_map$assay <- "acgh"
sample_map_df <- dplyr::bind_rows(sample_map_df, acgh_sample_map)

# Create MultiAssayExp object
multiAssay <- MultiAssayExperiment(
  assay_list, 
  samples_df, 
  sample_map_df
)

saveRDS(multiAssay, file.path(output_dir, "TCL38_MultiAssayExp.rds"))
