library(stringr)
library(dplyr)
library(MultiAssayExperiment)

load("./Data/Ensembl.v75.annotation.RData") #for hg17 (Ensembl v75)

in_dir <- "./data/rnaseq_abundance/"

get_abundance_df <- function(dir, samples, type){
  abundance_df <- data.frame(matrix(nrow=0, ncol=0))
  for(sample in samples){
    read_in <- read.delim(paste0(dir, sample, "/abundance.tsv"), header=TRUE, sep='\t')
    # formatted_rownames <- str_remove(read_in$target_id, "\\..+$")
    rownames(read_in) <- read_in$target_id
    colnames(read_in)[colnames(read_in) == type] <- sample
    
    if(dim(abundance_df)[1] == 0){
      abundance_df <- data.frame(matrix(nrow=dim(read_in[1]), ncol=0))
      rownames(abundance_df) <- read_in$target_id
    }
    abundance_df <- cbind(abundance_df, read_in[, sample][match(rownames(abundance_df), rownames(read_in))])
  }
  colnames(abundance_df) <- samples
  return(abundance_df)
}

get_rowdata_df <- function(abundance_df, transcript_annotation){
  rowdata <- data.frame(matrix(nrow=length(rownames(abundance_df)), ncol=length(colnames(transcript_annotation))))
  rownames(rowdata) <- rownames(abundance_df)
  colnames(rowdata) <- colnames(transcript_annotation)
  
  for(col_name in colnames(rowdata)){
    print(col_name)
    values <- lapply(rownames(rowdata), function(row_name){
      value <- transcript_annotation[[col_name]][transcript_annotation$transcript_id == str_remove(row_name, "\\..+$")]
      return(ifelse(length(value) == 0, NA, value))
    })
    rowdata[[col_name]] <- unlist(values)
  }
  return(rowdata)
}

samples <- list.dirs(in_dir, full.names=FALSE, recursive = FALSE)
counts_df <- get_abundance_df(in_dir, samples, "est_counts")
tpm_df <- get_abundance_df(in_dir, samples, "tpm")

# coldata
summary_df <- read.csv("./data/rnaseq_sample_summary.csv")
coldata <- summary_df[c("Kallisto.Sample.Name", "Sample.Name", "Cell")]
colnames(coldata) <- c("kallisto_sample_name", "sample_name", "cell")
rerun <- summary_df[!is.na(summary_df$Kallisto.Rerun.Sample), c("Kallisto.Rerun.Sample", "Sample.Name", "Cell")]
colnames(rerun) <- c("kallisto_sample_name", "sample_name", "cell")
rerun$rerun <- TRUE
coldata <- dplyr::bind_rows(coldata, rerun)
coldata <- coldata[order(coldata$kallisto_sample_name), ]
rownames(coldata) <- coldata$kallisto_sample_name
coldata$healthy <- unlist(lapply(coldata$cell, function(cell){if(cell == "BC"){return(TRUE)} else {return(NA)}}))

# rowdata
rowdata_counts <- rowdata
rowdata_tpm <- get_rowdata_df(tpm_df, features_transcript)
# write.csv(rowdata_tpm, "./data/rowdata_tpm.csv")

