library(stringr)
library(dplyr)
library(MultiAssayExperiment)

in_dir <- "./data/rnaseq_abundance/"

get_abundance_df <- function(dir, samples, type){
  abundance_df <- data.frame(matrix(nrow=0, ncol=0))
  for(sample in samples){
    read_in <- read.delim(paste0(dir, sample, "/abundance.tsv"), header=TRUE, sep='\t')
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

samples <- list.dirs(in_dir, full.names=FALSE, recursive = FALSE)
counts_df <- get_abundance_df(in_dir, samples, "est_counts")
tpm_df <- get_abundance_df(in_dir, samples, "tpm")