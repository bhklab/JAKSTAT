library(stringr)
library(dplyr)

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

get_rowdata_df <- function(abundance_df, annotation_df){
  rowdata <- annotation_df[rownames(annotation_df) %in% rownames(abundance_df), ]
  missing <- setdiff(rownames(abundance_df), rownames(features_transcript))
  added_df <- data.frame(matrix(NA, ncol = ncol(rowdata), nrow = length(missing), dimnames = list(missing, names(rowdata))))
  added_df$transcript_id <- rownames(added_df)
  rowdata <- dplyr::bind_rows(rowdata, added_df)
  return(rowdata)
}

get_df_with_averaged_duplicates <- function(abundance_df){
  abundance_transcripts <- str_remove(rownames(abundance_df), "\\..+$")
  multiple <- data.frame(table(abundance_transcripts))
  multiple <- multiple[multiple$Freq > 1, ]
  abundance_df$transcript_id <- str_remove(rownames(abundance_df), "\\..+$")
  
  transcripts <- as.list(as.character(multiple$abundance_transcripts))
  names(transcripts) <- multiple$abundance_transcripts
  averages <- sapply(transcripts, function(transcript){
    duplicate <- abundance_df[abundance_df$transcript_id == transcript, ]
    duplicate$transcript_id <- NULL
    average <- lapply(colnames(duplicate), function(col_name){
      return(as.numeric(formatC(mean(duplicate[[col_name]]), format="e", digits=5)))
    })
    average <- unlist(average)
    return(average)
  }, simplify = FALSE,USE.NAMES = TRUE)
  
  average_df <- data.frame(matrix(nrow=length(averages), ncol=length(colnames(abundance_df)[colnames(abundance_df) != "transcript_id"])))
  rownames(average_df) <- names(averages)
  colnames(average_df) <- colnames(abundance_df)[colnames(abundance_df) != "transcript_id"]
  for(row_name in rownames(average_df)){
    average_df[rownames(average_df) == row_name, ] <- averages[[row_name]]
  }
  
  abundance_df <- dplyr::bind_rows(abundance_df, average_df)
  abundance_df <- abundance_df[!(abundance_df$transcript_id %in% rownames(average_df)), ]
  rownames(abundance_df) <- str_remove(rownames(abundance_df), "\\..+$")
  abundance_df$transcript_id <- NULL
  return(abundance_df)
}

samples <- list.dirs(in_dir, full.names=FALSE, recursive = FALSE)
counts_df <- get_abundance_df(in_dir, samples, "est_counts")
tpm_df <- get_abundance_df(in_dir, samples, "tpm")

# remove duplicate transcripts by taking mean of the values, so that the transcript ids can be used row names.
counts_df_average <- get_df_with_averaged_duplicates(counts_df)
tpm_df_average <- get_df_with_averaged_duplicates(tpm_df)

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
rowdata_counts <- get_rowdata_df(counts_df_average, features_transcript)
rowdata_tpm <- get_rowdata_df(tpm_df_average, features_transcript)

# output the curated data
write.csv(counts_df, "../create_multiassay/data/counts.csv")
write.csv(tpm_df, "../create_multiassay/data/tpm.csv")
write.csv(counts_df_average, "../create_multiassay/data/counts_averaged.csv")
write.csv(tpm_df_average, "../create_multiassay/data/tpm_averaged.csv")
write.csv(rowdata_counts, "../create_multiassay/data/counts_annotation.csv")
write.csv(rowdata_tpm, "../create_multiassay/data/tpm_annotation.csv")
write.csv(coldata, "../create_multiassay/data/rnaseq_samples.csv")
