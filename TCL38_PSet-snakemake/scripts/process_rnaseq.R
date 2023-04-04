library(stringr)
library(dplyr)
library(data.table)
library(tximport)
library(utils)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

# input_dir <- '/Users/minoru/Code/bhklab/JAKSTAT/TCL38_PSet-snakemake/data/rnaseq'

# Annotation data for the transcripts
load(file.path(input_dir, "Gencode.v33.annotation.RData"))

# The abundance.tsv files for each samples were obtained from /repository/ukoeln_herling/procdata/human/rna_seq/202111/output
sample_dir <- file.path(input_dir, 'rnaseq')
dir.create(sample_dir)
unzip(zipfile=file.path(input_dir, 'rnaseq.zip'), exdir = sample_dir)

samples <- list.dirs(sample_dir, full.names=FALSE, recursive = FALSE)

files <- file.path(sample_dir, samples, "abundance.h5")
names(files) <- samples

expr_tx <- tximport(files, type = "kallisto", txOut = TRUE)
expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene)

expr_list <- list()
expr_list[['expr_gene_tpm']] <- log2(expr_gene$abundance + 0.001)
expr_list[['expr_gene_counts']] <- log2(expr_gene$counts + 1)
expr_list[['expr_isoform_tpm']] <- log2(expr_tx$abundance + 0.001)
expr_list[['expr_isoform_counts']] <- log2(expr_tx$counts + 1)

# coldata
summary_df <- read.csv(file.path(input_dir, "rnaseq_sample_summary.csv"))
coldata <- summary_df[c("Kallisto.Sample.Name", "Sample.Name", "Cell")]
colnames(coldata) <- c("kallisto_sample_name", "sample_name", "cell")
rerun <- summary_df[!is.na(summary_df$Kallisto.Rerun.Sample), c("Kallisto.Rerun.Sample", "Sample.Name", "Cell")]
colnames(rerun) <- c("kallisto_sample_name", "sample_name", "cell")
rerun$rerun <- TRUE
coldata <- dplyr::bind_rows(coldata, rerun)
coldata <- coldata[order(coldata$kallisto_sample_name), ]
rownames(coldata) <- str_remove_all(coldata$kallisto_sample_name, "[\\s\\t\\n]")
coldata$healthy <- unlist(lapply(coldata$cell, function(cell){if(cell == "BC"){return(TRUE)} else {return(NA)}}))

# output the curated data
saveRDS(expr_list, file.path(output_dir, 'expr_list.rds'))
write.csv(features_gene, file.path(output_dir, "gene_annotation.csv"), row.names = TRUE)
write.csv(features_transcript, file.path(output_dir, "tx_annotation.csv"), row.names = TRUE)
write.csv(coldata, file.path(output_dir, "rnaseq_samples.csv"), row.names = TRUE)

unlink(file.path(input_dir, 'rnaseq'), recursive=TRUE)
unlink(file.path(input_dir, '__MACOSX'), recursive=TRUE)