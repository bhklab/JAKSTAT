library(PharmacoGx)
library(SummarizedExperiment)
library(stringr)
library(data.table)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- paste0(args[[1]], "processed")
output_dir <- args[[1]]

# input_dir <- '/Users/minoru/Code/bhklab/JAKSTAT/TCL38_PSet-snakemake/data/processed'
output_dir <- '/Users/minoru/Code/bhklab/JAKSTAT/TCL38_PSet-snakemake/data/'

# The csv file is created with rnaseq_processing/process_abundance.R
rnaseq_samples <- read.csv(file.path(input_dir, "rnaseq_samples.csv"), row.names = 1)
rnaseq_samples$cellid <- str_replace_all(rnaseq_samples$cell, "\\W", "_")
rnaseq_samples$sampleid <- paste0(rnaseq_samples$cellid, "_", rnaseq_samples$kallisto_sample_name)
rnaseq_samples <- rnaseq_samples[order(rownames(rnaseq_samples)), ]
rownames(rnaseq_samples) <- rnaseq_samples$sampleid

# read RNA-seq data
# create SummarizedExp object for the RNA-seq data
expr_list <- readRDS(file.path(input_dir, "expr_list.rds"))
expr_annot <- read.csv(file.path(input_dir, "gene_annotation.csv"), row.names = 1)
expr_annot <- expr_annot[sort(rownames(expr_annot)), ]
expr_tpm <- as.matrix(expr_list[["expr_gene_tpm"]])
expr_tpm <- expr_tpm[, order(colnames(expr_tpm))]
colnames(expr_tpm) <- rnaseq_samples$sampleid

# add everything in SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(exprs = expr_tpm), # gene expression data in matrix
  rowData = expr_annot, # gene features Data
  colData = rnaseq_samples # patient/sample info
)
se@metadata$annotation <- "rnaseq"

raw.sensitivity <- readRDS(file.path(input_dir, "raw_sensitivity.rds"))
sensitivity_profile <- readRDS(file.path(input_dir, "sensitivity_profile.rds"))
sensitivity_info <- readRDS(file.path(input_dir, "sensitivity_info.rds"))

# create data.frame with unique sample names
cells <- sort(unique(rnaseq_samples$cellid))
curationCell <- data.frame(unique.cellid = cells)
rownames(curationCell) <- curationCell$unique.cellid

##### Create curationTissue data.frame #####
curationTissue <- data.frame(cellid = curationCell$unique.cellid, unique.tissueid = "Lymphoid")
curationTissue$tissue_type <- ifelse(grepl("cd3", curationTissue$cellid, ignore.case = TRUE) | grepl("BC", curationTissue$cellid), "CD3 Pan T cells", "T-PLL cells")
rownames(curationTissue) <- curationCell$unique.cellid

##### Create curationDrug data.frame #####
curationDrug <- data.frame(matrix(ncol = 0, nrow = (length(unique(sensitivity_info$drugid)))))
curationDrug$unique.drugid <- sort(unique(sensitivity_info$drugid))
rownames(curationDrug) <- curationDrug$unique.drugid

cell <- rnaseq_samples[, c("cellid", "sampleid")]
cells <- unique(sensitivity_info$cellid)
df <- data.frame(cellid = cells, sampleid = cells)
rownames(df) <- df$cellid
cell <- rbind(cell, df)
cell <- cell[order(rownames(cell)), ]
cell$tissueid <- "Lymphoid"
cell$tissue_type <- "Lympdoid"
rownames(cell) <- cell$sampleid
cell$cellid <- str_replace_all(cell$cellid, "_wo|_w", "")

drug <- data.frame(matrix(ncol = 0, nrow = (length(unique(sensitivity_info$drugid)))))
drug$drugid <- curationDrug$unique.drugid
rownames(drug) <- curationDrug$unique.drugid

acgh_se <- readRDS(file.path(input_dir, 'acgh_se.rds'))

tcl38 <- PharmacoGx::PharmacoSet(
  molecularProfiles = list("rnaseq" = se, 'acgh'=acgh_se),
  name = "TCL38",
  cell = cell,
  drug = drug,
  sensitivityInfo = sensitivity_info,
  sensitivityRaw = raw.sensitivity,
  sensitivityProfiles = sensitivity_profile,
  sensitivityN = as.numeric(length(colnames(raw.sensitivity[, , "Dose"]))),
  curationCell = curationCell,
  curationDrug = curationDrug,
  curationTissue = curationTissue,
  datasetType = c("both")
)

saveRDS(tcl38, paste0(output_dir, "PSet_TCL38.rds"))
