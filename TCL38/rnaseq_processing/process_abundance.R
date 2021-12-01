library(stringr)
library(dplyr)
library(MultiAssayExperiment)

in_dir <- "./data/"

sample <- readRDS(paste0(in_dir, "multiassay_sample/TCGA_ACC.rds"))