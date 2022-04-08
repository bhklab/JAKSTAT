###load packages##
library(PharmacoGx)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

# Before running this script, run the following scripts to get prerequisite data objects:
# 1. get_old_sensitivity_data.R and get_new_sensitivity_data.R to get raw drig sensitivity data.
# 2. get_molecular_data.R to get all molecular data summarized experiment and patient metadata.

source(file='functions_drug_screen.R') # functions used to process drug sensitivity data

#######################################
##### PROCESS DRUG SCREENING DATA #####
#######################################

# load data objects created by running get_old_sensitivity_data.R and get_new_sensitivity_data.R
old_sensitivity <- readRDS('./Data/curated_data/old_raw_sensitivity.rds')
new_sensitivity <- readRDS('./Data/curated_data/new_raw_sensitivity.rds')

# remove duplicate sensitivity data from the old data.
# identical_rows <- intersect(rownames(old_sensitivity), rownames(new_sensitivity))

rows_to_keep <- setdiff(rownames(old_sensitivity), rownames(new_sensitivity))
old_sensitivity <- old_sensitivity[rownames(old_sensitivity) %in% rows_to_keep, ]
raw.sensitivity <- rbind(new_sensitivity, old_sensitivity)
raw.sensitivity <- raw.sensitivity[order(rownames(raw.sensitivity)),]

conc_tested <- 5

# Make the merged dataframe into an array
raw.sensitivity <- array(
  c(
    as.matrix(raw.sensitivity[ ,1:conc_tested]), 
    as.matrix(raw.sensitivity[ ,(conc_tested + 1):(2 * conc_tested)])
  ), 
  c(nrow(raw.sensitivity), conc_tested , 2), 
  dimnames=list(
    rownames(raw.sensitivity), 
    colnames(raw.sensitivity[ ,1:conc_tested]), 
    c("Dose", "Viability")
  )
)

sensitivity_info <- get_sensitivity_info(raw.sensitivity)
sensitivity_profile <- get_sensitivity_profile(raw.sensitivity)

#######################################
##### CREATE CURATION DATAFRAMES ######
#######################################

# load objects that have been created by running get_molecular_data.R
# sample metadata
clinical_sample_data <- readRDS('./Data/curated_data/clinical_sample_data.rds')

# molecular data
TPLL_GEP <- readRDS('./Data/curated_data/TPLL_GEP.rds')
TPLL_mRNA <- readRDS('./Data/curated_data/TPLL_mRNA.rds')
TPLL_miRNA <- readRDS('./Data/curated_data/TPLL_miRNA.rds')
TPLL_SNP <- readRDS('./Data/curated_data/TPLL_SNP.rds')
TPLL_WGS_MUT <- readRDS('./Data/curated_data/TPLL_WGS_MUT.rds')
TPLL_WES_MUT_SINGLE <- readRDS('./Data/curated_data/TPLL_WES_MUT_SINGLE.rds')
TPLL_WES_MUT_PAIRED <- readRDS('./Data/curated_data/TPLL_WES_MUT_PAIRED.rds')

# NOTE: Each "cell" represents each patient
# Concatenate all patient(TP) numbers from molecular profiles + drug sensitivity data into one vector
total_samples <- unique(c(
  unique(clinical_sample_data$cellid),
  unique(sensitivity_info$cellid))
)
total_samples <- sort(total_samples)

# create data.frame with unique sample names
curationCell <- data.frame(matrix(ncol=1, nrow=(length(total_samples))))
colnames(curationCell) <- c("unique.cellid")
curationCell$unique.cellid <- total_samples
rownames(curationCell) <- curationCell$unique.cellid

##### Create curationTissue data.frame #####
curationTissue <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) 
curationTissue$cellid <- curationCell$unique.cellid
curationTissue$unique.tissueid <- rep('Lymphoid', length(total_samples))
# curationTissue$tissue_type <- ifelse(startsWith(curationTissue$cellid, 'TP'), 'T-PLL cells', 'CD3 Pan T cells')
curationTissue$tissue_type <- ifelse(grepl('cd3', curationTissue$cellid, ignore.case=TRUE) | grepl('BC', curationTissue$cellid), 'CD3 Pan T cells', 'T-PLL cells')
rownames(curationTissue) <- curationCell$unique.cellid

##### Create Cell data.frame #####
# rownames must be unique sample names. All samples names from molecular profiles + drug sensitivity data must be included here.
# this is the reference dataframe used to dictate which samples are included in the PSet. Any data samples that are not in this data frame will be removed.
cell <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) #ncol depends on number of sample features to include.
cell$cellid <- curationCell$unique.cellid
cell$tissueid <- curationTissue$unique.tissueid
cell$tissue_type <- curationTissue$tissue_type
rownames(cell) <- curationCell$unique.cellid

# Prepare the patient_data dataframe to be added to cells
patient_data <- read.csv("./Data/clinical_data_2020_10_23.csv") 
rownames(patient_data) <- patient_data[,1]
patient_data <- patient_data[-c(1)]
patient_data <- data.frame(t(patient_data)) # transpose patient_data with patient id as rows and attirbutes as columns.

# add any missing samples (patients and healthy individuals) into patient_data and populate the columns with NA
# this is done to ensure that all the samples in each SummarizedExperiment are included and not ommitted.
patient_missing <- total_samples[!total_samples %in% rownames(patient_data)]
patient_data[patient_missing,] <- NA

# make sure that the patient_data is ordered alphabetically by rowname
patient_data <- patient_data[order(rownames(patient_data)), ]

# add clinical data to cell dataframe
for (col in colnames(patient_data)){
  cell[, col] <- patient_data[, col][match(row.names(cell), row.names(patient_data))]
}

##### Create curationDrug data.frame #####
curationDrug <- data.frame(matrix(ncol=0, nrow=(length(unique(sensitivity_info$drugid)))))
curationDrug$unique.drugid <- sort(unique(sensitivity_info$drugid))
rownames(curationDrug) <- curationDrug$unique.drugid

##### Create Drug data.frame #####
drug <- data.frame(matrix(ncol=0, nrow=(length(unique(sensitivity_info$drugid)))))
drug$drugid <- curationDrug$unique.drugid
rownames(drug) <- curationDrug$unique.drugid

all(rownames(sensitivity_info) == rownames(sensitivity_profile) & 
      rownames(sensitivity_info) == dimnames(raw.sensitivity)[[1]])



#######################
##### CREATE PSET #####
#######################

PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=list(
    "gep"=TPLL_GEP, 
    "rnaseq"=TPLL_mRNA,
    "mi_rnaseq"=TPLL_miRNA,
    "snp"=TPLL_SNP,
    "mut_wes_single"=TPLL_WES_MUT_SINGLE, 
    "mut_wes_paired"=TPLL_WES_MUT_PAIRED,
    "mut_wgs"=TPLL_WGS_MUT),
  name="TPLL",
  cell=cell,
  drug=drug,
  sensitivityInfo=sensitivity_info,
  sensitivityRaw=raw.sensitivity,
  sensitivityProfiles=sensitivity_profile,
  sensitivityN=as.numeric(length(colnames(raw.sensitivity[,,'Dose']))),
  curationCell=curationCell,
  curationDrug=curationDrug,
  curationTissue=curationTissue,
  datasetType=c("both"))

PSet@annotation$version <- 2	

saveRDS(PSet, file="./Data/curated_data/TPLL_PSet.rds")
