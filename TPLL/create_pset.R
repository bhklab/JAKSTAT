###load packages##
library(PharmacoGx)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

# import all the functions used in the script
source(file='functions_molecular_data.R') # functions used to process molecular data
source(file='functions_drug_screen.R') # functions used to process drug sensitivity data

##################################
##### PROCESS MOLECULAR DATA #####
##################################

##### 1: Read in all data #####
# Read all data to be included in the PSet and format patient names.

# load Ensembl reference
load("./Data/Ensembl.v75.annotation.RData") #for hg17 (Ensembl v75)

gep <- read_and_format_data("./Data/GEP_matrix_69_TPLL_vs_10_healthy_CD3_pan_T.csv")
snp <- read_and_format_data("./Data/sCNA_of_82_TPLL_and_6_FU_samples.csv")
wgs_mut <- read_and_format_data("./Data/WGS_mutations_3_TPLL_tg_pairs.csv")
wes_mut_paired <- read_and_format_data("./Data/WES_mutations_17_TPLL_tg_pairs.csv")
wes_mut_single <- read_and_format_data("./Data/WES_mutations_32_TPLL_tumor_single.csv")
wes_mut_single_followup <- read_and_format_data("./Data/WES_mutations_5_patients_with_FU.csv")

# for gene expression profiles, average out the gene expression values of duplicated ensembl id, and set ensembl gene id as rown name.
gep <- gep[!is.na(gep$Ensembl_ID),] # remove rows where Ensembl_ID is NA
gep_avg <- get_avg_gep(gep)

# remove unnecessary columns and assign gene names as row names for SNP, WES and WGS data.

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

# Order the samples by column name
gep_avg <- order_dataframes(gep_avg)
snp <- order_dataframes(snp)
wes_mut_paired <- order_dataframes(wes_mut_paired)
wes_mut_single_merged <- order_dataframes(wes_mut_single_merged)
wgs_mut <- order_dataframes(wgs_mut)


#######################################
##### PROCESS DRUG SCREENING DATA #####
#######################################

##### 1. Parse the drug screen data into raw.sensitivity array #####
raw.sensitivity <- get_raw_sensitivity('./Data/drug_screen/', c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5', 'Dose6', 'Dose7'))

##### 2. Create sensitivity_info #####
sensitivity_info <- get_sensitivity_info(raw.sensitivity)

##### 3. Create sensitivity_profile #####
sensitivity_profile <- get_sensitivity_profile(raw.sensitivity)


#######################################
##### CREATE CURATION DATAFRAMES ######
#######################################

# NOTE: Each "cell" represents each sample
# Concatenate all sample names from molecular profiles + drug sensitivity data into one vector
total_samples <- unique(c(colnames(gep_avg), colnames(snp), colnames(wes_mut_paired), colnames(wes_mut_single_merged), colnames(wgs_mut), unique(sensitivity_info$cellid)))
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
curationTissue$tissue_type <- ifelse(grepl('cd3', curationTissue$cellid, ignore.case=TRUE), 'CD3 Pan T cells', 'T-PLL cells')
rownames(curationTissue) <- curationCell$unique.cellid

##### Create Cell data.frame #####
# rownames must be unique sample names. All samples names from molecular profiles + drug sensitivity data must be included here.
# this is the reference dataframe used to dictate which smaples are included in the PSet. Any data samples that are not in this data frame will be removed.
cell <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) #ncol depends on number of sample features to include.
cell$cellid <- curationCell$unique.cellid
cell$tissueid <- curationTissue$unique.tissueid
cell$tissue_type <- curationTissue$tissue_type
rownames(cell) <- curationCell$unique.cellid

##### Create curationDrug data.frame #####
curationDrug <- data.frame(matrix(ncol=0, nrow=(length(unique(sensitivity_info$drugid)))))
curationDrug$unique.drugid <- sort(unique(sensitivity_info$drugid))
rownames(curationDrug) <- curationDrug$unique.drugid

##### Create Drug data.frame #####
drug <- data.frame(matrix(ncol=0, nrow=(length(unique(sensitivity_info$drugid)))))
drug$drugid <- curationDrug$unique.drugid
rownames(drug) <- curationDrug$unique.drugid

################################################################
##### CREATE SummarizedExperiments FOR MOLECULAR PROFILES ######
################################################################

##### 1 Prepare the clinical_data data frame to be used as colData in each SummarizedExperiment object #####
clinical_data <- read_and_format_clinical_data(
  path = "./Data/clinical_data_2020_10_23.csv", 
  all_samples = rownames(cell),
  data_colnames = list(colnames(snp), colnames(wes_mut_single_merged)), 
  col_names = c('timeseries_snp', 'timeseries_wes_mut_single'))

##### 2 Create Gene Expression Data SummarizedExperiment #####
TPLL_GEP <- get_summarized_experiment(gep_avg, clinical_data, features_gene, 'gene_id', 'rnaseq')

##### 3 Create SNP Data SummarizedExperiment #####
TPLL_SNP <- get_summarized_experiment(snp, clinical_data, features_gene, 'gene_name', 'mutation')

##### 4 Create WGS (Whole Genome Sequencing) mutation data SummarizedExperiment #####
wgs_mut[wgs_mut == 'unmutated'] <- 'wt'
TPLL_WGS_MUT <- get_summarized_experiment(wgs_mut, clinical_data, features_gene, 'gene_name', 'mutation')

##### 5 Create WES (Whole Exome Sequencing) single-end and paired-end mutation data SummarizedExperiments#####
### Single-end Mutation Data ###
wes_mut_single_merged[wes_mut_single_merged == 'unmutated'] <- 'wt'
TPLL_WES_MUT_SINGLE <- get_summarized_experiment(wes_mut_single_merged, clinical_data, features_gene, 'gene_name', 'mutation')

##### 6 Paied-end Mutation Data #####
wes_mut_paired[wes_mut_paired == 'unmutated'] <- 'wt'
TPLL_WES_MUT_PAIRED <- get_summarized_experiment(wes_mut_paired, clinical_data, features_gene, 'gene_name', 'mutation')


#######################
##### CREATE PSET #####
#######################

PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=list(
    "rnaseq"=TPLL_GEP, 
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

PSet@annotation$version <- 1	

saveRDS(PSet, file="TPLL_PSet.RDS")
