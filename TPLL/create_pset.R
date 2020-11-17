###load packages##
library(PharmacoGx)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

# import all the functions used in the script
source(file='functions.R')

##################################
##### Prep: Read in all data #####
##################################

#1. Read all data to be included in the PSet and format patient names.

# load Ensembl reference
load("./Data/Ensembl.v75.annotation.RData") #for hg17 (Ensembl v75)

gep <- read_and_format_data("./Data/GEP_matrix_69_TPLL_vs_10_healthy_CD3_pan_T.csv")
snp <- read_and_format_data("./Data/sCNA_of_82_TPLL_and_6_FU_samples.csv")
wgs_mut <- read_and_format_data("./Data/WGS_mutations_3_TPLL_tg_pairs.csv")
wes_mut_paired <- read_and_format_data("./Data/WES_mutations_17_TPLL_tg_pairs.csv")
wes_mut_single <- read_and_format_data("./Data/WES_mutations_32_TPLL_tumor_single.csv")
wes_mut_single_followup <- read_and_format_data("./Data/WES_mutations_5_patients_with_FU.csv")

#2. for gene expression profiles, average out the gene expression values of duplicated ensembl id, and set ensembl gene id as rown name.
gep <- gep[!is.na(gep$Ensembl_ID),] # remove rows where Ensembl_ID is NA
gep_avg <- get_avg_gep(gep)

#3. remove unnecessary columns and assign gene names as row names for SNP, WES and WGS data.

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


#############################################
##### 1: Create curationCell data.frame #####
#############################################

# NOTE: Each "cell" represents each smaple
# Concatenate all sample names from molecular profiles + drug sensitivity data into one vector
total_samples <- unique(c(colnames(gep_avg), colnames(snp), colnames(wes_mut_paired), colnames(wes_mut_single_merged), colnames(wgs_mut)))
total_samples <- sort(total_samples)

# create data.frame with unique sample names
curationCell <- data.frame(matrix(ncol=1, nrow=(length(total_samples))))
colnames(curationCell) <- c("unique.cellid")
curationCell$unique.cellid <- total_samples
rownames(curationCell) <- curationCell$unique.cellid


#############################################
##### 2: Create curationTissue data.frame #####
#############################################
curationTissue <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) 
curationTissue$cellid <- curationCell$unique.cellid
curationTissue$unique.tissueid <- rep('Lymphoid', length(total_samples))
curationTissue$tissue_type <- ifelse(startsWith(curationTissue$cellid, 'TP'), 'T-PLL cells', 'CD3 Pan T cells')
rownames(curationTissue) <- curationCell$unique.cellid


#####################################
##### 3: Create Cell data.frame #####
#####################################

# rownames must be unique sample names. All samples names from molecular profiles + drug sensitivity data must be included here.
# this is the reference dataframe used to dictate which smaples are included in the PSet. Any data samples that are not in this data frame will be removed.
cell <- data.frame(matrix(ncol=0, nrow=(length(total_samples)))) #ncol depends on number of sample features to include.
cell$cellid <- curationCell$unique.cellid
cell$tissueid <- curationTissue$unique.tissueid
cell$tissue_type <- curationTissue$tissue_type
rownames(cell) <- curationCell$unique.cellid


################################################################
##### 4: Create summarizedExperiments (molecular profiles) #####
################################################################

##### 4.0.1 Prepare the clinical_data data frame to be used as colData in each SummarizedExperiment object #####
clinical_data <- read_and_format_clinical_data(
  path = "./Data/clinical_data_2020_10_23.csv", 
  all_samples = rownames(cell),
  data_colnames = list(colnames(snp), colnames(wes_mut_single_merged)), 
  col_names = c('timeseries_snp', 'timeseries_wes_mut_single'))

##### 4.1 Create Gene Expression Data SummarizedExperiment #####
TPLL_GEP <- get_summarized_experiment(gep_avg, clinical_data, features_gene, 'gene_id', 'rnaseq')

##### 4.2 Create SNP Data SummarizedExperiment #####
TPLL_SNP <- get_summarized_experiment(snp, clinical_data, features_gene, 'gene_name', 'mutation')

##### 4.3 Create WGS (Whole Genome Sequencing) mutation data SummarizedExperiment #####
wgs_mut[wgs_mut == 'unmutated'] <- 'wt'
TPLL_WGS_MUT <- get_summarized_experiment(wgs_mut, clinical_data, features_gene, 'gene_name', 'mutation')

##### 4.4 Create WES (Whole Exome Sequencing) single-end and paired-end mutation data SummarizedExperiments#####
### Single-end Mutation Data ###
wes_mut_single_merged[wes_mut_single_merged == 'unmutated'] <- 'wt'
TPLL_WES_MUT_SINGLE <- get_summarized_experiment(wes_mut_single_merged, clinical_data, features_gene, 'gene_name', 'mutation')

### Paied-end Mutation Data ###
wes_mut_paired[wes_mut_paired == 'unmutated'] <- 'wt'
TPLL_WES_MUT_PAIRED <- get_summarized_experiment(wes_mut_paired, clinical_data, features_gene, 'gene_name', 'mutation')

##########################
##### 4: Create PSet #####
##########################

PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=list(
    "rnaseq"=TPLL_GEP, 
    "snp"=TPLL_SNP,
    "mut_wes_single"=TPLL_WES_MUT_SINGLE, 
    "mut_wes_paired"=TPLL_WES_MUT_PAIRED,
    "mut_wgs"=TPLL_WGS_MUT),
  name="TPLL",
  cell=cell,
  drug=NULL,
  sensitivityInfo=NULL,
  sensitivityRaw=NULL,
  sensitivityProfiles=NULL,
  sensitivityN=NULL,
  curationCell=curationCell,
  curationDrug=NULL,
  curationTissue=curationTissue,
  datasetType=c("perturbation"))

PSet@annotation$version <- 1	

saveRDS(PSet, file="TPLL_PSet.RDS")


# load('./Data/drug_screen_processed/201810_p1383/p1383_response data.Rdata')
# drug_meta <- read.csv('./Data/drug_screen_processed/201810_p1383/p1383_drug_metadata.csv')
# patient_meta <- read.csv('./Data/drug_screen_processed/201810_p1383/p1383_patient_metadata.csv')
