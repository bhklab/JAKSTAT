library(SummarizedExperiment)
library(dplyr)
library(stringr)

source(file='functions_molecular_data.R') # functions used to process molecular data

##################################
##### PROCESS MOLECULAR DATA #####
##################################

##### 1: Read in all data #####
# Read all data to be included in the PSet and format patient names.

# load Ensembl reference
load("./Data/Ensembl.v75.annotation.RData") #for hg17 (Ensembl v75)

gep <- read_and_format_data("./Data/GEP_matrix_69_TPLL_vs_10_healthy_CD3_pan_T.csv")
rnaseq <- read_and_format_data("./Data/RNA_sequencing_FPKM_values_48_TPLL_and_2_FU_samples_vs_6_healthy_CD3_pan_T.csv")
sm_rnaseq <- read_and_format_data("./Data/small_RNA_sequencing_CPM_values_46_TPLL_vs_6_healthy_CD3_pan_T.csv")
snp <- read_and_format_data("./Data/sCNA_of_82_TPLL_and_6_FU_samples.csv")
wgs_mut <- read_and_format_data("./Data/WGS_mutations_3_TPLL_tg_pairs.csv")
wes_mut_paired <- read_and_format_data("./Data/WES_mutations_17_TPLL_tg_pairs.csv")
wes_mut_single <- read_and_format_data("./Data/WES_mutations_32_TPLL_tumor_single.csv")
wes_mut_single_followup <- read_and_format_data("./Data/WES_mutations_5_patients_with_FU.csv")

# for gene expression profiles, average out the gene expression values of duplicated ensembl id, and set ensembl gene id as rown name.
gep <- gep[!is.na(gep$Ensembl_ID),] # remove rows where Ensembl_ID is NA
gep_avg <- get_avg_gep(gep)
# gep_avg <- readRDS('./Data/gep_avg.rds')

# remove unnecessary columns and assign gene names as row names for SNP, WES and WGS data.

# RNA Seq
rownames(rnaseq) <- rnaseq$...1
rnaseq <- rnaseq[-c(1:5)]

# Sm-RNA Seq
rownames(sm_rnaseq) <- sm_rnaseq$...1
sm_rnaseq <- sm_rnaseq[-c(1:9)]

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
wes_mut_single_merged <- wes_mut_single_merged %>% rename(TP079 = TP029)

# WGS Mutation 
rownames(wgs_mut) <- wgs_mut$gene
wgs_mut <- wgs_mut[-c(1)]

# replace patient colnames with patient sample ids
gep_avg <- assign_p_number(gep_avg, 'gep')
rnaseq <- assign_p_number(rnaseq, 'mRNAseq')
sm_rnaseq <- assign_p_number(sm_rnaseq, 'miRNA')
snp <- assign_p_number(snp, 'snp')
wes_mut_paired <- assign_p_number(wes_mut_paired, 'wes_paired')
wes_mut_single_merged <- assign_p_number(wes_mut_single_merged, 'wes_single')
wgs_mut <- assign_p_number(wgs_mut, 'wgs')

# Order the samples by column name
gep_avg <- gep_avg[,order(colnames(gep_avg))]
rnaseq <- rnaseq[,order(colnames(rnaseq))]
sm_rnaseq <- sm_rnaseq[,order(colnames(sm_rnaseq))]
snp <- snp[,order(colnames(snp))]
wes_mut_paired <- wes_mut_paired[,order(colnames(wes_mut_paired))]
wes_mut_single_merged <- wes_mut_single_merged[,order(colnames(wes_mut_single_merged))]
wgs_mut <- wgs_mut[,order(colnames(wgs_mut))]

# Prepare clinical_sample dataframe to be used as colData in each molecular profile SummarizedExperiment
# Concatenate all sample names from molecular profiles + drug sensitivity data into one vector
total_samples <- unique(c(
  colnames(gep_avg), 
  colnames(rnaseq),
  colnames(sm_rnaseq),
  colnames(snp), 
  colnames(wes_mut_paired), 
  colnames(wes_mut_single_merged), 
  colnames(wgs_mut)
))
total_samples <- sort(total_samples)

clinical_sample_data <- get_clinical_sample_data(
  total_samples, 
  './Data/clinical_data_sample_based_2020_12_07.csv',
  './Data/p_number_moldata.xlsx',
  'wes_single_followup'
)

# remove columns that are not publicly shared - USED ONLY WHEN SUBMITTED TO ORCESTRA
# clinical_sample_data <- clinical_sample_data %>% select(p_number, cellid, time_series)

################################################################
##### CREATE SummarizedExperiments FOR MOLECULAR PROFILES ######
################################################################

# Create Gene Expression Data SummarizedExperiment
TPLL_GEP <- get_summarized_experiment(gep_avg, clinical_sample_data, features_gene, 'gene_id', 'rnaseq', datatype='gep')

# Create mRNASeq data SummarizedExperiment
TPLL_mRNA <- get_summarized_experiment(rnaseq, clinical_sample_data, features_gene, 'gene_id', 'rnaseq', datatype='mrna')

# Create miRNASeq data SummarizedExperiment
TPLL_miRNA <- get_summarized_experiment(sm_rnaseq, clinical_sample_data, features_gene, 'id', 'rnaseq', datatype='mirna')

# Create SNP Data SummarizedExperiment
TPLL_SNP <- get_summarized_experiment(snp, clinical_sample_data, features_gene, 'gene_name', 'mutation', datatype='snp')

# Create WGS (Whole Genome Sequencing) mutation data SummarizedExperiment
wgs_mut[wgs_mut == 'unmutated'] <- 'wt'
TPLL_WGS_MUT <- get_summarized_experiment(wgs_mut, clinical_sample_data, features_gene, 'gene_name', 'mutation', datatype='wgs')

# Create WES (Whole Exome Sequencing) single-end and paired-end mutation data SummarizedExperiments
### Single-end Mutation Data ###
wes_mut_single_merged[wes_mut_single_merged == 'unmutated'] <- 'wt'
TPLL_WES_MUT_SINGLE <- get_summarized_experiment(wes_mut_single_merged, clinical_sample_data, features_gene, 'gene_name', 'mutation', datatype='wes')

# Paied-end Mutation Data
wes_mut_paired[wes_mut_paired == 'unmutated'] <- 'wt'
TPLL_WES_MUT_PAIRED <- get_summarized_experiment(wes_mut_paired, clinical_sample_data, features_gene, 'gene_name', 'mutation', datatype='wes')

saveRDS(TPLL_GEP, './Data/curated_data/TPLL_GEP.rds')
saveRDS(TPLL_mRNA, './Data/curated_data/TPLL_mRNA.rds')
saveRDS(TPLL_miRNA, './Data/curated_data/TPLL_miRNA.rds')
saveRDS(TPLL_SNP, './Data/curated_data/TPLL_SNP.rds')
saveRDS(TPLL_WGS_MUT, './Data/curated_data/TPLL_WGS_MUT.rds')
saveRDS(TPLL_WES_MUT_SINGLE, './Data/curated_data/TPLL_WES_MUT_SINGLE.rds')
saveRDS(TPLL_WES_MUT_PAIRED, './Data/curated_data/TPLL_WES_MUT_PAIRED.rds')
saveRDS(clinical_sample_data, './Data/curated_data/clinical_sample_data.rds')