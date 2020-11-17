library(PharmacoGx)

#########################
##### Using PSet ########
#########################

PSet <- readRDS('TPLL_PSet.RDS')

# 1. get summary of molecular profiles:
PSet@molecularProfiles

# 2. get summary of each molecular profile:
## for GEP
summarizeMolecularProfiles(PSet, mData='rnaseq')

## for SNP, WES and WGS mutation data (change 'mData' parameter to corresponding profile name)
summarizeMolecularProfiles(PSet, mData='snp', summary.stat='and')

# 3. summary of a molecular profile with a specific sample and  specific gene
# 'cell.lines' accespts one or more sample names
# for GEP, 'features' parameter accepts Ensembl ID(s)
gep_summary <- summarizeMolecularProfiles(PSet, mData='rnaseq', cell.lines = 'TP047', features='ENSG00000223972')
gep_summary@assays@data$exprs # displays the summary

# for SNP, WES and WGS mutation data, 'features' parameter accepts gene names
snp_summary <- summarizeMolecularProfiles(PSet, mData='snp', summary.stat='and', cell.lines='TP047', features='ZXDB')
snp_summary@assays@data$exprs

wes_single_summary <- summarizeMolecularProfiles(PSet, mDataType='mut_wes_single', summary.stat='or', cell.lines=c('TP014'), features=c('A1CF', 'ABCB6'))
wes_single_summary@assays@data$exprs

# to summarize the time series data,
# 1. obtain samples that are in the time series
samples <- data.frame(PSet@molecularProfiles$mut_wes_single@colData)
timeseries_samples <- samples[!is.na(samples$timeseries_wes_mut_single), c('cellid', 'timeseries_wes_mut_single')]
timeseries_cells <- rownames(timeseries_samples[timeseries_samples$timeseries_wes_mut_single == 'TP092', ])
timeseries_summary <- summarizeMolecularProfiles(PSet, mData='mut_wes_single', summary.stat='and', cell.lines=timeseries_cells, features=c('ABCB5', 'ABCB6'))
timeseries_summary@assays@data$exprs