library(PharmacoGx)
library(plotly)

##### Using PSet ########

# load PSet
PSet <- readRDS('./TPLL_PSet.rds')

# 1. get list of molecular profiles:
molecularProfiles  <- molecularProfilesSlot(PSet)

# To access sample data for each molecular profile:
gep_samples <- data.frame(phenoInfo(PSet, 'gep')) # replace "gep" with other molecular data types

# To access patient data:
patient <- cellInfo(PSet)

# 2. get summary of each molecular profile:
## for GEP
summarizeMolecularProfiles(PSet, mData='gep')

## For SNP, WES and WGS mutation data (change 'mData' parameter to corresponding profile name)
summarizeMolecularProfiles(PSet, mData='mut_wgs', summary.stat='and')

# 3. summary of a molecular profile with a specific sample and  specific gene
# 'cell.lines' accespts one or more sample names
# for GEP, 'features' parameter accepts Ensembl ID(s)

gep_samples <- data.frame(phenoInfo(PSet, 'gep')) # sample data included in the gep molecular profile.

gep_summary <- summarizeMolecularProfiles(
  PSet, 
  mData='gep', 
  cell.lines = c('TP001', 'TP002', 'TP003', 'TP004'), 
  features='ENSG00000223972'
)

molecularProfiles()

# For SNP, WES and WGS mutation data, 'features' parameter accepts gene names
wes_single_summary <- summarizeMolecularProfiles(
  PSet, 
  mDataType='mut_wes_single', 
  summary.stat='or', 
  cell.lines=c('TP014'), 
  features=c('A1CF', 'ABCB6')
)

# To summarize the time series data,

# 1. obtain samples that are in the time series
wes_mut_single_samples <- data.frame(phenoInfo(PSet, 'mut_wes_single'))
timeseries_samples <- wes_mut_single_samples[!is.na(wes_mut_single_samples$time_series), c('cellid', 'time_series')]
timeseries_samples <- timeseries_samples[with(timeseries_samples, order(cellid, time_series)), ]

# 2. display the genes that are recorded in the timeseries samples
wes_single <- molecularProfiles(PSet, 'mut_wes_single')
timeseries_rows <- rownames(timeseries_samples[timeseries_samples$cellid == 'TP092', ])
wes_single <- wes_single[, timeseries_rows]
wes_single <- wes_single[rowSums(is.na(wes_single)) != ncol(wes_single), ]

# 3. Drug Sensitivity Data Summary
# summmarize drug sensitivity data
drugs <- drugInfo(PSet)

sens_info <- sensitivityInfo(PSet)

sens_profile <- sensitivityProfiles(PSet)

drug_sens_summary_ic50 <- summarizeSensitivityProfiles(PSet, cell.lines=unique(sens_info$cellid), sensitivity.measure="ic50_recomputed")

# dose response curve
drugDoseResponseCurve('Tipifarnib', 'TP050', list(PSet), plot.type = c("Both"))

