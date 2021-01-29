library(PharmacoGx)
library(plotly)

#########################
##### Using PSet ########
#########################

# load PSet
PSet <- readRDS('TPLL_PSet.RDS')

# 1. get summary of molecular profiles:
molecularProfiles  <- PSet@molecularProfiles

# 2. get summary of each molecular profile:
## for GEP
summarizeMolecularProfiles(PSet, mData='rnaseq')

## for SNP, WES and WGS mutation data (change 'mData' parameter to corresponding profile name)
summarizeMolecularProfiles(PSet, mData='snp', summary.stat='and')

# 3. summary of a molecular profile with a specific sample and  specific gene
# 'cell.lines' accespts one or more sample names
# for GEP, 'features' parameter accepts Ensembl ID(s)

gep_samples <- data.frame(molecularProfiles[["rnaseq"]]@colData) # sample data included in the rnaseq molecular profile.

gep_summary <- summarizeMolecularProfiles(PSet, mData='rnaseq', cell.lines = c('TP001', 'TP002', 'TP003', 'TP004'), features='ENSG00000223972')

gep_summary_table <- gep_summary@assays@data$exprs # displays the summary

# for SNP, WES and WGS mutation data, 'features' parameter accepts gene names
wes_single_summary <- summarizeMolecularProfiles(PSet, mDataType='mut_wes_single', summary.stat='or', cell.lines=c('TP014'), features=c('A1CF', 'ABCB6'))

wes_single_summary@assays@data$exprs

# to summarize the time series data,
# 1. obtain samples that are in the time series
wes_mut_single_samples <- data.frame(PSet@molecularProfiles$mut_wes_single@colData)

timeseries_samples <- wes_mut_single_samples[!is.na(wes_mut_single_samples$timeseries_wes_mut_single), c('cellid', 'timeseries_wes_mut_single')]

timeseries_summary <- summarizeMolecularProfiles(
  PSet, 
  mData='mut_wes_single', 
  summary.stat='or', 
  cell.lines=rownames(timeseries_samples[timeseries_samples$timeseries_wes_mut_single == 'TP092', ]), 
)

TP092_timeseries_summary <- timeseries_summary@assays@data$exprs


# 4. Drug Sensitivity Data Summary
# summmarize drug sensitivity data
drugs <- PSet@drug

sens_info <- PSet@sensitivity$info

sens_profile <- PSet@sensitivity$profiles

drug_sens_summary_ic50 <- summarizeSensitivityProfiles(PSet, cell.lines=unique(sens_info$cellid), sensitivity.measure="ic50_recomputed")

# dose response curve
PharmacoGx::drugDoseResponseCurve('Dinaciclib', 'TP111', list(PSet), plot.type = c("Both"))




# visualization of AAC in sample p1332, treated with different drugs
aac_p1332 <- drug_sens_summary_aac[,'p1332']
aac_p1332 <- sort(aac_p1332, decreasing=TRUE)

# create a barplot using the subset data.
table <- data.frame(drugs=names(aac_p1332), aac=(aac_p1332))
table$drugs <- factor(table$drugs, levels=names(aac_p1332))

fig <- plot_ly(
  data=table,
  x = ~drugs,
  y = ~aac,
  name = "AAC p1332",
  type = "bar"
)

fig <- fig %>% layout(
  title='AAC p1332',
  xaxis=list(title='Drugs'), 
  yaxis=list(title='AAC')
)

fig # displays the figure.