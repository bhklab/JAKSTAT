library(stringr)
library(dplyr)
library(readxl)
library(PharmacoGx)

get_sample_raw_sensitivity <- function(sample_name, corrected_name, datatype, path){
  sample_data <- read_excel(path, sheet=datatype)
  sample_data <- data.frame(sample_data[, c('DRUG_NAME', 'Max.Conc.tested', 'Min.Conc.tested', 'D1', 'D2', 'D3', 'D4', 'D5')])
  row_names <- paste(corrected_name, sample_data$DRUG_NAME, sep='_')
  doses <- c('D1', 'D2', 'D3', 'D4', 'D5')
  
  dose_df <- sample_data
  num <- 0
  for(d in doses){
    dose_df[, d] <- (dose_df$Min.Conc.tested * (10^num))/1000
    num <- num + 1
  }
  dose_df <- dose_df[, doses]
  rownames(dose_df) <- row_names
  
  viability_df <- sample_data
  viability_df <- viability_df[, doses]
  rownames(viability_df) <- row_names
  viability_df <- 100 - viability_df
  return(list('dose'=dose_df, 'viability'=viability_df))
}

get_sensitivity_info <- function(raw.sensitivity, min_col, max_col) {
  sensitivity_info <- data.frame(matrix(data=NA, ncol=4, nrow=0))
  colnames(sensitivity_info) <- c('cellid', 'drugid', 'chosen.min.range', 'chosen.max.range')
  doseDF <- raw.sensitivity[,,'Dose']
  for(row in rownames(doseDF)){
    split <- strsplit(row, '_')
    cell <- if(stringr::str_detect(split[[1]][1], '^CD*') || stringr::str_detect(split[[1]][1], '^Helsinki*')) paste(split[[1]][1], split[[1]][2], sep='_') else split[[1]][1]
    drug <- if(stringr::str_detect(split[[1]][1], '^CD*') || stringr::str_detect(split[[1]][1], '^Helsinki*')) split[[1]][3] else split[[1]][2]
    sensitivity_info[row, ] <- c(cell, drug, doseDF[row, min_col], doseDF[row, max_col], FALSE)
  }
  return(sensitivity_info[order(rownames(sensitivity_info)),])
}

get_sensitivity_profile <- function(raw.sensitivity) {
  sensitivity_profile <- data.frame(matrix(data=NA, ncol=5, nrow=0))
  colnames(sensitivity_profile) <- c(
    'aac_recomputed', 
    'ic50_recomputed', 
    'HS', 
    'E_inf', 
    'EC50'
  )
  calculated_profiles <- PharmacoGx:::.calculateFromRaw(raw.sensitivity)
  
  for(row in rownames(raw.sensitivity[,,'Dose'])){
    sensitivity_profile[row, ] <- c(
      calculated_profiles$AUC[row], 
      calculated_profiles$IC50[row],
      calculated_profiles$pars[[row]]$HS,
      calculated_profiles$pars[[row]]$E_inf,
      calculated_profiles$pars[[row]]$EC50
    )
  }
  return(sensitivity_profile[order(rownames(sensitivity_profile)),])
}

tcl38_cells <- read.csv('./data/rnaseq_samples.csv')
tcl38_cells <- unique(tcl38_cells$cell)

filepath <- './data/41467_2018_3987_MOESM10_ESM.xlsx'

cell_names <- excel_sheets(filepath)

cells_map_df <- data.frame(matrix(ncol=0, nrow=(length(tcl38_cells))))
cells_map_df['cells'] <- tcl38_cells
cells_map_df['dufva_cells'] <- unlist(lapply(cells_map_df$cells, function(tcl38_cell){
  index <- match(str_replace(tolower(tcl38_cell), '-', ''), str_replace(tolower(cell_names), '-', ''))
  if(is.na(index)) NA else cell_names[index]
}))
cells_map_df <- cells_map_df[!is.na(cells_map_df$dufva_cells), ]

doses <- c('D1', 'D2', 'D3', 'D4', 'D5')
dose_df <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
viability_df <- data.frame(matrix(data=NA, ncol=length(doses), nrow=0))
colnames(dose_df) <- doses 
colnames(viability_df) <- doses

for(cell_name in cells_map_df$dufva_cells){
  raw_sens <- get_sample_raw_sensitivity(cell_name, cells_map_df[cells_map_df$dufva_cells == cell_name, ]$cells[1], cell_name, filepath)
  dose_df <- rbind(dose_df, raw_sens[['dose']])
  viability_df <- rbind(viability_df, raw_sens[['viability']])
}
colnames(dose_df) <- c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5') 
colnames(viability_df) <- c('Dose1', 'Dose2', 'Dose3', 'Dose4', 'Dose5') 
dose_df <- dose_df[order(rownames(dose_df)), ]
viability_df <- viability_df[order(rownames(viability_df)), ]

conc_tested <- 5

# Make the merged dataframe into an array
raw.sensitivity <- array(
  c(
    as.matrix(dose_df), 
    as.matrix(viability_df)
  ), 
  c(nrow(dose_df), conc_tested , 2), 
  dimnames=list(
    rownames(dose_df), 
    colnames(dose_df), 
    c("Dose", "Viability")
  )
)

sensitivity_info <- get_sensitivity_info(raw.sensitivity)
sensitivity_profile <- get_sensitivity_profile(raw.sensitivity)

# curate cell info
curationCell <- data.frame(matrix(ncol=1, nrow=(length(rownames(cells_map_df)))))
colnames(curationCell) <- c("unique.cellid")
curationCell$unique.cellid <- cells_map_df$cells
rownames(curationCell) <- curationCell$unique.cellid

# curate tissue info
curationTissue <- data.frame(matrix(ncol=0, nrow=(length(rownames(cells_map_df))))) 
curationTissue$cellid <- curationCell$unique.cellid
curationTissue$unique.tissueid <- 'Tissue'
curationTissue$tissue_type <- 'Tissue'
rownames(curationTissue) <- curationCell$unique.cellid

cell <- data.frame(matrix(ncol=0, nrow=(length(rownames(cells_map_df)))))
cell$cellid <- curationCell$unique.cellid
cell$tissueid <- curationTissue$unique.tissueid
cell$tissue_type <- curationTissue$tissue_type
rownames(cell) <- curationCell$unique.cellid

# extract drug info
drug <- read_excel(filepath, sheet='DSS_all')
drug <- drug[, c('DRUG_NAME', 'Mechanism.Targets', 'Class.explained')]
colnames(drug) <- c('drugid', 'mechanism_targets', 'class_explained')

curationDrug <- data.frame(matrix(ncol=0, nrow=(length(drug$drugid))))
curationDrug$unique.drugid <- sort(drug$drugid)
rownames(curationDrug) <- curationDrug$unique.drugid

PSet <- PharmacoGx::PharmacoSet(
  name="Dufva",
  cell=cell,
  drug=drug,
  sensitivityInfo=sensitivity_info,
  sensitivityRaw=raw.sensitivity,
  sensitivityProfiles=sensitivity_profile,
  sensitivityN=as.numeric(length(colnames(raw.sensitivity[,,'Dose']))),
  curationCell=curationCell,
  curationDrug=curationDrug,
  curationTissue=curationTissue,
  datasetType=c("both")
)

PSet@annotation$version <- 1	

saveRDS(PSet, file="./data/Dufva_PSet.rds")
