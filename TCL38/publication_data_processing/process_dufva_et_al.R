library(stringr)
library(dplyr)
library(readxl)

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

tcl38_cells <- read.csv('../create_multiassay/data/rnaseq_samples.csv')
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
drugs <- read_excel(filepath, sheet='DSS_all')
drugs <- drugs[, c('DRUG_NAME', 'Mechanism.Targets', 'Class.explained')]
colnames(drugs) <- c('drugid', 'mechanism_targets', 'class_explained')

curationDrug <- data.frame(matrix(ncol=0, nrow=(length(drugs$drugid))))
curationDrug$unique.drugid <- sort(drugs$drugid)
rownames(curationDrug) <- curationDrug$unique.drugid


