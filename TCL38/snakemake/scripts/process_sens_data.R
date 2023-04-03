library(stringr)
library(readxl)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

work_dir <- "/Users/minoru/Code/bhklab/JAKSTAT/Data"

files <- list.files(file.path(work_dir, "sensitivity"))
dfs <- list()
for (file in files) {
  if (str_detect(file, "CTG")) {
    df <- read_excel(file.path(work_dir, paste0("sensitivity/", file)))
    df <- df[!df$ProductId %in% c("empty", "BzCl", "DMSO"), ]
    name <- str_replace(file, "Raw_Data_", "")
    name <- str_replace(name, ".xlsx", "")
    dfs[[name]] <- data.frame(df)
  }
}

all_drugs <- c()
all_cells <- c()
merged_df <- data.frame(matrix(ncol=length(colnames(dfs[[1]]))))
colnames(merged_df) <- colnames(dfs[[1]])
for (name in names(dfs)) {
  merged_df <- rbind(merged_df, dfs[[name]])
}
merged_df <- merged_df[!is.na(merged_df$DWell), ]
all_drugs <- sort(unique(merged_df$ProductId))
all_cells <- sort(unique(merged_df$screen_id))

cols <- c('cell', 'drug', "Dose1", "Dose2", "Dose3", "Dose4", "Dose5")
dose_df <- data.frame(matrix(ncol=length(cols), nrow=0))
viability_df <- data.frame(matrix(ncol=length(cols), nrow=0))
colnames(dose_df) <- cols
colnames(viability_df) <- cols

for(cell in all_cells){
  for(drug in all_drugs){
    subset <- merged_df[merged_df$screen_id == cell & merged_df$ProductId == drug, ]
    subset <- subset[order(subset$Concentration), ]
    
    added_dose <- setNames(as.list(c(cell, drug, subset$Concentration)), cols)
    dose_df <- rbind(dose_df, added_dose)
    
    added_viability <- setNames(as.list(c(cell, drug, subset$inhibition_percent)), cols) 
    viability_df <- rbind(viability_df, added_viability)
  }
}

# Convert screen_id into standardized cell names
multiassay <- readRDS(file.path(work_dir, 'MultiAssayExp/TCL38_MultiAssayExp.rds'))
cells_standardized <- sort(rownames(data.frame(colData(multiassay))))

modified <- str_to_lower(all_cells)
modified <- str_replace_all(modified, '[\\W_]', '')
cell_map<- as.data.frame(list(screenid=all_cells, modified=modified))
cell_map$standard_name <- NA
for(cell in cells_standardized){
  comp <- str_to_lower(cell)
  comp <- str_replace_all(comp, '[\\W_]', '')
  found <- modified[str_detect(modified, comp)]
  print(found)
  if(length(found) > 0){
    cell_map[cell_map$modified %in% found, ]$standard_name <- cell
  }
}
# manually add missing standard cell name and make modifications:
cell_map[cell_map$screenid == 'SUH-DHL1', ]$standard_name <- 'Su_DHL_1'
cell_map[cell_map$screenid == 'CTG_Phera_NK_YS_w', ]$standard_name <- 'NK_YS_w'
cell_map[cell_map$screenid == 'CTG_Phera_NK_YS_wo', ]$standard_name <- 'NK_YS_wo'
cell_map[cell_map$screenid == 'CTG_Phera_TLBR2_w', ]$standard_name <- 'TLBR2_w'
cell_map[cell_map$screenid == 'CTG_Phera_TLBR2_wo', ]$standard_name <- 'TLBR2_wo'

dose_df$cell <- unlist(lapply(dose_df$cell, function(cell){
  cell_map[cell_map$screenid == cell, c('standard_name')]
}))
viability_df$cell <- unlist(lapply(viability_df$cell, function(cell){
  cell_map[cell_map$screenid == cell, c('standard_name')]
}))

dose_df$drug <- str_replace_all(dose_df$drug, '\\s', '_')
viability_df$drug <- str_replace_all(viability_df$drug, '\\s', '_')

