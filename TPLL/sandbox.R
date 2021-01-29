### place to test out code before changing other files ###

col_names <- c(
  'p_number_gep', 
  'p_number_miRNA', 
  'p_number_mRNASeq', 
  'p_number_snp', 
  'p_number_wes_paired',
  'p_number_wes_single',
  'p_number_wes_single_followup',
  'p_number_wgs'
)
clinical_data[, col_names] <- NA

clin_rows <- rownames(clinical_data)

annotation <- read.csv(paste0('./Data/p-numbers/', 'p_number_gep', '.csv'))
row_names <- annotation$TP_number
str_replace(row_names, '^TP\\d{3}_', format_name)

format_name <- function(match){
  return(str_replace(match, '_', '.'))
}

for(col_name in col_names){
  annotation <- read.csv(paste0('./Data/p-numbers/', col_name, '.csv'))
  row_names <- annotation$TP_number
  row_names <- str_replace(row_names, '^TP\\d{3}_', format_name)
  rownames(annotation) <- row_names
  annotation_rows <- rownames(annotation)
  for(clin_row_name in clin_rows){
    if(clin_row_name %in% annotation_rows){
      clinical_data[clin_row_name, col_name] <- if(annotation[clin_row_name, 'p_number'] != '') annotation[clin_row_name, 'p_number'] else NA
    }
  }
}

clinical_data <- add_p_numbers_to_clinical_data(clinical_data, './Data/p-numbers/')
