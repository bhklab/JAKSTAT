library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(data.table)

raw_data_dir <- './Data/additional_drug_screen_he_et_al/combination_screen/TPLL_data'
data_dir <- './Data/combinatory_drugs_data'
samples <- list.dirs(raw_data_dir, full.names = FALSE, recursive = FALSE)

for(sample in samples){
  files <- list.files(file.path(raw_data_dir, sample))
  for(file in files){
    raw_sheet <- read_excel(file.path(raw_data_dir, sample, file), sheet = 1)
    dir.create(file.path(data_dir, sample, 'csv'), recursive = TRUE, showWarnings = FALSE)
    write.csv(raw_sheet, file.path(data_dir, sample, 'csv', str_replace(file, '.xlsx', '.csv')), row.names=FALSE)
  }
}

df_doses_total <- data.frame(matrix(nrow=0, ncol=0))
for(sample in samples){
  print(sample)
  dir.create(file.path(data_dir, sample, 'rdata'), recursive = TRUE, showWarnings = FALSE)
  files <- list.files(file.path(data_dir, sample, 'csv'))
  for(file in files){
    print(file)
    lines <- readLines(file.path(data_dir, sample, 'csv', file), warn = FALSE)
    
    start <- str_which(lines, "Inhibition Data") + 2
    end <- start + 15
    joined <- paste(lines[start:end], collapse='\n')
    df_inhibition <- read.csv(text=joined, header=F)
    df_inhibition$V1 <- NULL
    colnames(df_inhibition) <- 1:24

    start <- str_which(lines, "Drug Pair") + 2
    end <- start + 15
    joined <- paste(lines[start:end], collapse='\n')
    df_drugpair <- read.csv(text=joined, header=F)
    df_drugpair$V1 <- NULL
    colnames(df_drugpair) <- 1:24

    drugpair_inhibition_dfs <- list()
    for(sect in 1: 6){
      row_drug1 <- row_drug2 <- col_drug1 <- col_drug2 <- 0
      if(sect %% 2 == 0){
        row_drug1 <- 10
        row_drug2 <- 9
        col_drug1 <- (sect / 2 - 1) * 8 + 1
        col_drug2 <- (sect / 2 - 1) * 8 + 2
      }else{
        row_drug1 <- 2
        row_drug2 <- 1
        col_drug1 <- (floor(sect / 2)) * 8 + 1
        col_drug2 <- (floor(sect / 2)) * 8 + 2
      }
      drug1=df_drugpair[row_drug1, col_drug1]
      drug2=df_drugpair[row_drug2, col_drug2]
      section_df <- df_inhibition[c(row_drug2:(row_drug2+7)), c(col_drug1:(col_drug1+7))]
      rownames(section_df) <- sprintf(paste(drug1, 'dose%d', sep='_'), 1:8)
      colnames(section_df) <- sprintf(paste(drug2, 'dose%d', sep='_'), 1:8)
      # drugpair_inhibition_dfs[[paste(drug1, drug2, sep='_')]] <- section_df
      saveRDS(section_df, file.path(data_dir, sample, 'rdata', paste0(sample, '_', drug1, '_', drug2, '.rds')))
    }
    
    start <- str_which(lines, "Summary") + 2
    end <- start + 5
    joined <- paste(lines[start:end], collapse='\n')
    df_doses <- read.csv(text=joined, header=F)
    df_doses <- df_doses[, c('V2', 'V3', 'V4', 'V5')]
    colnames(df_doses)[colnames(df_doses) %in% c('V2', 'V3')] <- c('drug1', 'drug2')
    df_doses <- separate(df_doses, col='V4', into=sprintf("drug1_dose%d", 1:8), sep=',')
    df_doses <- separate(df_doses, col='V5', into=sprintf("drug2_dose%d", 1:8), sep=',')
    df_doses <- add_column(df_doses, sample=sample, .before = "drug1")
    
    if(dim(df_doses_total)[1] == 0){
      df_doses_total <- data.frame(matrix(nrow=0, ncol=length(names(df_doses))))
      colnames(df_doses_total) <- names(df_doses)
      df_doses_total <- rbind(df_doses_total, df_doses)
    }else{
      df_doses_total <- rbind(df_doses_total, df_doses)
    }
  }
}
