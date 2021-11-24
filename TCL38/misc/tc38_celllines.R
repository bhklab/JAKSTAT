setwd("~/Documents/github/jakstat-pset/TCL38/misc")
library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
library(rlist)
library(PharmacoGx)

source(file='functions.R') # import functions 

##### Identify datasets that contains selected cell lines #####
tcl_cells <- data.frame(read.delim('cells.txt'))
tcl_cells$Cell.line <- gsub("[*]|[[:space:]+]", "", tcl_cells$Cell.line)
pharmacodb_cells <- read.csv("cells_pharmacodb.csv")

filtered <- pharmacodb_cells[tolower(pharmacodb_cells$cell.name) %in% tolower(tcl_cells$Cell.line) | tolower(pharmacodb_cells$synonym) %in% tolower(tcl_cells$Cell.line), ]
filtered$tcl.cell.name <- ""
filtered <- filtered[, c("tcl.cell.name", "cell.name", "synonym", "dataset")]
for(cell in tcl_cells$Cell.line){
  filtered$tcl.cell.name[tolower(filtered$cell.name) == tolower(cell) | tolower(filtered$synonym) == tolower(cell)] <- cell
}
filtered <- filtered[filtered$synonym != "", ]

df_list <- split(filtered, f = filtered$dataset)

##### download PSets that contains matched cell lines #####
available <- PharmacoGx::availablePSets()
psets <- list()
for(dataset in names(df_list)){
  pset_name <- ""
  if(str_detect(dataset, "GDSC")){
    gdsc_ver <- regmatches(dataset, gregexpr("[[:digit:]]+", dataset))
    available_versions <- available$version[available$`Dataset Name` == 'GDSC']
    version <- grep(paste0("v", gdsc_ver), available_versions, value=TRUE)
    print(paste0("GDSC ", version))
    pset_name <- paste0("GDSC_", version)
  }else{
    print(dataset)
    pset_name <- available$`PSet Name`[available$`Dataset Name` == dataset]
  }
  psets <- append(psets, pset_name)
  # PharmacoGx::downloadPSet(name=pset_name, saveDir="./data", timeout=2000) # downloading may generate errors if done all at once.
}

##### Get a list of PSet subsets that contains the selected cell lines only. #####
PSet_Subsets <- list()
for(dataset in names(df_list)){
  pset_name <- get_pset_name(dataset, available)
  if(length(pset_name) > 0){
    pset <- readRDS(paste0('./data/', pset_name, '.rds'))
    pset_cells <- cellNames(pset)
    dataset_cells <- intersect(df_list[[dataset]][["cell.name"]], pset_cells)
    subset <- PharmacoGx::subsetTo(pset, cells=dataset_cells)
    PSet_Subsets[dataset] <- subset
    rm(pset, pset_cells, dataset_cells, subset)
  }
}


# dataset <- "GDSC2"
# pset_name <- get_pset_name(dataset, available)
# pset <- readRDS(paste0('./data/', pset_name, '.rds'))
# pset_cells <- cellNames(pset)
# dataset_cells <- intersect(df_list[[dataset]][["cell.name"]], pset_cells)
# subset <- PharmacoGx::subsetTo(pset, cells=dataset_cells)
# PSet_Subsets[dataset] <- subset
# rm(pset, pset_cells, dataset_cells, subset)

# saveRDS(PSet_Subsets, "./data/PSet_Subsets.rds")
PSet_Subsets <- readRDS("./data/PSet_Subsets.rds")

##### Create merged cell and drug dataframes. #####
pset_cell <- get_slot_data(PSet_Subsets, "cell")
pset_drug <- get_slot_data(PSet_Subsets, "drug")

##### Merge Molecular Profiles #####
molecular_prof_names <- c()
for(dataset.name in names(PSet_Subsets)){
  molecular_prof_names <- c(molecular_prof_names, names(PSet_Subsets[[dataset.name]]@molecularProfiles))
}
molecular_prof_names <- unique(molecular_prof_names)
molecular_profiles <- list()
for(mol.prof in molecular_prof_names){
  se <- merge_molecular_profile(PSet_Subsets, mol.prof)
  molecular_profiles[mol.prof] <- se
}

##### Merge drug sensitivity data #####
sensitivity_info_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
sensitivity_profile_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))

for(dataset.name in names(PSet_Subsets)){
  sensitivity_info_total <- merge_sensitivity_data(
    PSet_Subsets[[dataset.name]]@sensitivity[["info"]],
    dataset.name,
    sensitivity_info_total,
    "info"
  )
  sensitivity_profile_total <- merge_sensitivity_data(
    PSet_Subsets[[dataset.name]]@sensitivity[["profiles"]],
    dataset.name,
    sensitivity_profile_total,
    "profiles"
  )
}

# Standardize the dose/viability col names
dose_colnames_total <- c()
for(dataset.name in names(PSet_Subsets)){
  tmp <- data.frame(PSet_Subsets[[dataset.name]]@sensitivity[["raw"]])
  dose_colnames_total <- c(dose_colnames_total, colnames(tmp))
}
dose_colnames_total <- format_dose_col(dose_colnames_total)
dose_colnames_total <- unique(dose_colnames_total)

doses <- dose_colnames_total[grepl(".Dose", dose_colnames_total)]
doses <- str_sort(doses, numeric=TRUE)
viabilities <- dose_colnames_total[grepl(".Viability", dose_colnames_total)]
viabilities <- str_sort(viabilities, numeric=TRUE)
dose_colnames_total <- c(doses, viabilities)

# Merge raw sensitivity data
sensitivity_raw_total <- data.frame(matrix(data=NA, ncol=length(dose_colnames_total), nrow=0))
colnames(sensitivity_raw_total) <- dose_colnames_total
for(dataset.name in names(PSet_Subsets)){
  tmp <- data.frame(PSet_Subsets[[dataset.name]]@sensitivity[["raw"]])
  rownames(tmp) <- paste0(rownames(tmp), ".", dataset.name)
  colnames_tmp <- colnames(tmp)
  colnames(tmp) <- format_dose_col(colnames_tmp)
  tmp[setdiff(colnames(sensitivity_raw_total), colnames(tmp))] <- NA
  sensitivity_raw_total <- dplyr::bind_rows(sensitivity_raw_total, tmp)
}

##### Create curation objects #####
curationCell <- get_slot_data(PSet_Subsets, "cell", isCuration=TRUE)
curationDrug <- get_slot_data(PSet_Subsets, "drug", isCuration=TRUE)
curationTissue <- get_slot_data(PSet_Subsets, "tissue", isCuration=TRUE)

##### Create the meta PSet #####
PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=molecular_profiles,
  name="TCL38_Meta",
  cell=pset_cell,
  drug=pset_drug,
  sensitivityInfo=sensitivity_info_total,
  sensitivityRaw=sensitivity_raw_total,
  sensitivityProfiles=sensitivity_profile_total,
  sensitivityN=as.numeric(length(colnames(sensitivity_raw_total)[grepl(".Dose", colnames(sensitivity_raw_total))])),
  curationCell=curationCell,
  curationDrug=curationDrug,
  curationTissue=curationTissue,
  datasetType=c("both")
)

psetList <- list()
tmp <- intersectPSet(PSet_Subsets, intersectOn=c("drugs", "cell.lines"))
tmp <- list.remove(tmp, "gCSI")
drugDoseResponseCurve(drug="Crizotinib", cellline="HH [Human lymphoma]", pSets=list(PSet_Subsets[["CTRPv2"]]))
drugDoseResponseCurve(drug="Crizotinib", cellline="SU-DHL-1", pSets=tmp)

sensitivity_dataset <- expand.grid(cellid=pset_cell$cellid, drugid=pset_drug$drugid)
for(dataset.name in names(PSet_Subsets)){
  info <- PSet_Subsets[[dataset.name]]@sensitivity[["info"]]
  sensitivity_dataset[dataset.name] <- ifelse(
    is.na(match(
      paste0(sensitivity_dataset$cellid, sensitivity_dataset$drugid),
      paste0(info$cellid, info$drugid),
    )),
    NA, "Yes"
  )
}


cellid <- "ALL-SIL"
drugid <- "Nilotinib"

sensitivity_dataset_subset <- sensitivity_dataset[sensitivity_dataset$cellid == cellid & sensitivity_dataset$drugid == drugid, ]
datasets <- names(sensitivity_dataset_subset)[which(sensitivity_dataset_subset == "Yes", arr.ind=T)[, "col"]]

sens_list <- list()
for(dataset.name in datasets){
  info <- PSet_Subsets[[dataset.name]]@sensitivity[["info"]]
  info <- info[info$cellid == cellid & info$drugid == drugid, ]
  profiles <- PSet_Subsets[[dataset.name]]@sensitivity[["profiles"]]
  profiles <- profiles[rownames(profiles) %in% rownames(info), ]
  raw <- data.frame(PSet_Subsets[[dataset.name]]@sensitivity[["raw"]])
  raw <- raw[rownames(raw) %in% rownames(info), ]
  sens_list[[dataset.name]] <- list(info=info, profiles=profiles, raw=raw)
}

sensitivity_info_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
sensitivity_profile_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))

for(dataset.name in names(sens_list)){
  sensitivity_info_total <- merge_sensitivity_data(
    sens_list[[dataset.name]][["info"]],
    dataset.name,
    sensitivity_info_total,
    "info"
  )
  sensitivity_profile_total <- merge_sensitivity_data(
    sens_list[[dataset.name]][["profiles"]],
    dataset.name,
    sensitivity_profile_total,
    "profiles"
  )
}

# Standardize the dose/viability col names
dose_colnames_total <- c()
for(dataset.name in names(sens_list)){
  tmp <- sens_list[[dataset.name]][["raw"]]
  dose_colnames_total <- c(dose_colnames_total, colnames(tmp))
}
dose_colnames_total <- format_dose_col(dose_colnames_total)
dose_colnames_total <- unique(dose_colnames_total)

doses <- dose_colnames_total[grepl(".Dose", dose_colnames_total)]
doses <- str_sort(doses, numeric=TRUE)
viabilities <- dose_colnames_total[grepl(".Viability", dose_colnames_total)]
viabilities <- str_sort(viabilities, numeric=TRUE)
dose_colnames_total <- c(doses, viabilities)

# Merge raw sensitivity data
sensitivity_raw_total <- data.frame(matrix(data=NA, ncol=length(dose_colnames_total), nrow=0))
colnames(sensitivity_raw_total) <- dose_colnames_total
for(dataset.name in names(sens_list)){
  tmp <- sens_list[[dataset.name]][["raw"]]
  rownames(tmp) <- paste0(rownames(tmp), ".", dataset.name)
  colnames_tmp <- colnames(tmp)
  colnames(tmp) <- format_dose_col(colnames_tmp)
  tmp[setdiff(colnames(sensitivity_raw_total), colnames(tmp))] <- NA
  sensitivity_raw_total <- dplyr::bind_rows(sensitivity_raw_total, tmp)
}

subset_sensitivity <- list(info=sensitivity_info_total, profiles=sensitivity_profile_total, raw=sensitivity_raw_total)
