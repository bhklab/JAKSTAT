setwd("~/Documents/github/jakstat-pset/TCL38")
library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
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
    dataset_cells <- intersect(df_list[dataset]$cell.name, pset_cells)
    subset <- PharmacoGx::subsetTo(pset, cells=dataset_cells)
    PSet_Subsets[[dataset]] <- c(PSet_Subsets[[dataset]], subset)
    rm(pset, pset_cells, dataset_cells, subset)
  }
}
# saveRDS(PSet_Subsets, "./data/PSet_Subsets.rds")
PSet_Subsets <- readRDS("./data/PSet_Subsets.rds")

##### Create merged cell and drug dataframes. #####
pset_cell <- data.frame(matrix(ncol=5, nrow=0))
colnames(pset_cell) <- c('cellid', 'original.cellid', 'dataset', 'pset.name', 'tissueid')

pset_drug <- data.frame(matrix(ncol=3, nrow=0))
colnames(pset_drug) <- c("drugid", "cid", "FDA")

for(dataset.name in names(PSet_Subsets)){
  tmp_cell <- PSet_Subsets[[dataset.name]][[1]]@cell[c("cellid", "tissueid")]
  tmp_cell$original.cellid <- tmp_cell$cellid
  tmp_cell$dataset <- dataset.name
  tmp_cell$pset.name <- get_pset_name(dataset.name, available)
  tmp_cell$cellid <- paste0(tmp_cell$cellid, '-', dataset.name) 
  row.names(tmp_cell) <- tmp_cell$cellid
  pset_cell <- rbind(pset_cell, tmp_cell)
  
  tmp_drug <- PSet_Subsets[[dataset.name]][[1]]@drug[c("drugid", "cid", "FDA")]
  rownames(tmp_drug) <- NULL
  pset_drug <- rbind(pset_drug, tmp_drug)
}

pset_drug <- pset_drug[!duplicated(pset_drug[, c("drugid", "cid")]),]
rownames(pset_drug) <- pset_drug$drugid

##### Merge Molecular Profiles #####
molecular_prof_names <- c()
for(dataset.name in names(PSet_Subsets)){
  molecular_prof_names <- c(molecular_prof_names, names(PSet_Subsets[[dataset.name]][[1]]@molecularProfiles))
}
molecular_prof_names <- unique(molecular_prof_names)
molecular_profiles <- list()
for(mol.prof in molecular_prof_names){
  se <- merge_molecular_profile(PSet_Subsets, mol.prof)
  molecular_profiles[[mol.prof]] <- c(molecular_profiles[[mol.prof]], se)
}

##### Merge drug sensitivity data #####
sensitivity_info_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))
sensitivity_profile_total <- data.frame(matrix(data=NA, ncol=0, nrow=0))

for(dataset.name in names(PSet_Subsets)){
  sensitivity_info_total <- merge_sensitivity_data(
    PSet_Subsets[[dataset.name]][[1]]@sensitivity[["info"]],
    dataset.name,
    sensitivity_info_total,
    "info"
  )
  sensitivity_profile_total <- merge_sensitivity_data(
    PSet_Subsets[[dataset.name]][[1]]@sensitivity[["profiles"]],
    dataset.name,
    sensitivity_profile_total,
    "profiles"
  )
}



##### Create the meta PSet #####
PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=molecular_profiles,
  name="TCL38_Meta",
  cell=cell,
  drug=drug,
  sensitivityInfo=sensitivity_info,
  sensitivityRaw=raw.sensitivity,
  sensitivityProfiles=sensitivity_profile,
  sensitivityN=as.numeric(length(colnames(raw.sensitivity[,,'Dose']))),
  curationCell=curationCell,
  curationDrug=curationDrug,
  curationTissue=curationTissue,
  datasetType=c("both"))

# # Access PhamracoDB's GraphQL API
# link <- 'https://www.pharmacodb.ca/graphql'
# # link <- 'http://localhost:5000/graphql'
# conn <- ghql::GraphqlClient$new(url = link)
# query <- 'query getExperimentsQuery($cellLineName: String){
#   experiments(cellLineName: $cellLineName, all: true){
#     id
#     cell_line {
#       name
#     }
#     compound {
#       name
#     }
#     dataset {
#       name
#     }
#   }
# }'
# 
# for(cellLine in cells){
#   variable <- list(cellLineName = cellLine)
#   experimentsQuery <- Query$new()$query('query', query)
#   # Execute query
#   result <- conn$exec(experimentsQuery$query, variables = variable) %>% jsonlite::fromJSON(flatten = T)
# 
#   # Modify the returned data
#   experiments <- result$data$experiments
#   
#   profile <- experiments[, -which(names(experiments) == "dose_response")]
#   
#   experiments <- tidyr::unnest(experiments, dose_response)
#   experiments <- experiments[, c("id", "cell_line.name", "compound.name", "dataset.name", "dose", "response")]
#   colnames(experiments)[which(names(experiments) == "id")] <- "experiment.id"
# 
#   dir.create(str_glue("./out/{cellLine}"))
# 
#   # Output data by compound
#   compounds <- unique(experiments$compound.name)
#   for(compound in compounds){
#     prof_subset <- profile[profile$compound.name == compound, ]
#     exp_subset <- experiments[ experiments$compound.name == compound, ]
#     compound <- gsub("[:]|[/]", "_", compound)
#     dir.create(str_glue("./out/{cellLine}/{compound}"))
#     write.csv(prof_subset, file=str_glue("./out/{cellLine}/{compound}/profile_{cellLine}_{compound}.csv"), row.names = FALSE)
#     write.csv(exp_subset, file=str_glue("./out/{cellLine}/{compound}/dose_response_{cellLine}_{compound}.csv"), row.names = FALSE)
#   }
# 
#   links <- character(length(compounds))
#   for(i in 1:length(links)){
#     links[i] <- str_glue("https://www.pharmacodb.ca/search?compound={compounds[i]}&cell_line={cellLine}")
#   }
#   pharmacodb_links <- data.frame(
#     cell_line = cellLine,
#     compound = compounds,
#     pharmacodb_link = links
#   )
#   write.csv(pharmacodb_links, str_glue("./out/{cellLine}/_pharmacodb_links.csv"), row.names = FALSE)
# }
