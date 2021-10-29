setwd("~/Documents/github/jakstat-pset/TCL38")
library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
library(PharmacoGx)

get_pset_name <- function(dataset, available.psets){
  if(str_detect(dataset, "GDSC")){
    gdsc_ver <- regmatches(dataset, gregexpr("[[:digit:]]+", dataset))
    available_versions <- available.psets$version[available.psets$`Dataset Name` == 'GDSC']
    version <- grep(paste0("v", gdsc_ver), available_versions, value=TRUE)
    return(paste0("GDSC_", version))
  }else{
    return(available.psets$`PSet Name`[available.psets$`Dataset Name` == dataset])
  }
}

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

# down load PSets that contains matched cell lines
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
  PharmacoGx::downloadPSet(name=pset_name, saveDir="./data", timeout=2000) # downloading may generate errors if done all at once.
}

subset_list <- list()
for(dataset in names(df_list)){
  pset_name <- available$`PSet Name`[available$`Dataset Name` == dataset]
  if(pset_name != ""){
    pset <- readRDS(paste0('./data/', pset_name, '.rds'))
    pset_cells <- cellNames(pset)
    dataset_cells <- intersect(df_list[dataset]$cell.name, pset_cells)
    subset <- PharmacoGx::subsetTo(pset, cells=dataset_cells)
    subset_list[[dataset]] <- c(subset_list[[dataset]], subset)
  }
}

# saveRDS(subset_list, "./data/PSet_Subsets.rds")

PSet_Subsets <- readRDS("./data/PSet_Subsets.rds")

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

assay <- data.frame(PSet_Subsets[["CCLE"]][[1]]@molecularProfiles[["rna"]]@assays@data@listData[["exprs"]])
coldata <- data.frame(PSet_Subsets[["CCLE"]][[1]]@molecularProfiles[["rna"]]@colData)
rowData <- data.frame(PSet_Subsets[["CCLE"]][[1]]@molecularProfiles[["rna"]]@elementMetadata)

Sum_Exp <- SummarizedExperiment(
  assays=list(exprs=as.matrix(filtered_data)), #gene expression data in matrix
  rowData=filtered_gene_info, #gene features Data
  colData=filtered_clinical_sample_data #patient/sample info
) 
Sum_Exp@metadata$annotation <- annotation

#######################
##### CREATE PSET #####
#######################

PSet <- PharmacoGx::PharmacoSet(
  molecularProfiles=list(
    "gep"=TPLL_GEP, 
    "rnaseq"=TPLL_mRNA,
    "mi_rnaseq"=TPLL_miRNA,
    "snp"=TPLL_SNP,
    "mut_wes_single"=TPLL_WES_MUT_SINGLE, 
    "mut_wes_paired"=TPLL_WES_MUT_PAIRED,
    "mut_wgs"=TPLL_WGS_MUT),
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

# Access PhamracoDB's GraphQL API
link <- 'https://www.pharmacodb.ca/graphql'
# link <- 'http://localhost:5000/graphql'
conn <- ghql::GraphqlClient$new(url = link)
query <- 'query getExperimentsQuery($cellLineName: String){
  experiments(cellLineName: $cellLineName, all: true){
    id
    cell_line {
      name
    }
    compound {
      name
    }
    dataset {
      name
    }
  }
}'

# experiments(cellLineName: $cellLineName, all: true){
#   id
#   cell_line {
#     name
#   }
#   compound {
#     name
#   }
#   dataset {
#     name
#   }
#   dose_response {
#     dose
#     response
#   }
#   profile {
#     AAC
#     IC50
#     EC50
#     Einf
#   }
# }

for(cellLine in cells){
  variable <- list(cellLineName = cellLine)
  experimentsQuery <- Query$new()$query('query', query)
  # Execute query
  result <- conn$exec(experimentsQuery$query, variables = variable) %>% jsonlite::fromJSON(flatten = T)

  # Modify the returned data
  experiments <- result$data$experiments
  
  profile <- experiments[, -which(names(experiments) == "dose_response")]
  
  experiments <- tidyr::unnest(experiments, dose_response)
  experiments <- experiments[, c("id", "cell_line.name", "compound.name", "dataset.name", "dose", "response")]
  colnames(experiments)[which(names(experiments) == "id")] <- "experiment.id"

  dir.create(str_glue("./out/{cellLine}"))

  # Output data by compound
  compounds <- unique(experiments$compound.name)
  for(compound in compounds){
    prof_subset <- profile[profile$compound.name == compound, ]
    exp_subset <- experiments[ experiments$compound.name == compound, ]
    compound <- gsub("[:]|[/]", "_", compound)
    dir.create(str_glue("./out/{cellLine}/{compound}"))
    write.csv(prof_subset, file=str_glue("./out/{cellLine}/{compound}/profile_{cellLine}_{compound}.csv"), row.names = FALSE)
    write.csv(exp_subset, file=str_glue("./out/{cellLine}/{compound}/dose_response_{cellLine}_{compound}.csv"), row.names = FALSE)
  }

  links <- character(length(compounds))
  for(i in 1:length(links)){
    links[i] <- str_glue("https://www.pharmacodb.ca/search?compound={compounds[i]}&cell_line={cellLine}")
  }
  pharmacodb_links <- data.frame(
    cell_line = cellLine,
    compound = compounds,
    pharmacodb_link = links
  )
  write.csv(pharmacodb_links, str_glue("./out/{cellLine}/_pharmacodb_links.csv"), row.names = FALSE)
}
