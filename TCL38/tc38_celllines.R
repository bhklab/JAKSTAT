setwd("~/Documents/github/jakstat-pset/TCL38")
library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
library(PharmacoGx)

tcl_cells <- data.frame(read.delim('cells.txt'))
tcl_cells$Cell.line <- gsub("[*]|[[:space:]+]", "", tcl_cells$Cell.line)
pharmacodb_cells <- read.csv("cells_pharmacodb.csv")
filtered <- pharmacodb_cells[tolower(pharmacodb_cells$cell.name) %in% tolower(tcl_cells$Cell.line) | tolower(pharmacodb_cells$synonym) %in% tolower(tcl_cells$Cell.line), ]
filtered$tcl.cell.name <- ""
filtered$dataset.cell.name <- ""
filtered <- filtered[, c("tcl.cell.name", "cell.name", "synonym", "dataset.cell.name", "dataset")]
for(cell in tcl_cells$Cell.line){
  filtered$tcl.cell.name[tolower(filtered$cell.name) == tolower(cell) | tolower(filtered$synonym) == tolower(cell)] <- cell
}
filtered$dataset.cell.name <- ifelse(filtered$synonym == "", filtered$cell.name, filtered$synonym)

df_list <- split(filtered, f = filtered$dataset)

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
  PharmacoGx::downloadPSet(name=pset_name, saveDir="./data")
}

psets <- unlist(psets)

cells <- unique(filtered$name)
dir.create('./out')

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
