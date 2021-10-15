setwd("~/Documents/github/jakstat-pset/TCL38")
library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)

tcl_cells <- data.frame(read.delim('./data/cells.txt'))
tcl_cells$cell.lines <- gsub("[*]|[[:space:]+]", "", tcl_cells$cell.lines)
pharmacodb_cells <- read.csv("./data/cells_pharmacodb.csv")
filtered <- pharmacodb_cells[tolower(pharmacodb_cells$cell_name) %in% tolower(tcl_cells$cell.lines), ]
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
    dose_response {
      dose
      response
    }
    profile {
      AAC
      IC50
      EC50
      Einf
    }
  }
}'

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
