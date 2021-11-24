library(dplyr) 
library(tidyr)
library(stringr)
library(PharmacoGx)
library(plotly)
library(ggplot2)
library(reshape2)

source(file='../tcl38_functions.R')

pset.list <- readRDS("../misc/data/tcl38_pset_list.rds")

cell_drug_combinations <- get_cell_drug_combinations(pset.list)

cells <- levels(cell_drug_combinations$cellid)
drugs <- levels(cell_drug_combinations$drugid)

current.compounds <- read.delim("current_compounds.txt", header=FALSE)
current.compounds <- current.compounds$V1
available.current.compounds <- drugs[tolower(drugs) %in% tolower(current.compounds)]

experiments.pharmacodb <- read.csv("experiments_pharmacodb.csv")
experiments.pharmacodb <- experiments.pharmacodb[experiments.pharmacodb$cell %in% cells, ]

heatmap.aac <- data.frame(matrix(data=NA, nrow=length(current.compounds), ncol=length(cells)))
colnames(heatmap.aac) <- cells
rownames(heatmap.aac) <- current.compounds

heatmap.stdev <- heatmap.aac

# cell.drug.datasets <- heatmap.aac
# for(cell in colnames(cell.drug.datasets)){
#   cell.drug.dataset <- lapply(rownames(cell.drug.datasets), function(compound){
#     subsetted = experiments.pharmacodb[experiments.pharmacodb$cell == cell & experiments.pharmacodb$compound == compound, ]
#     if(dim(subsetted)[1] > 0){
#       aac <- subsetted$AAC[subsetted$AAC != "NULL"]
#       if(length(aac) == 0){
#         return(NA)
#       }
#       datasets <- unique(subsetted$dataset)
#       return(paste(datasets, collapse=" "))
#     }
#     return(NA)
#   })
#   cell.drug.dataset <- unlist(cell.drug.dataset)
#   cell.drug.datasets[, cell] <- cell.drug.dataset
# }
# cell.drug.datasets <- melt(as.matrix(cell.drug.datasets))
# cell.drug.datasets <- cell.drug.datasets[!is.na(cell.drug.datasets$value), ]

for(cell in colnames(heatmap.aac)){
  mdn <- lapply(rownames(heatmap.aac), function(compound){
    subsetted = experiments.pharmacodb[experiments.pharmacodb$cell == cell & experiments.pharmacodb$compound == compound, ]
    if(dim(subsetted)[1] > 0){
      aac <- subsetted$AAC[subsetted$AAC != "NULL"]
      if(length(aac) == 0){
        return(NA)
      }
      aac <- as.numeric(aac)
      return(median(aac, na.rm=TRUE))
    }
    return(NA)
  })
  mdn <- unlist(mdn)
  heatmap.aac[, cell] <- mdn
  
  stdev <- lapply(rownames(heatmap.stdev), function(compound){
    subsetted = experiments.pharmacodb[experiments.pharmacodb$cell == cell & experiments.pharmacodb$compound == compound, ]
    if(dim(subsetted)[1] > 0){
      aac <- subsetted$AAC[subsetted$AAC != "NULL"]
      if(length(aac) == 0){
        return(NA)
      }
      aac <- as.numeric(aac)
      if(length(aac) == 1){
        return(0)
      }
      return(sd(aac, na.rm=TRUE))
    }
    return(NA)
  })
  stdev <- unlist(stdev)
  heatmap.stdev[, cell] <- stdev
}

aac.long <- melt(as.matrix(heatmap.aac))
colnames(aac.long) <- c("compound", "cell", "AAC")
ggplot(aac.long, aes(x = cell, y = compound, fill = AAC)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "yellow", high = "red", na.value="white") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

stdev.long <- melt(as.matrix(heatmap.stdev))
colnames(stdev.long) <- c("compound", "cell", "Std.Dev")
ggplot(stdev.long, aes(x = cell, y = compound, fill = Std.Dev)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "red", high = "yellow", na.value="white") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))