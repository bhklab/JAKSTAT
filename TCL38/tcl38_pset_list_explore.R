library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
library(rlist)
library(PharmacoGx)

source(file='tcl38_functions.R') # Functions used to help extract information from the list of PSets

# list of PSets that are subset to contain all TCL38 cell lines data that could be found. (data pertaining to 16 of the 38 TCL38 cell lines were found in 5 PSets).
pset.list <- readRDS("./tcl38_pset_list.rds")

# Obtain the cell-drug combinations tested in each PSet.
cell_drug_combinations <- get_cell_drug_combinations(pset.list)

# Specify the cellid and drugid of interest, obtained from the "cell_drug_combinations" dataframe.
cellid <- "ALL-SIL"
drugid <- "Nilotinib"

# Filter the pset.list, only to keep the PSets that have the drug sensitivity data for the given cell line and drug combination.
cell_drug_subset <- cell_drug_combinations[cell_drug_combinations$cellid == cellid & cell_drug_combinations$drugid == drugid, ]
datasets <- names(cell_drug_subset)[which(cell_drug_subset == "Yes", arr.ind=T)[, "col"]]
filtered.pset.list <- pset.list[datasets]

# You can obtain a merged sensitivity data from all the PSets in the filtered list for the given combination of cellid and drugid
sensitivity_data_subset <- subset_sensitivity_data(filtered.pset.list, cellid, drugid)

# You can plot a dose response curve for the given combination of cellid and drugid.
PharmacoGx::drugDoseResponseCurve(drug=drugid, cellline=cellid, pSets=filtered.pset.list)