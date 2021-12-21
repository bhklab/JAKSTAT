#' A series of script examples for the TCL38 MultiAssayExperiment object use case.
#' Documentation can be found at https://www.bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html
library(MultiAssayExperiment)

#' Read in the data object.
tcl38_multiassay <- readRDS("./data/TCL38_MultiAssayExp.rds")

#' 1. View the list of Experiments
#' The data object contains three data slots in `tcl38_multiassay@ExperimentList`:
#' 1. expr_counts: RNA-seq counts data.
#' 2. expr_tpm: RNA-seq tpm data.
#' 3: acgh: aCGH data
#' Summary of each data slots can be viewed with:
experiments(tcl38_multiassay)


#' 2. View the Cell Line Data
#' Information about each cell line used in the experiments are in `tcl38_multiassay@colData`.
#' They can be viewed with: 
colData(tcl38_multiassay)
#' Or they can be viewed in a data frame:
cells <- as.data.frame(colData(tcl38_multiassay))
#' There were differences in cell names between the ACGH data (the acgh_cells column) and the RNA-seq data (the rnaseq_cells column). 
#' They were standardized to the cell names assigned as row names in each row of the cells data frame.
#' For exmaple, "SMZ_1" referes to "SMZ1" in aCGH data and "SMZ-1" in the RNA-seq data.


#' 3. View the Mapping between the Cell Lines and Experiments
#' `tcl38_multiassay@sampleMap` contains the mapping table between the standardized cell names and sample names in each data slots.
#' They can be accessed with:
sampleMap(tcl38_multiassay)
#' Or viewed in a data frame:
sample_map <- as.data.frame(sampleMap(tcl38_multiassay))
#' `colname` refers to the sample name in each assay.
#' `primary` refers to the standardized name of the cell line that corresponds to the sample.
#' `assay` refers to the name of the data slot in the `ExperimentList`.


#' 4. Extracting the Whole Data
#' To get all RNA-seq counts data:
tcl38_multiassay[, , "expr_counts"]
#' `assays()` or `assay()` function can be used to extract the actual data:
counts <- assay(tcl38_multiassay[, , "expr_counts"])

#' Sample specific metadata for the RNA-seq data can be obtained with:
counts_rowdata <- tcl38_multiassay@ExperimentList@listData[["expr_counts"]]@phenoData@data
#' The sample specific data contains information such as which cell line the sample belongs to
#' and whether the sample was a re-run and/or healthy cell.

#' The aCGH data are stored in RaggedExperiment object.
#' To extract the whole aCGH data:
acgh <- as.data.frame(tcl38_multiassay@ExperimentList@listData[["acgh"]]@assays@unlistData)
#' The rownames are formatted as "cellname_abberationNum". 
#' "strand" column is set to "*" for all rows.


#' 5. Subsetting the Data
#' The data can be subset in the following format: tcl38_multiassay[rownames, colnames, assay].
#' rownames: transcript ID(s) for the RNA-seq data, "cellname_abberationNum" for the aCGH data.
#' colnames: the standardized cell line name(s) from `tcl38_multiassay@colData`.
#' assay: the experiment name(s).
#' For example, to get all RNA-seq counts data for the cell lines "ALL_SIL" and "DEL":
subset <- tcl38_multiassay[, c("ALL_SIL", "DEL") , "expr_counts"]
#' The subset object can be used to perform the operations as described in previous sections.
