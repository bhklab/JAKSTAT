##### Using TPLL MultiassayExperiment object ########
# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html

# load PSet
TPLL_MultiAssayExp <- readRDS('./PSet/TPLL_MultiAssay.RDS')

# list of assay data in the object
TPLL_MultiAssayExp@ExperimentList

# list of sample and patient Ids (p number and TP number)
TPLL_MultiAssayExp@colData

# SampleMap: Identify which sample exists in which assay
smapleMap <- data.frame(sampleMap(TPLL_MultiAssayExp))

# Get experiments in the MultiAssayExperiment object
experiments <- experiments(TPLL_MultiAssayExp)

# Subset the object by assays
subset_data <- TPLL_MultiAssayExp[, , "rnaseq"]

# Subset the object by column (sample ids)
subset_data <- subsetByColumn(TPLL_MultiAssayExp, c("p1095", "p1376"))

# Subset the object by row (gene ids) ex. TP53
subset_data <- subsetByRow(TPLL_MultiAssayExp, c("ENSG00000141510"))


