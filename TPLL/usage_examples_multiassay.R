library(MultiAssayExperiment)
library(S4Vectors)

##### Using TPLL MultiassayExperiment object ########
# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/QuickStartMultiAssay.html

# load PSet
TPLL_MultiAssayExp <- readRDS('./TPLL_MultiAssay.RDS')

# list of assay data in the object
experiments(TPLL_MultiAssayExp)

# list of sample and patient Ids (p number and TP number)
cols <- data.frame(colData(TPLL_MultiAssayExp))

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


