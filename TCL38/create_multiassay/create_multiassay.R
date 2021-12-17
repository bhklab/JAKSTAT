library(stringr)
library(dplyr)
library(MultiAssayExperiment)

# create sample dataframe
acgh_samples <- read.csv("./data/acgh_samples.csv", row.names=1)
rnaseq_samples <- read.csv("./data/rnaseq_samples.csv", row.names=1)

acgh_cells <- unique(acgh_samples$cell)
rnaseq_cells <- unique(rnaseq_samples$cell)

# create sample info df
samples_df <- as.data.frame(matrix(nrow=39, ncol=2, dimnames=list(NULL, c("acgh_cells", "rnaseq_cells"))))
samples_df$rnaseq_cells <- rnaseq_cells

# map cell names in aCGH and RNA-seq data, and come up with a common name as the row name.
cellname_map <- list(`CCRF-H-SB2`="hsb2", `OCI-Ly13.2`="oci13", `Karpas299`="k299", `BC`=NA)
for(cellname in rnaseq_cells){
  found <- acgh_cells[str_detect(tolower(str_replace_all(acgh_cells, "[-_\\.]", "")), tolower(str_replace_all(cellname, "[-_\\.]", "")))]
  samples_df[samples_df$rnaseq_cells == cellname, ]$acgh_cells <- ifelse(length(found) > 0, found, cellname_map[[cellname]])
}
rownames(samples_df) <- str_replace_all(samples_df$rnaseq_cells, "[-\\.]", "_")

# add RNA-seq samples data
samples_df$rnaseq_sample_names <- lapply(samples_df$rnaseq_cells, function(cell){
  return(paste(rnaseq_samples[rnaseq_samples$cell == cell, ]$sample_name, collapse=", "))
})
samples_df$rnaseq_kallisto_sample_names <- lapply(samples_df$rnaseq_cells, function(cell){
  return(paste(rnaseq_samples[rnaseq_samples$cell == cell, ]$kallisto_sample_name, collapse=", "))
})
samples_df$rnaseq_healthy <- lapply(samples_df$rnaseq_cells, function(cell){
  ifelse(cell == "BC", "Yes", "")
})

# merge aCGH samples data
samples_df[colnames(acgh_samples)[colnames(acgh_samples) != "cell"]] <- NA
for(cell in acgh_cells){
  values <- as.list(acgh_samples[acgh_samples$cell == cell, -which(colnames(acgh_samples) == "cell")])
  for(name in names(values)){
    samples_df[which(samples_df$acgh_cells == cell), which(colnames(samples_df) == name)] <- values[[name]]
  }
}
colnames(samples_df)[colnames(samples_df) %in% colnames(acgh_samples)] <- paste0("acgh_", colnames(samples_df)[colnames(samples_df) %in% colnames(acgh_samples)])

# read RNA-seq data
rnaseq_counts <- read.csv("./data/counts_averaged.csv", row.name=1)