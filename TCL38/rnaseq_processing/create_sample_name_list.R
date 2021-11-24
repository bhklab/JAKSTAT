library(stringr)

extract_sample <- function(sample_list) {
  samples <- lapply(sample_list, function(row){
    id <- str_extract(row, "^\\d+(?=\\s)")
    sample.name <- str_extract(row, "(?<=\\s)[^\\s]+")
    return(c(id, sample.name))
  })
  return(samples)
}

rnaseq.samples <- read.table("./data/samples.tsv", header=TRUE, sep="\t")
sample.names.1 <- read.table("./data/sample_names_1.tab", header=TRUE, sep="\t")
sample.names.2 <- read.table("./data/sample_names_2.tab", header=TRUE, sep="\t")
sample.names.rerun <- read.table("./data/sample_names_rerun.tab", header=TRUE, sep="\t")
rnaseq.rerun.sample.names <- c("157085_S53", "157098_S54", "157104_S55", "157110_S56", "157124_S57", "157128_S58", "157132_S59")
tcl38 <- data.frame(read.delim('../misc/data/cells.txt'))
tcl38$Cell.line <- gsub("[*]|[[:space:]+]", "", tcl38$Cell.line)
  
samples.list <- extract_sample(sample.names.1$CCG.Sample.ID....Sample.Name)
samples.list <- c(samples.list, extract_sample(sample.names.2$CCG.Sample.ID.....Sample.Name))

id <- lapply(samples.list, function(vec){
  return(vec[1])
})
id <- unlist(id)
sample.names <- lapply(samples.list, function(vec){
  return(vec[2])
})
sample.names <- unlist(sample.names)

col.names <- c("CCG Sample ID", "Sample Name", "Cell", "RNA-seq Sample", "Rerun Sample")
samples.df <- data.frame(matrix(data=NA, ncol=length(col.names), nrow=length(id)))
colnames(samples.df) <- col.names
samples.df$`CCG Sample ID` <- id
samples.df$`Sample Name` <- sample.names
samples.df <- samples.df[order(samples.df$`CCG Sample ID`), ]
rownames(samples.df) <- samples.df$`CCG Sample ID`

rnaseq.samples.rerun <- rnaseq.samples[rnaseq.samples$sample %in% rnaseq.rerun.sample.names, ]
rnaseq.samples <- rnaseq.samples[!(rnaseq.samples$sample %in% rnaseq.rerun.sample.names), ]

rnaseq.sample.names <- lapply(rownames(samples.df), function(id){
  filtered <- rnaseq.samples[str_detect(rnaseq.samples$sample, id), ]
  return(filtered$sample)
})
samples.df$`RNA-seq Sample` <- unlist(rnaseq.sample.names)

rnaseq.sample.names <- lapply(rownames(samples.df), function(id){
  filtered <- rnaseq.samples.rerun[str_detect(rnaseq.samples.rerun$sample, id), ]
  if(dim(filtered)[1] == 0){
    return(NA)
  }
  return(filtered$sample)
})
samples.df$`Rerun Sample` <- unlist(rnaseq.sample.names)

for(cellname in tcl38$Cell.line){
  samples.df[str_detect(tolower(str_replace_all(samples.df$`Sample Name`, "-", "")), tolower(str_replace_all(cellname, "-", ""))), ]$Cell <- cellname
}

samples.df[str_detect(samples.df$`Sample Name`, "(?<=\\d)BC(?=\\d+)"), ]$Cell <- "BC"

write.csv(samples.df, "./data/rnaseq_sample_summary.csv")