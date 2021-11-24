library(stringr)

files <- read.delim("filenames.txt", header=FALSE)
sample.names <- str_extract(files$V1, "\\d+_S\\d+")
sample.names <- unique(sample.names)

read1 <- lapply(sample.names, function(sample){
  return(files$V1[str_which(files$V1, paste0(sample, "_L\\d+_R1"))])
})
read1 <- unlist(read1)

read2 <- lapply(sample.names, function(sample){
  return(files$V1[str_which(files$V1, paste0(sample, "_L\\d+_R2"))])
})
read2 <- unlist(read2)

tsv <- data.frame(matrix(data=NA, ncol=3, nrow=length(sample.names)))
colnames(tsv) <- c("sample", "read1", "read2")
tsv$sample <- sample.names
tsv$read1 <- read1
tsv$read2 <- read2

write.table(tsv, file='samples.tsv', quote=FALSE, sep='\t', row.names = FALSE)