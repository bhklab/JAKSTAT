counts_annotation <- read.csv("./data/counts_annotation.csv", row.names=1)

cd247 <- counts_annotation[counts_annotation$gene_name == 'CD247', ]