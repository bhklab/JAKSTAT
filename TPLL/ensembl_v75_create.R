library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)

options(stringsAsFactors = FALSE)

txdb <- makeTxDbFromGFF("~/Desktop/Homo_sapiens.GRCh37.75.gtf")
saveDb(txdb, "~/Desktop/ensembl.v75.annotation.txdb")

txdb <- loadDb("~/Desktop/ensembl.v75.annotation.txdb")
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene) <- c("transcripts", "genes")

features <- rtracklayer::import('~/Desktop/Homo_sapiens.GRCh37.75.gtf')
features_gene = as.data.frame(features)
features_gene = features_gene[which(features_gene$type == "gene"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","tag","transcript_name","transcript_version","transcript_source","havana_transcript","transcript_biotype" ,"ccds_id","protein_version","ont", "protein_id","exon_id", "exon_number", "exon_version")
features_gene <- features_gene[,!names(features_gene) %in% remove]
rownames(features_gene) <- features_gene$gene_id


features_transcript = as.data.frame(features)
features_transcript = features_transcript[which(features_transcript$type == "transcript"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("exon_id", "exon_number", "exon_version")
features_transcript <- features_transcript[,!names(features_transcript) %in% remove]
rownames(features_transcript) <- features_transcript$transcript_id
save(tx2gene, features_gene, features_transcript, file="~/Desktop/Ensembl.v75.annotation.RData")