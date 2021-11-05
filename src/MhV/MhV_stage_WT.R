library(data.table)
library(DESeq2)
library(rtracklayer)

mh_gff <- data.frame(readGFF("data/Mh_prodigal/gene_predictions.gff"))
mh_gff_ids <- data.table(mh_gff[,c(1,9)])
mh_gff_ids$gene_no <- tstrsplit(mh_gff_ids$ID, "_", keep=c(2))
mh_gff_ids$gene_id <- paste(mh_gff_ids$seqid, mh_gff_ids$gene_no, sep="_")
gene_to_id <- mh_gff_ids[,c(2,4)]

mh_dds <- readRDS("output/deseq2/mh_dds.rds")
##estimate size factors on whole table - differences in library size between samples
mh_dds <- estimateSizeFactors(mh_dds)


##keep only viral genes
counts_table <- data.table(counts(mh_dds), keep.rownames=TRUE)
viral <- subset(counts_table, grepl("MHPROD", rn))
viral_genes <- viral$rn
mh_dds_viral <- mh_dds[viral_genes,]

##re-level factors
mh_dds_viral$Stage <- factor(mh_dds_viral$stage)
mh_dds_viral$Rep <- factor(mh_dds_viral$rep)
design(mh_dds_viral) <- ~Rep+Stage
mh_dds_viral <- DESeq(mh_dds_viral, test="LRT", reduced=~Rep)


res_group <- results(mh_dds_viral, alpha = 0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
ordered_sig_res_group_table$ID <- tstrsplit(ordered_sig_res_group_table$rn, "MHPROD_", keep=c(2))
sig_genes <- merge(ordered_sig_res_group_table, gene_to_id, by="ID")

