library(data.table)
library(DESeq2)

mhv1_dds <- readRDS("output/05_deseq2/MhV1/MhV1_dds.rds")

##re-level factors
mhv1_dds$Stage <- factor(mhv1_dds$stage)
mhv1_dds$Flowcell <- factor(mhv1_dds$Flowcell)
design(mhv1_dds) <- ~Stage
mhv1_dds <- DESeq(mhv1_dds)

res_group <- results(mhv1_dds, contrast=c("Stage", "Pupa", "Adult"), alpha = 0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
