library(tximport)
library(data.table)
library(DESeq2)
library(VennDiagram)

mh_dds <- readRDS("output/deseq2/mh_dds.rds")

##select only adult tissue samples
mh_dds_adult <- mh_dds[,mh_dds$stage=="Adult"]
##re-level tissue factor
mh_dds_adult$Tissue <- factor(mh_dds_adult$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_adult$Rep <- factor(mh_dds_adult$rep)

design(mh_dds_adult) <- ~Rep+Tissue
mh_dds_adult <- DESeq(mh_dds_adult, test="LRT", reduced=~Rep)
saveRDS(mh_dds_adult, "output/deseq2/tissue_LRT/mh_tissue_LRT.rds")

res_group <- results(mh_dds_adult, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, "output/deseq2/tissue_LRT/sig_degs.csv")

plotCounts(mh_dds_adult, "TRINITY_DN5452_c0_g1", intgroup=("tissue"))
