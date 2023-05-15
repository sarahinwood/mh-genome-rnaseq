#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

#############
# LIBRARIES #
#############

library(tximport)
library(data.table)
library(DESeq2)
library(pheatmap)

###########
# GLOBALS #
###########

MhV1_dds <- snakemake@input[["MhV1_dds"]]

########
# MAIN #
########

mh_dds <- readRDS(MhV1_dds)

##select only adult tissue samples
mh_dds_adult <- mh_dds[,mh_dds$stage=="Adult"]
##re-level tissue factor
mh_dds_adult$Tissue <- factor(mh_dds_adult$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mh_dds_adult$Rep <- factor(mh_dds_adult$rep)

design(mh_dds_adult) <- ~Rep+Tissue
mh_dds_adult <- DESeq(mh_dds_adult, test="LRT", reduced=~Rep)
res_group <- results(mh_dds_adult, alpha = 0.05)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, snakemake@output[["degs"]])

# write log
sessionInfo()