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
mh_dds$Stage <- factor(mh_dds$stage, levels=c("Pupa", "Adult"))
mh_dds$Rep <- factor(mh_dds$rep)

design(mh_dds) <- ~Rep+Stage
mh_dds <- DESeq(mh_dds)
res_group <- results(mh_dds, alpha = 0.05, lfcThreshold = 1)
summary(res_group)
##Order based of padj
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)
fwrite(ordered_sig_res_group_table, snakemake@output[["degs"]])

# write log
sessionInfo()