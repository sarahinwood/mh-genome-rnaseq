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
library(tidyverse)
library(data.table)
library(DESeq2)
library(viridis)
library(pheatmap)

###########
# GLOBALS #
###########

gene2tx_file <- snakemake@input[["gene2tx_file"]]
sample_data_file <- snakemake@input[["sample_data_file"]]

########
# MAIN #
########

tx2gene <- fread(gene2tx_file, header = F)

##Find all salmon quant files
quant_files <- list.files(path="output/03_salmon/", pattern = "quant.sf", full.names=TRUE, recursive = TRUE)
##assign names to quant files from folder name
names(quant_files) <- gsub(".*/(.+)_quant/.*", "\\1", quant_files)
##import the salmon quant files (tx2gene links transcript ID to Gene ID
##required for gene-level summarisation for methods that only provide transcript level estimates e.g. salmon)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, dropInfReps=TRUE, abundanceCol="TPM")
txi_abundance <- data.table(txi$abundance, keep.rownames=TRUE)
fwrite(txi_abundance, snakemake@output[["salmon_tpm"]])

##Import table describing samples
sample_data <- fread(sample_data_file, header=TRUE)
setkey(sample_data, sample_name)

##create dds object and link to sample data
dds <- DESeqDataSetFromTximport(txi, colData = sample_data[colnames(txi$counts)], design = ~1)
dds <- estimateSizeFactors(dds)
##save dds object
saveRDS(dds, snakemake@output[["all_dds"]])

## split into Mh and MhV1

##mh gene list
mh_tx <- subset(tx2gene, grepl("MHPGA_", V1))
mh_gene <- unique(mh_tx$V1)
##subset mh dds
mh_dds <- dds[mh_gene,]
saveRDS(mh_dds, snakemake@output[["Mh_dds"]])

##MhV1 gene list
mhv1_tx <- subset(tx2gene, grepl("ORF", V1))
mhv1_gene <- unique(mhv1_tx$V1)
##subset MhV dds
mhv1_dds <- dds[mhv1_gene,]
saveRDS(mhv1_dds, snakemake@output[["MhV1_dds"]])

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mhv1_dds, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)
mh_vst_all <- mh_vst_assay_dt %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_all_plot <- mh_vst_all[,c(10,11,12,4,5,6,19,13,14,15,1,2,3,18,7,8,9,20,21,16,17)]

## get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mhv1_dds)[,c("tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

# sample colour labels
tissue_colours <- list(tissue = c(Pupa = '#fcfdbf', Head='#2D1160FF', Thorax='#721F81FF', Abdomen='#B63679FF', Ovaries='#F1605DFF', Venom='#FEAF77FF'))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

## plot
## not clustered by samples
pdf(snakemake@output[["MhV1_heatmap"]], height=12, width=7)
pheatmap(mh_vst_all_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=T, clustering_callback=callback,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50), fontsize=10, fontsize_row=6)
dev.off()

# write log
sessionInfo()