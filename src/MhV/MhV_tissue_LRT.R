library(data.table)
library(DESeq2)
library(viridis)
library(pheatmap)
library(tidyverse)

mhv1_dds <- readRDS("output/05_deseq2/MhV1/MhV1_dds.rds")

##only adult
mhv1_dds_tissues <- mhv1_dds[,mhv1_dds$stage=="Adult"]

##re-level factors
mhv1_dds_tissues$Tissue <- factor(mhv1_dds_tissues$tissue, levels=c("Head", "Thorax", "Abdomen", "Ovaries", "Venom"))
mhv1_dds_tissues$Batch <- factor(mhv1_dds_tissues$rep)
mhv1_dds_tissues$Flowcell <- factor(mhv1_dds_tissues$Flowcell)
design(mhv1_dds_tissues) <- ~Flowcell+Tissue
mhv1_dds_tissues <- DESeq(mhv1_dds_tissues, test="LRT", reduced=~Flowcell)

res_group <- results(mhv1_dds_tissues, alpha = 0.05)
summary(res_group)
ordered_res_group <- res_group[order(res_group$padj),]
##Make data table and write to output
ordered_res_group_table <- data.table(data.frame(ordered_res_group), keep.rownames = TRUE)
ordered_sig_res_group_table <- subset(ordered_res_group_table, padj < 0.05)

#############
## heatmap ##
#############

## vst transform
mh_vst <- varianceStabilizingTransformation(mhv1_dds_tissues, blind=TRUE)
mh_vst_assay_dt <- data.table(assay(mh_vst), keep.rownames=TRUE)

mh_vst_all <- mh_vst_assay_dt %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_all_plot <- mh_vst_all[,c(4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18,13,14)]

## subset for DEGs
mh_vst_degs <- subset(mh_vst_assay_dt, rn %in% ordered_sig_res_group_table$rn)
## turn first row back to row name
mh_vst_degs <- mh_vst_degs %>% remove_rownames %>% column_to_rownames(var="rn")
## reorder tissues for plot - tissue of interest first
mh_vst_degs_plot <- mh_vst_degs[,c(4,5,6,16,10,11,12,1,2,3,15,7,8,9,17,18,13,14)]

## get tissue label info
sample_to_tissue <- data.table(data.frame(colData(mhv1_dds_tissues)[,c("Tissue", "sample_name")]))
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")
## for plot label
sample_to_tissue <- as.data.frame(colData(mhv1_dds_tissues)[,c("Tissue", "sample_name")])
sample_to_tissue <- sample_to_tissue %>% remove_rownames %>% column_to_rownames(var="sample_name")

# sample colour labels
tissue_colours <- list(Tissue = c(Head='#2D1160FF', Thorax='#721F81FF', Abdomen='#B63679FF', Ovaries='#F1605DFF', Venom='#FEAF77FF'))

## plot
## not clustered by samples
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=FALSE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

pheatmap(mh_vst_all_plot, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))

## clustered by samples
pheatmap(mh_vst_degs_plot, cluster_rows=TRUE, cluster_cols=F, show_rownames=T,
         annotation_col=sample_to_tissue, annotation_colors=tissue_colours, annotation_names_col=FALSE,
         show_colnames = FALSE, border_color=NA, color=viridis(50))
