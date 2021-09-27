library(data.table)
library(stringr)
library(DESeq2)

dds_stage <- readRDS("output/deseq2/stage_WT/mh_stage_WT.rds")

stage_degs <- fread("output/deseq2/stage_WT/sig_degs.csv")

##filter prodigal genes
prodigal_DEGs <- stage_degs %>% filter(str_detect(stage_degs$rn, 'MHPROD'))

###############################
## plot multiple gene counts ##
###############################

##get gene counts
counts_table <- data.table(counts(dds_stage, normalized=TRUE), keep.rownames = TRUE)

#### filter by list of interest ####
annot_counts <- filter(counts_table, rn %in% prodigal_DEGs$rn)

##melt for plotting
plot_annots_counts <- annot_counts %>% gather(colnames(annot_counts)[2:18], key="sample_name", value="normalized_counts")
##sample group information
sample_table <- fread("data/sample_table.csv")
name_vs_group <- sample_table[,c(1,2)]
plotting_counts <- inner_join(plot_annots_counts, name_vs_group)
tissue_order <- c("Head", "Thorax", "Abdomen", "Ovaries", "Venom", "Pupa")
plotting_counts$tissue <- factor(plotting_counts$tissue, levels=tissue_order)
##add alphabetical label to each plot
plotting_counts$gene_label <- paste(plotting_counts$rn)
##plot all annot DEGs using ggplot2
ggplot(plotting_counts) +
  geom_point(aes(x = tissue, y = normalized_counts, colour=tissue)) +
  labs(colour="Tissue", y="Normalized counts", x="")+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  facet_wrap(~gene_label, scales="free")
