library(data.table)
library(DESeq2)
library(VennDiagram)

dds_tissue <- readRDS("output/deseq2/itWT/mh_itWT.rds")

###############################
## iterative pairwise comp.s ##
###############################

##head venom
head_venom <- results(dds_tissue, contrast=c("Tissue", "Head", "Venom"), alpha = 0.05, lfcThreshold = 1)
summary(head_venom)
##Order based of padj
ordered_head_venom <- head_venom[order(head_venom$padj),]
##Make data table and write to output
ordered_head_venom_table <- data.table(data.frame(ordered_head_venom), keep.rownames = TRUE)
venom <- subset(ordered_head_venom_table, padj < 0.05)

##head abdo
head_abdo <- results(dds_tissue, contrast=c("Tissue", "Head", "Abdomen"), alpha = 0.05, lfcThreshold = 1)
summary(head_abdo)
##Order based of padj
ordered_head_abdo <- head_abdo[order(head_abdo$padj),]
##Make data table and write to output
ordered_head_abdo_table <- data.table(data.frame(ordered_head_abdo), keep.rownames = TRUE)
abdo <- subset(ordered_head_abdo_table, padj < 0.05)

##head thorax
head_thorax <- results(dds_tissue, contrast=c("Tissue", "Head", "Thorax"), alpha = 0.05, lfcThreshold = 1)
summary(head_thorax)
##Order based of padj
ordered_head_thorax <- head_thorax[order(head_thorax$padj),]
##Make data table and write to output
ordered_head_thorax_table <- data.table(data.frame(ordered_head_thorax), keep.rownames = TRUE)
thorax <- subset(ordered_head_thorax_table, padj < 0.05)

##head ovaries 
head_ovaries <- results(dds_tissue, contrast=c("Tissue", "Head", "Ovaries"), alpha = 0.05, lfcThreshold = 1)
summary(head_ovaries)
##Order based of padj
ordered_head_ovaries <- head_ovaries[order(head_ovaries$padj),]
##Make data table and write to output
ordered_head_ovaries_table <- data.table(data.frame(ordered_head_ovaries), keep.rownames = TRUE)
ovaries <- subset(ordered_head_ovaries_table, padj < 0.05)

################################
## overlap for head-specific ##
################################

head_specific_DEGs <- intersect(intersect(intersect(venom$rn, thorax$rn), abdo$rn), ovaries$rn)


##venn diagram
vd1 <- venn.diagram(x = list("Venom"=venom$rn, "Thorax"=thorax$rn, "Abomen"=abdo$rn, "Ovaries"=ovaries$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd1)

fwrite(list(head_specific_DEGs), "output/deseq2/itWT/head/mh_head_specific_degs.txt")
