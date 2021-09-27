library(data.table)
library(VennDiagram)

##LRT results
ordered_sig_res_group_table <- fread("output/deseq2/tissue_LRT/sig_degs.csv")
##reduce to 0.01 and see if overlap is maintained
ordered_sig_res_group_table_01 <- subset(ordered_sig_res_group_table, padj<0.01)

##compare genes found in tissue itWT with LRT DEGs

###########
## ovary ##
###########
ovary_specific_DEGs <- fread("output/deseq2/itWT/ovary/mh_ovary_specific_degs.txt", header=FALSE)
vd_ov <- venn.diagram(x = list("itWT Ovary DEGs"=ovary_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ov)
LRT_ovary_degs <- subset(ordered_sig_res_group_table, rn %in% ovary_specific_DEGs$V1)
fwrite(LRT_ovary_degs, "output/deseq2/tissue_LRT/ovary/ovary_sp_LRT_annots.csv")
fwrite((list(LRT_ovary_degs$rn)), "output/deseq2/tissue_LRT/ovary/ovary_degs_ids.txt")

###########
## venom ##
###########
venom_specific_DEGs <- fread("output/deseq2/itWT/venom/mh_venom_specific_degs.txt", header=FALSE)
vd_ven <- venn.diagram(x = list("itWT Venom DEGs"=venom_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_ven)
LRT_venom_degs <- subset(ordered_sig_res_group_table, rn %in% venom_specific_DEGs$V1)
fwrite(LRT_venom_degs, "output/deseq2/tissue_LRT/venom/venom_sp_LRT_annots.csv")
fwrite((list(LRT_venom_degs$rn)), "output/deseq2/tissue_LRT/venom/venom_degs_ids.txt")

##########
## head ##
##########
head_specific_DEGs <- fread("output/deseq2/itWT/head/mh_head_specific_degs.txt", header=FALSE)
vd_head <- venn.diagram(x = list("itWT head DEGs"=head_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_head)
LRT_head_degs <- subset(ordered_sig_res_group_table, rn %in% head_specific_DEGs$V1)
fwrite(LRT_head_degs, "output/deseq2/tissue_LRT/head/head_sp_LRT_annots.csv")
fwrite((list(LRT_head_degs$rn)), "output/deseq2/tissue_LRT/head/head_degs_ids.txt")

############
## thorax ##
############
thorax_specific_DEGs <- fread("output/deseq2/itWT/thorax/mh_thorax_specific_degs.txt", header=FALSE)
vd_thor <- venn.diagram(x = list("itWT thorax DEGs"=thorax_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_thor)
LRT_thorax_degs <- subset(ordered_sig_res_group_table, rn %in% thorax_specific_DEGs$V1)
fwrite(LRT_thorax_degs, "output/deseq2/tissue_LRT/thorax/thorax_sp_LRT_annots.csv")
fwrite((list(LRT_thorax_degs$rn)), "output/deseq2/tissue_LRT/thorax/thorax_degs_ids.txt")

##########
## abdo ##
##########
abdo_specific_DEGs <- fread("output/deseq2/itWT/abdomen/mh_abdo_specific_degs.txt", header=FALSE)
vd_abdo <- venn.diagram(x = list("itWT abdo DEGs"=abdo_specific_DEGs$V1, "LRT"=ordered_sig_res_group_table_01$rn), filename=NULL, alpha=0.7, cex = 1, cat.cex=1, lwd=1.5)
grid.newpage()
grid.draw(vd_abdo)
LRT_abdo_degs <- subset(ordered_sig_res_group_table, rn %in% abdo_specific_DEGs$V1)
fwrite(LRT_abdo_degs, "output/deseq2/tissue_LRT/abdomen/abdo_sp_LRT_annots.csv")
fwrite((list(LRT_abdo_degs$rn)), "output/deseq2/tissue_LRT/abdomen/abdomen_degs_ids.txt")

