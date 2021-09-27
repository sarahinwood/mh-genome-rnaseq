library(data.table)
library(dplyr)
library(VennDiagram)
library(stringr)

#Get annotation report in right format for venn diagram
annotation.report <- fread('output/Mh_venom2_transcriptome/trinotate_Mh_venom2/trinotate/trinotate_annotation_report.txt', na.strings = ".")
pfam <- annotation.report[!is.na(Pfam), unique(`#gene_id`)]
blastx <- annotation.report[!is.na(sprot_Top_BLASTX_hit), unique(`#gene_id`)]
blastp <- annotation.report[!is.na(sprot_Top_BLASTP_hit), unique(`#gene_id`)]
kegg <- annotation.report[!is.na(Kegg), unique(`#gene_id`)]
pfam_go <- annotation.report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
bx_go <- annotation.report[!is.na(gene_ontology_BLASTX), unique(`#gene_id`)]
bp_go <- annotation.report[!is.na(gene_ontology_BLASTP), unique(`#gene_id`)]

number.genes <- annotation.report[!is.na(`#gene_id`),length(unique(`#gene_id`))]

#Draw Venn Diagram
vd <- venn.diagram(x = list("Pfam"=pfam, "BlastX"=blastx, "Kegg"=kegg), filename=NULL,
                   fill=c("#440154FF", "#21908CFF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5,
                   main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd)

vd2 <- venn.diagram(x = list("Pfam GO"=pfam_go, "BlastX GO"=bx_go, "BlastP GO"=bp_go), filename=NULL,
                    fill=c("#440154FF", "#21908CFF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5,
                    main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd2)

#Sum of genes with any annotation
long.annotationreport <- melt(annotation.report,id.vars = "#gene_id", measure.vars = c("sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "Pfam",  "eggnog", "Kegg"))
any.annotations <- long.annotationreport[,.(any_annotations = any(!is.na(value))),by=`#gene_id`]
##percentage with any annotation
(any.annotations[,sum(any_annotations)]/any.annotations[,length(unique(`#gene_id`))])*100

#Sum of genes with blast annotation
long.blastreport <- melt(annotation.report,id.vars = "#gene_id", measure.vars = c("sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit"))
any.blast <- long.blastreport[,.(any_annotations = any(!is.na(value))),by=`#gene_id`]
any.blast[,length(unique(`#gene_id`))]
any.blast[,sum(any_annotations)]
any_blast_ids <- subset(any.blast, any.blast$any_annotations==TRUE)

#Sum of genes with any Transdecoder predicted protein
transdecoder.report <- melt(annotation.report, id.vars="#gene_id", measure.vars = c("prot_id"))
any.transdecoder <- transdecoder.report[,.(transdecoder_annotation = any(!is.na(value))),by=`#gene_id`]
any.transdecoder[,length(unique(`#gene_id`))]
any.transdecoder[,sum(transdecoder_annotation)]
sum(any.transdecoder$transdecoder_annotation==TRUE)

any_annot_ids <- subset(any.annotations, any.annotations$any_annotations==TRUE)
any_transdecoder_ids <- subset(any.transdecoder, any.transdecoder$transdecoder_annotation==TRUE)

vd3 <- venn.diagram(x = list("Transdecoder"=any_transdecoder_ids$`#gene_id`, "Blast"=any_blast_ids$`#gene_id`), filename=NULL,
                    fill=c("#440154FF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5, height=3, width=6,
                    main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd3)

##split blastX results
blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
split.first.blastx$e_value <- tstrsplit(split.first.blastx$V5, "E:", keep=c(2))
split.first.blastx$pid <- tstrsplit(split.first.blastx$V4, "%ID", keep=c(1))
gene_pid_eval <- split.first.blastx[!is.na(e_value),.(`#gene_id`, pid, e_value)]
gene_pid_eval$pid <- as.numeric(gene_pid_eval$pid)
gene_pid_eval$e_value <- as.numeric(gene_pid_eval$e_value)

mean(gene_pid_eval$pid)

#Annotations per genus
genes.per.genus <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
setkey(genes.per.genus, V1)
print(genes.per.genus)
fwrite(genes.per.genus, "output/trinotate_Mh_venom2/genes_per_genus.csv")
