library(rtracklayer)
library(data.table)
library(tidyverse)

####################################
## reformat Prodigal GFF for STAR ##
####################################
mh_gff <- data.frame(readGFF("data/final_microctonus_assemblies_annotations/Mh.gff3"))
viral_contig_ids <- fread("data/Mh_viral_contig_ids.txt", header=FALSE)
prodigal_gff <- data.frame(readGFF("data/Mh_prodigal/gene_predictions.gff"))

##prodigal doesn't have parent column - add in, drop unnecessary col.s
red_prod_gff <- prodigal_gff[,c(1:10)]
red_prod_gff$ID <- paste("MHPROD", red_prod_gff$ID, sep="_")
red_prod_gff$Parent <- red_prod_gff$ID

##CDS rows
red_prod_gff_CDS <- red_prod_gff
red_prod_gff_CDS$ID <- paste(red_prod_gff_CDS$ID, "-T1.cds", sep="")
red_prod_gff_CDS$Parent <- paste(red_prod_gff_CDS$Parent, "-T1", sep="")

##exon rows
red_prod_gff_exon <- red_prod_gff
red_prod_gff_exon$type <- paste("exon")
red_prod_gff_exon$ID <- paste(red_prod_gff_exon$ID, "-T1.exon1", sep="")
red_prod_gff_exon$Parent <- paste(red_prod_gff_exon$Parent, "-T1", sep="")

##gene rows
red_prod_gff_gene <- red_prod_gff
red_prod_gff_gene$type <- paste("gene")
red_prod_gff_gene$Parent <- paste("character(0)")

##mRNA rows
red_prod_gff_mrna <- red_prod_gff
red_prod_gff_mrna$type <- paste("mRNA")
red_prod_gff_mrna$ID <- paste(red_prod_gff_mrna$ID, "-T1", sep="")
red_prod_gff_mrna$product <- paste("hypothetical protein")

##join all prodigal gff tables
full_prodigal_gff <- list(red_prod_gff_CDS, red_prod_gff_exon, red_prod_gff_gene, red_prod_gff_mrna) %>% reduce(full_join)
##write full prodigal GFF
export(full_prodigal_gff, "output/mh_prodigal_full.gff3", format="GFF3")

##################################
## join with rest of genome GFF ##
##################################

##remove euk-coded viral predictions
mh_gff_nonviral <- subset(mh_gff, !(seqid %in% viral_contig_ids$V1))
mh_gff_nonviral <- mh_gff_nonviral %>% mutate(Parent = as.character(Parent))

##add prodigal predictions instead
mh_viral_gff <- full_join(full_prodigal_gff, mh_gff_nonviral)
##write GFF3 file
export(mh_viral_gff, "output/mh_genome_viral.gff3", format="GFF3")
