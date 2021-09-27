library(data.table)
library(DESeq2)

#finding the count files from star
count_files<- list.files(path="output/star/star_pass2", pattern="ReadsPerGene.out.tab", full.names = TRUE)

#importing count files
starCounts<- do.call(cbind,
                     lapply(count_files, read.table, header = FALSE, sep = "\t", row.names = 1,
                            colClasses = c("character", rep("NULL", 2), "integer")))
##4th col = reverse stranded counts, don't need middle 2

#naming columns with full file names
colnames(starCounts)<- gsub(".*/.*/.*/(.+).*", "", basename(count_files))
#renaming it wothout the extra stuff at the end
colnames(starCounts)<-gsub("(.+).ReadsPerGene.out.tab", "\\1", colnames(starCounts))

#discard top 4vrows that aren't genes
remove_row_names <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous", "MissingGeneID")
starCounts_final <- starCounts[!(row.names(starCounts) %in% remove_row_names), ]

#read in files with sample names and all of conditions
sample_data<- fread("data/full_sample_table.csv", header=TRUE)
setkey(sample_data, sample_name)

#Create dds object and link to sample data info- get table with counts- like table with sample names and treatments in a way Deseq undertands
dds <-DESeqDataSetFromMatrix(countData=starCounts_final, colData = sample_data, design = ~1)
colData(dds)

#save dds as file
saveRDS(dds, "output/deseq2/mh_dds_all.rds")

##remove venom3 sample
mh_dds_ven <- dds[,-c(18)]
saveRDS(mh_dds_ven, "output/deseq2/mh_dds.rds")
