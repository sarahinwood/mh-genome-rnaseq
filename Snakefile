#!/usr/bin/env python3
import peppy
import pathlib2

###########
# GLOBALS #
###########

##needed to get BUSCO running in new folder
def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

#this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
#can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

#containers
star_container = 'https://github.com/deardenlab/container-star/releases/download/0.0.1/container-star.2.7.9a.sif'
bbduk_container = 'https://github.com/deardenlab/container-bbmap/releases/download/0.0.3/container-bbmap.bbmap_38.90.sif'
tidyverse_container = 'docker://rocker/tidyverse:4.1.0'
multiqc_container = 'docker://ewels/multiqc:1.9'
kraken_container = 'docker://staphb/kraken2:2.1.2-no-db'
trinity_container = 'docker://trinityrnaseq/trinityrnaseq:2.12.0'
trinotate_container = 'shub://TomHarrop/trinotate_pipeline:v0.0.12'
busco_container = 'docker://ezlabgva/busco:v5.2.1_cv1'

###########
## RULES ##
###########

rule target:
    input:
        'output/multiqc_filtermatch/multiqc_report.html',
        'output/star/star_pass2/Mh_venom3.ReadsPerGene.out.tab',
        expand('output/kraken_unmapped/reports/kraken_{sample}_report.txt', sample=all_samples),
        expand('output/star_filtermatch/star_pass2/{sample}.ReadsPerGene.out.tab', sample=all_samples),
        'output/Mh_annots.gtf',
        'output/star/star_reference/Genome',
        ##venom2 unmapped
        'output/Mh_venom2_transcriptome/trinotate_Mh_venom2/trinotate/trinotate_annotation_report.txt',
        'output/kraken_unmapped/reports/kraken_Mh_venom3_report.txt'

###################################
## what are unmapped reads from? ##
###################################

rule trinotate_v2_unmapped:
    input:
        fasta = 'output/Mh_venom2_transcriptome/trinity_Mh_venom2/Trinity.fasta',
        blastdb = 'bin/trinotate_db/uniprot_sprot.pep',
        hmmerdb = 'bin/trinotate_db/Pfam-A.hmm',
        sqldb = 'bin/trinotate_db/Trinotate.sqlite'
    output:
        'output/Mh_venom2_transcriptome/trinotate_Mh_venom2/trinotate/trinotate_annotation_report.txt',
        'output/Mh_venom2_transcriptome/trinotate_Mh_venom2/trinotate/Trinotate.sqlite'
    params:
        wd = 'output/Mh_venom2_transcriptome/trinotate_Mh_venom2'
    threads:
        50
    log:
        'output/logs/trinotate.log'
    singularity:
        trinotate_container
    shell:
        'trinotate_pipeline '
        '--trinity_fasta {input.fasta} '
        '--blast_db {input.blastdb} '
        '--hmmer_db {input.hmmerdb} '
        '--sqlite_db {input.sqldb} '
        '--outdir {params.wd} '
        '--threads {threads} '
        '&> {log}'

rule Trinity_v2_unmapped:
    input:
        r2 = 'output/star/star_pass2/Mh_venom2.Unmapped.out.mate1',
        r1 = 'output/star/star_pass2/Mh_venom2.Unmapped.out.mate2',
    output:
        'output/Mh_venom2_transcriptome/trinity/Trinity.fasta',
        'output/Mh_venom2_transcriptome/trinity/Trinity.fasta.gene_trans_map'
    params:
        outdir = 'output/Mh_venom2_transcriptome/trinity'
    singularity:
        trinity_container
    threads:
        50
    log:
        'output/logs/trinity.log'
    shell:
        'Trinity '
        '--SS_lib_type RF '
        '--max_memory 300G '
        '--CPU {threads} '
        '--output {params.outdir} '
        '--left {input.r1} '
        '--right {input.r2} '
        '--seqType fq '
        '&> {log}'

rule kraken_unmapped_v3:
    input:
        r1 = 'output/star/star_pass2/Mh_venom3.Unmapped.out.mate1',
        r2 = 'output/star/star_pass2/Mh_venom3.Unmapped.out.mate2',
        db = 'bin/db/kraken_std'
    output:
        out = 'output/kraken_unmapped/out/kraken_Mh_venom3_out.txt',
        report = 'output/kraken_unmapped/reports/kraken_Mh_venom3_report.txt'
    log:
        'output/logs/kraken_unmapped/kraken_Mh_venom3.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

##run kraken on unmapped reads from STAR - are unmapped reads something prokaryotic or unclassified?
rule kraken_unmapped:
    input:
        r1 = 'output/star/star_pass2/{sample}.Unmapped.out.mate1',
        r2 = 'output/star/star_pass2/{sample}.Unmapped.out.mate2',
        db = 'bin/db/kraken_std'
    output:
        out = 'output/kraken_unmapped/out/kraken_{sample}_out.txt',
        report = 'output/kraken_unmapped/reports/kraken_{sample}_report.txt'
    log:
        'output/logs/kraken_unmapped/kraken_{sample}.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

#################
## mapping etc ##
#################

rule multiqc:
    input:
        expand('output/star_filtermatch/star_pass2/{sample}.ReadsPerGene.out.tab', sample=all_samples)
    output:
        'output/multiqc_filtermatch/multiqc_report.html'
    params:
        outdir = 'output/multiqc_filtermatch',
        indirs = ['output/star_filtermatch/star_pass2', 'output/fastqc']
    log:
        'output/logs/multiqc.log'
    container:
        multiqc_container
    shell:
        'multiqc '
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'

##index so can view in igv
rule index_star_bam:
    input:
        bam = 'output/star_filtermatch/star_pass2/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        bai = 'output/star_filtermatch/star_pass2/{sample}.Aligned.sortedByCoord.out.bam.bai'
    threads:
        20
    log:
        'output/logs/{sample}_index.log'
    shell:
        'samtools '
        'index '
        '{input.bam} '

rule star_second_pass_v3:
    input:
        left = 'output/bbduk_trim/Mh_venom3_r1.fq.gz',
        right = 'output/bbduk_trim/Mh_venom3_r2.fq.gz',
        v3_junctions = 'output/star/star_pass1/Mh_venom3.SJ.out.tab',
        junctions = expand('output/star/star_pass1/{sample}.SJ.out.tab', sample=all_samples)
    output:
        bam = 'output/star/star_pass2/Mh_venom3.Aligned.sortedByCoord.out.bam',
        reads_per_gene = 'output/star/star_pass2/Mh_venom3.ReadsPerGene.out.tab',
        unmapped = 'output/star/star_pass2/Mh_venom3.Unmapped.out.mate1'
    threads:
        20
    params:
        genome_dir = 'output/star/star_reference',
        prefix = 'output/star/star_pass2/Mh_venom3.'
    log:
        'output/logs/star/star_pass2_Mh_venom3.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outBAMcompression 10 '
        '--quantMode GeneCounts '
        '--readFilesIn {input.left} {input.right} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '--outReadsUnmapped Fastx '
        '&> {log}'

rule star_first_pass_v3:
    input:
        left = 'output/bbduk_trim/Mh_venom3_r1.fq.gz',
        right = 'output/bbduk_trim/Mh_venom3_r2.fq.gz',
        star_reference = 'output/star/star_reference/Genome'
    output:
        sjdb = 'output/star/star_pass1/Mh_venom3.SJ.out.tab'
    params:
        genome_dir = 'output/star/star_reference',
        prefix = 'output/star/star_pass1/Mh_venom3.'
    threads:
        20
    log:
        'output/logs/star/star_pass1_Mh_venom3.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '
        '--readFilesIn {input.left} {input.right} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_second_pass:
    input:
        left = 'output/bbduk_trim/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim/{sample}_r2.fq.gz',
        junctions = expand('output/star_filtermatch/star_pass1/{sample}.SJ.out.tab', sample=all_samples)
    output:
        bam = 'output/star_filtermatch/star_pass2/{sample}.Aligned.sortedByCoord.out.bam',
        reads_per_gene = 'output/star_filtermatch/star_pass2/{sample}.ReadsPerGene.out.tab',
        unmapped = 'output/star_filtermatch/star_pass2/{sample}.Unmapped.out.mate1'
    threads:
        20
    params:
        genome_dir = 'output/star/star_reference',
        prefix = 'output/star_filtermatch/star_pass2/{sample}.'
    log:
        'output/logs/star_filtermatch/star_pass2_{sample}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outFilterScoreMinOverLread 0.5 '
        '--outFilterMatchNminOverLread 0.5 '
        '--outSAMtype BAM SortedByCoordinate '
        '--outBAMcompression 10 '
        '--quantMode GeneCounts '
        '--readFilesIn {input.left} {input.right} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '--outReadsUnmapped Fastx '
        '&> {log}'

rule star_first_pass:
    input:
        left = 'output/bbduk_trim/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim/{sample}_r2.fq.gz',
        star_reference = 'output/star/star_reference/Genome'
    output:
        sjdb = 'output/star_filtermatch/star_pass1/{sample}.SJ.out.tab'
    params:
        genome_dir = 'output/star/star_reference',
        prefix = 'output/star_filtermatch/star_pass1/{sample}.'
    threads:
        20
    log:
        'output/logs/star/star_pass1_{sample}.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outFilterScoreMinOverLread 0.5 '
        '--outFilterMatchNminOverLread 0.5 '
        '--outSAMtype None '
        '--readFilesIn {input.left} {input.right} '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'

rule star_reference:
    input:
        mh_genome = 'data/final_microctonus_assemblies_annotations/Mh.scaffolds.fa',
        gtf = 'output/Mh_annots.gtf'
    output:
        'output/star/star_reference/Genome'
    params:
        genome_dir = 'output/star/star_reference'
    threads:
        20
    log:
        'output/logs/star/star_reference.log'
    singularity:
    	star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeSAindexNbases 12 '
        '--genomeDir {params.genome_dir} '
        '--genomeFastaFiles {input.mh_genome} '
        '--sjdbGTFfile {input.gtf} '
        '2> {log} '

rule convert_gff_to_gtf:
    input:
        'output/mh_genome_viral.gff3'
    output:
        'output/Mh_annots.gtf'
    log:
        'output/logs/convert_gff_to_gtf.log'
    shell:
        'gffread '
        '{input} '
        '-T '
        '-o {output} '
        '&> {log}'

