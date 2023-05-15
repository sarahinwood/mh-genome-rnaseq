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
salmon_container = 'docker://combinelab/salmon:1.5.1'
multiqc_container = 'docker://ewels/multiqc:1.9'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.14:0.0.1'

all_samples_minus_v2 = ["Mh_head1", "Mh_head2", "Mh_head3", "Mh_thorax1", "Mh_thorax2", "Mh_thorax3", "Mh_abdo1", "Mh_abdo2", "Mh_abdo3", 
"Mh_ovary1", "Mh_ovary2", "Mh_ovary3", "Mh_venom1", "Mh_pupa1", "Mh_pupa2", "Mh_pupa3", "Mhyp_old_head", "Mhyp_old_abdo", "Mhyp_old_ovary1", "Mhyp_old_ovary2"]

###########
## RULES ##
###########

rule target:
    input:
        'output/04_multiqc/multiqc_report.html',
        'output/05_deseq2/MhV1/stage_WT/sig_degs.csv',
        'output/05_deseq2/MhV1/tissue_LRT/sig_degs.csv'

############
## DESeq2 ##
############

rule stage_WT:
    input:
        MhV1_dds = 'output/05_deseq2/MhV1/MhV1_dds.rds'
    output:
        degs = 'output/05_deseq2/MhV1/stage_WT/sig_degs.csv'
    log:
        'output/logs/stage_WT_MhV1.log'
    singularity:
        bioconductor_container
    script:
        'src/stage_WT/stage_WT.R'

rule tissue_LRT:
    input:
        MhV1_dds = 'output/05_deseq2/MhV1/MhV1_dds.rds'
    output:
        degs = 'output/05_deseq2/MhV1/tissue_LRT/sig_degs.csv'
    log:
        'output/logs/tissue_LRT_MhV1.log'
    singularity:
        bioconductor_container
    script:
        'src/tissue_LRT/tissue_LRT.R'

rule make_dds:
    input:
        gene2tx_file = 'output/02_concat_Mh_MhV1_genomes/gene2tx.tsv',
        salmon_output = expand('output/03_salmon/{sample}_quant/quant.sf', sample=all_samples),
        sample_data_file = 'data/sample_table.csv'
    output:
        salmon_tpm = 'output/05_deseq2/salmon_tpm.csv',
        all_dds = 'output/05_deseq2/Mh_MhV1_dds.rds',
        Mh_dds = 'output/05_deseq2/Mh/Mh_dds.rds',
        MhV1_dds = 'output/05_deseq2/MhV1/MhV1_dds.rds',
        MhV1_heatmap = 'output/05_deseq2/MhV1/MhV1_heatmap.pdf'
    log:
        'output/logs/deseq2/make_dds.log'
    singularity:
        bioconductor_container
    script:
        'src/make_dds.R'

#################
## mapping etc ## mapping to hyp genome is terrible - can I even use this in paper because it might suggest bad annotations?
################# doesn't display multiple sets of salmon results in one file

rule multiqc:
    input:
        Mh_Mhv1 = expand('output/03_salmon/{sample}_quant/quant.sf', sample=all_samples),
        Mhv1_only = expand('output/03_salmon_MhV1_only/{sample}_quant/quant.sf', sample=all_samples_minus_v2)
    output:
        'output/04_multiqc/multiqc_report.html'
    params:
        outdir = 'output/04_multiqc',
        indirs = ['output/03_salmon/', 'output/03_salmon_MhV1_only']
    log:
        'output/logs/multiqc.log'
    container:
        multiqc_container
    shell:
        'multiqc '
        '-f '
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'

##########################
## salmon quant Mh MhV1 ##
##########################

rule salmon_quant:
    input:
        index = 'output/03_salmon/transcripts_index/refseq.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/03_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/03_salmon/transcripts_index',
        outdir = 'output/03_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/03_salmon/salmon_logs/salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.trimmed_r1} '
        '-2 {input.trimmed_r2} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule salmon_index:
    input:
        transcripts = 'output/02_concat_Mh_MhV1_genomes/Mh_MhV1_transcripts.fa'
    output:
        'output/03_salmon/transcripts_index/refseq.bin'
    params:
        outdir = 'output/03_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcripts} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

############################
## salmon quant MhV1 only ##
############################

rule Mhv1_salmon_quant:
    input:
        index = 'output/03_salmon_MhV1_only/transcripts_index/refseq.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/03_salmon_MhV1_only/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/03_salmon_MhV1_only/transcripts_index',
        outdir = 'output/03_salmon_MhV1_only/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/03_salmon_MhV1_only/salmon_logs/salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.trimmed_r1} '
        '-2 {input.trimmed_r2} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule Mhv1_salmon_index:
    input:
        transcripts = 'data/MhV1_prodigal/renamed_nucleotide_seq.fasta'
    output:
        'output/03_salmon_MhV1_only/transcripts_index/refseq.bin'
    params:
        outdir = 'output/03_salmon_MhV1_only/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon_index_Mhv1_only.log'
    shell:
        'salmon index '
        '-t {input.transcripts} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

####################
## concat genomes ##
####################

rule concat_Mh_MhV1_genomes:
    input:
        Mh_genome = 'data/final_microctonus_assemblies_annotations/Mh.scaffolds.virus_removed.fa',
        Mhv1_genome = 'data/MhV1_genome.fa',
        Mh_gtf = 'output/01_gtfs/Mh_annots.gtf',
        Mhv1_gtf = 'output/01_gtfs/MhV1_gene_predictions.gtf',
        Mh_transcripts = 'data/final_microctonus_assemblies_annotations/Mh.mrna-transcripts.virus_removed.fa',
        Mhv1_transcripts = 'data/MhV1_prodigal/renamed_nucleotide_seq.fasta'
    output:
        genome_concat = 'output/02_concat_Mh_MhFV_genomes/Mh_MhFV.fa',
        gtf_concat = 'output/02_concat_Mh_MhFV_genomes/Mh_MhFV.gtf',
        transcripts_concat = 'output/02_concat_Mh_MhFV_genomes/Mh_MhFV_transcripts.fa'
    shell:
        'cat {input.Mh_genome} {input.Mhv1_genome} > {output.genome_concat} & '
        'cat {input.Mh_gtf} {input.Mhv1_gtf} > {output.gtf_concat} & '
        'cat {input.Mh_transcripts} {input.Mhv1_transcripts} > {output.transcripts_concat} & ' 
        'wait'

rule convert_MhV1_gff_to_gtf:
    input:
        'data/MhV1_prodigal/renamed_gene_predictions.gff'
    output:
        'output/01_gtfs/MhFV_gene_predictions.gtf'
    log:
        'output/logs/convert_MhV1_gff_to_gtf.log'
    shell:
        'gffread '
        '{input} '
        '-T '
        '-o {output} '
        '&> {log}'

rule convert_Mh_gff_to_gtf:
    input:
        'data/final_microctonus_assemblies_annotations/Mh.gff3.virus_removed.gff3'
    output:
        'output/01_gtfs/Mh_annots.gtf'
    log:
        'output/logs/convert_gff_to_gtf.log'
    shell:
        'gffread '
        '{input} '
        '-T '
        '-o {output} '
        '&> {log}'

