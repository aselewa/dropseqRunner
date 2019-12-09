# Snakefile for 10X processing
#  Alan Selewa & Sebastian Pott, UChicago 2019

import glob
import os

#Scripts
scripts = config["scripts"]

# 10x barcodes
barcodes = config['10x_barcodes']

#genome_index
GenomeIndex = config["genome_index"]

#gene annotations
ref_flat=config["refFlat"]

pd = config["proj_dir"]
output = pd + "output/"
fastqc_dir = output + "fastqc/"
qc_data = output + "qc_data/"
reports = pd + "reports/"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)

samples = set(glob_wildcards(pd + "_{samples}.R1.fastq.gz").samples)
read_num = ['1','2']

localrules: index_bam, make_report

rule all:
    input:
        expand(fastqc_dir + "_{sample}.R{readn}_fastqc.html", sample=samples, readn=read_num),
        expand(output + "{sample}_Aligned.sortedByCoord.out.bam", sample = samples),
        expand(output + "{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples),
        expand(qc_data + "{sample}_RNAmetrics.picard.txt", sample = samples),
        expand(reports + "{sample}/{sample}_pipeline_report.html", sample = samples)

#fastqc will be run on both input files
rule fastqc:
    input:
        pd + "_{sample}.R{read_num}.fastq.gz"
    output:
        fastqc_dir + "_{sample}.R{read_num}_fastqc.html",
        fastqc_dir + "_{sample}.R{read_num}_fastqc.zip"
    params:
        outdir = fastqc_dir
    shell:
        "fastqc -o {params.outdir} {input}"

rule unzip_whitelist:
    input:
        barcodes
    output:
        output + "barcodes_for_star.txt"
    shell:
        "gunzip -c {input} > {output}"

rule align:
    input:
        bc_read = pd + "_{sample}.R1.fastq.gz",
        cDNA_read = pd + "_{sample}.R2.fastq.gz",
        ref_genome = GenomeIndex,
        whitelist = output + "barcodes_for_star.txt"
    output:
        bam = output + "{sample}_Aligned.sortedByCoord.out.bam",
        unsorted = output + "{sample}_Aligned.out.bam"
    params:
        CBstart = 1,
        CBlen = 16,
        UMIstart = 17,
        UMIlen = 12,
        multimap = 1,
        threads = 8,
        strand = "Forward"
    shell:
        """
        STAR --runThreadN {params.threads} --genomeDir {input.ref_genome} --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outStd BAM_SortedByCoordinate --soloType Droplet --soloCBwhitelist {input.whitelist} \
            --soloCBstart {params.CBstart} --soloCBlen {params.CBlen} --soloUMIstart {params.UMIstart} --soloUMIlen {params.UMIlen} \
            --soloStrand {params.strand} --soloFeatures Gene --soloUMIdedup 1MM_Directional \
            --soloOutFileNames Solo.out/ "genes.tsv" "barcodes.tsv" "matrix.mtx" "matrixGeneFull.mtx" \
            --readFilesIn {input.cDNA_read} {input.bc_read} --readFilesCommand zcat --outFilterMultimapNmax {params.multimap} --outFileNamePrefix output/{wildcards.sample}_ --limitBAMsortRAM 48000000000 > {output.bam}
        """

rule index_bam:
    input:
        output + "{sample}_Aligned.sortedByCoord.out.bam"
    output:
        output + "{sample}_Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule collect_rna_metrics:
    input:
        bam = output + "{sample}_Aligned.sortedByCoord.out.bam",
        ref_flat = ref_flat
    output:
        rna_metrics = qc_data + "{sample}_RNAmetrics.picard.txt"
    params:
        strand = "SECOND_READ_TRANSCRIPTION_STRAND",
        memory = "-Xmx12G"
    shell:
        "picard CollectRnaSeqMetrics {params.memory} I={input.bam} O={output.rna_metrics} REF_FLAT={input.ref_flat} STRAND={params.strand}"

rule make_report:
     input:
         qc_data + "{sample}_RNAmetrics.picard.txt",
     output:
         reports + "{sample}/{sample}_pipeline_report.html"
     shell:
         """ R -e "rmarkdown::render(input = '{scripts}report_10x.Rmd', knit_root_dir='{pd}', output_file='{output}', intermediates_dir='{reports}{wildcards.sample}/', params=list(sampleID='{wildcards.sample}'))" """

