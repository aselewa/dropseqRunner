#Snakefile for Dropseq analysis based on a fixed (i.e. not automatically inferred) number of cells

import glob
import os

# Configuration ----------------------------------------------------------------
#these things should go to the config file
#configfile: "config.yaml"

#Scripts
scripts = config["scripts"]
#expected number of cells (upper limit)
cell_num = config["cell_num"]

#cell barcode UMI structure
barcode = config["barcode"]

#genome_index
GenomeIndex = config["genome_index"]
#gene file
#txn_file = config["txn_file"]
ref_flat=config["refFlat"]

pd = config["proj_dir"]
fastq = pd + ".fastq/"

output = pd + "output/"
solo_out = output + "Solo.out/"
fastqc_dir = output + "fastqc/"
qc_data = output + "qc_data/"
cell_stats = pd + "cell_stats/"
reports = pd + "reports/"

#make sure the project directory actually exists
#assert os.path.exists(pd), "Project directory exists"

# Directory to send log files. Needs to be created manually since it
# is not a file created by a Snakemake rule.
dir_log = config["dir_log"]
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)
# input data # might need to be changed to be universal

samples = set(glob_wildcards(fastq + "{samples}_R1_001.fastq.gz").samples)
read_num = ['1','2']

rule all:
    input:
        expand(fastqc_dir + "{sample}_R{readn}_001_fastqc.html", sample=samples, readn=read_num),
        expand(cell_stats + "{sample}_whitelist.txt", sample=samples),
        expand(cell_stats + "{sample}_whitelist_for_solo.txt", sample=samples),
        expand(output + "{sample}_Aligned.sortedByCoord.out.bam", sample = samples),
        expand(output + "{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples),
        expand(qc_data + "{sample}_RNAmetrics.picard.txt", sample = samples),
        expand(reports + "{sample}_pipeline_report.html", sample = samples)

#fastqc will be run on both input files
rule fastqc:
    input:
        fastq + "{sample}_R{read_num}_001.fastq.gz"
    output:
        fastqc_dir + "{sample}_R{read_num}_001_fastqc.html",
        fastqc_dir + "{sample}_R{read_num}_001_fastqc.zip"
    params:
        outdir = fastqc_dir
    shell:
        "fastqc -o {params.outdir} {input}"

rule umi_create_whitelist:
    input:
        fastq + "{sample}_R1_001.fastq.gz"
    output:
        cell_stats + "{sample}_whitelist.txt"
    params:
        cell_num = cell_num,
        bc = barcode
    shell:
        "umi_tools whitelist --stdin {input} --bc-pattern={params.bc} --set-cell-number={params.cell_num} --extract-method=string --log2stderr > {output}"

rule whitelist_for_solo:
    input:
        cell_stats + "{sample}_whitelist.txt"
    output:
        cell_stats + "{sample}_whitelist_for_solo.txt"
    shell:
        "less {input} | cut -f1 > {output}"

rule align:
    input:
        bc_read = fastq + "{sample}_R1_001.fastq.gz",
        cDNA_read = fastq + "{sample}_R2_001.fastq.gz",
        ref_genome = GenomeIndex,
        whitelist = cell_stats + "{sample}_whitelist_for_solo.txt"
    output:
        bam = output + "{sample}_Aligned.sortedByCoord.out.bam"
    params:
        CBstart = 1,
        CBlen = 12,
        UMIstart = 13,
        UMIlen = 8,
        multimap = 1,
        threads = 8,
        strand = "Forward"
    shell:
        """
        STAR --runThreadN {params.threads} --genomeDir {input.ref_genome} --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outStd BAM_SortedByCoordinate --soloType Droplet --soloCBwhitelist {input.whitelist} \
            --soloCBstart {params.CBstart} --soloCBlen {params.CBlen} --soloUMIstart {params.UMIstart} --soloUMIlen {params.UMIlen} \
            --soloStrand {params.strand} --soloFeatures Gene GeneFull SJ Velocyto --soloUMIdedup 1MM_Directional \
            --soloOutFileNames Solo.out/ "genes.tsv" "barcodes.tsv" "matrix.mtx" "matrixGeneFull.mtx" \
            --readFilesIn {input.cDNA_read} {input.bc_read} --readFilesCommand zcat --outFilterMultimapNmax {params.multimap} --outFileNamePrefix output/{wildcards.sample}_ > {output.bam}
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
         reports + "{sample}_pipeline_report.html"
     resources:
         load=1
     shell:
         """ R -e "rmarkdown::render(input = '{scripts}dropRunner.Rmd', knit_root_dir='{pd}', output_file='{output}', params=list(sampleID='{wildcards.sample}'))" """

