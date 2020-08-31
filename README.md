## TLDR

```
git clone git@github.com:aselewa/dropseqRunner.git
cd dropseqRunner
conda env create -f environment.yaml
conda activate dropRunner

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $OUTDIR --genomeFastaFiles  $FASTA --sjdbGTFfile $ANNOTATION_GTF 
python dropRunner.py  --R1 path/to/{}.R1.fastq.gz \
                      --R2 path/to/{}.R2.fastq.gz \
                      --indices $OUTDIR \
                      --protocol drop/10x-v3/10x-v2
                      --sample pbmc_v3
```

If the above give you any trouble, run the demo to ensure everything is installed properly:
```
make run_test_workflow
```
Look for a message at the end that tells you whether the demo ran properly or not.

## Getting started

dropRunner is a Snakemake-based pipeline for processing single-cell RNA-seq data from the Drop-seq and 10x platform. We utilize STARsolo for alignment and constructing the digital expression matrix. We also supply a detailed report in HTML format that shows the sequencing statistics, as well as read distribution across the genome. The pipeline only works on Linux systems (excluding Windows linux subsystem).

**This pipeline is still under active development. If you have issues, please report them via GitHub.** 

### Setting up conda

You may **skip** this if you already have conda installed and configured.
 
`miniconda3` is a light version of `Anaconda`. To install on 64Linux, do the following:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Once it is done, initialize conda so it is in your path:

```
conda init bash
source ~/.bashrc
```

Check conda works:

```
conda --version
```

### 0. Set up and activate environment

Use the provided `environment.yaml` file to set up the conda environment.

```
git clone git@github.com:aselewa/dropseqRunner.git
cd dropseqRunner
conda env create -f environment.yaml
```
This may take some time depending on your environment. A fresh conda installation should take about 10-15 minutes. 

**Activate the environment before running the next steps**

```
conda activate dropRunner
```

### 1. Make reference genome indices

Use `STAR` to create indices for your reference genome of interest. You will need two things:

* fasta file of reference genome
* reference genome GTF annotations

You can get both of these from [GENCODE](https://www.gencodegenes.org/human/) for humans.

```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $OUTDIR --genomeFastaFiles  $FASTA --sjdbGTFfile $ANNOTATION_GTF 
```

### 2. Run the pipeline

Use `dropRunner.py` on your fastq files to generate count matrices. Use the `protocol` parameter to specify drop, 10x-v2, or 10x-v3. The last two are version 2 and version 3 10x platforms.

```
python dropRunner.py  --R1 path/to/{}.R1.fastq.gz
                      --R2 path/to/{}.R2.fastq.gz
                      --indices $OUTDIR
                      --cluster
                      --sample my_example_project
```

**Note 1**: You can supply multiple R1s and R2s by passing a comma-delimited list. I find this bash command useful:

```
R1=$(ls *.R1.fastq.gz | paste -sd,)
```

### 3. Output

There are two pieces of information that most users will need:

* html reports
* count matrices

The html report is in `reports/`. The count matrices are in `output/{project_name}_Solo.out`. There are two types of count matrices: `filtered` and `raw`. The raw matrix contains all valid barcodes, while filtered contains only barcodes with a certain number of UMI. This threshold is determined by `STARsolo` using a hueristic approach.

Please report any issues you run into via GitHub.
