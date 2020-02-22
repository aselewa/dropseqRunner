## TLDR

```
git clone git@github.com:aselewa/dropseqRunner.git
cd dropseqRunner
conda env create -f environment.yaml
conda activate dropRunner
python makeref.py --fasta path/to/fasta \
                  --gtf /path/to/gtf \
                  --outDir name_of_output_dir
python dropRunner.py  --R1 path/to/{}.R1.fastq.gz \
                      --R2 path/to/{}.R2.fastq.gz \
                      --indices path/to/makeref/indices \
                      --protocol drop
```

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

Use `makeref.py` to create indices for your reference genome of interest. You will need two things:

* fasta file of reference genome
* reference genome GTF annotations

You can get both of these from [GENCODE](https://www.gencodegenes.org/human/) for humans.

```
python makeref.py --fasta refgenome.fa
                  --gtf annots.gtf
                  --outDir myref_indices
                  --cluster
```

This will run on the RCC midway2 using the `broadwl` partition. If you are not at UChicago, do not include the cluster flag. 
This command will create a folder called `myref_indices`. You will need this folder in the next step.

### 2. Run the pipeline

Use `dropRunner.py` on your fastq files to generate count matrices. Use the `protocol` parameter to specify drop, 10x-v2, or 10x-v3. The last two are version 2 and version 3 10x platforms.

```
python dropRunner.py  --R1 path/to/{}.R1.fastq.gz
                      --R2 path/to/{}.R2.fastq.gz
                      --indices myref_indices
                      --cluster
                      --sample my_example_project
```

**Note 1**: You can supply multiple R1s and R2s by passing a comma-delimited list. I find this bash command useful:

```
R1=$(ls *.R1.fastq.gz | paste -sd,)
```

**NOTE 2**: Make sure your fastq files match the following pattern: **{project_name}_R1.fastq.gz** where {project_name} is a unique identifider.

### 3. Output

There are two pieces of information that most users will need:

* html reports
* count matrices

The html report is in `reports/`. The count matrices are in `output/{project_name}_Solo.out`. There are two types of count matrices: `filtered` and `raw`. The raw matrix contains all valid barcodes, while filtered contains only barcodes with a certain number of UMI. This threshold is determined by `STARsolo` using a hueristic approach.

Please report any issues you run into via GitHub.
