## Getting started

To use this pipeline, you will need `Anaconda` or `miniconda` installed. Python 3.x is also required but that should have been installed with `conda`. Note that this pipeline has only been tested on 64bit Linux. Use on MacOS will likely result in errors, while Windows is completely unsupported.

### Setting up conda

You may **skip** this if you already have conda installed.
 
`miniconda3` is a light version of `Anaconda`. To install on 64Linux, do the following:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Make sure that `miniconda3/bin/` and `miniconda3/condabin/` are in your `PATH`! They are usually added to your `.bashrc` by the conda installer automatically. To apply the changes, type:

```
source ~/.bashrc
```

Check your path:

```
echo $PATH
```

If you don't see something with `../miniconda3/bin/` nor `../miniconda3/condabin/`, then you need to add these manually by typing:

```
export PATH=path/to/miniconda3/bin/:$PATH
export PATH=path/to/miniconda4/condabin/:$PATH
```

### 0. Set up environment

Use the provided `environmeny.yaml` file to set up the conda environment.

```
git clone git@github.com:aselewa/dropseqRunner.git
cd dropseqRunner
conda env create -f environment.yaml
```

This may take some time depending on your environment. A fresh conda installation should take about 10-15 minutes. 

Once it is finished, proceed to the next step. You do not need to activate the environment.

### 1. Make reference genome indeces

Use `makeref.py` to create indeces for your reference genome of interest. You will need two things:

* fasta file of reference genome
* reference genome GTF annotations

You can get both of these from [GENCODE](https://www.gencodegenes.org/human/) for humans.

```
python makeref.py --fasta refgenome.fa 
                   --gtf annots.gtf 
                   --outDir myref_indeces 
                   --cluster yes
```

This will run on the RCC midway2 using the `broadwl` partition. If you are not at UChicago, pass no to the cluster flag: `--cluster no`
This command will create a folder called `myref_indeces`. You will need this folder in the next step.

### 2. Run the pipeline

Use `dropRunner.py` on your fastq files to generate count matrices.

```
python dropRunner.py  --R1 path/to/sample_001_R1.fastq.gz 
                      --R2 path/to/sample_001_R2.fastq.gz 
                      --indeces myref_indeces
                      --cluster yes
                      --sample my_example_project
```


Once again, this will run on the RCC midway2 using the `broadwl` partition. If you are not at UChicago, pass no to the cluster flag: `--cluster no`

**NOTE 1**: Paths to fastq files must be absolute paths. You may give multiple, comma-delimited fastq files for parallel processing. 

**NOTE 2**: Make sure your fastq files match the following pattern: **{project_name}_001_R1.fastq.gz** where {project_name} is a unique identifider.

### 3. Output

There are two pieces of information that most users will need:

* html reports
* count matrices

The html report is in `reports/`. The count matrices are in `output/{project_name}_Solo.out`. There are two types of count matrices: `filtered` and `raw`. The raw matrix contains all valid barcodes, while filtered contains only barcodes with a certain number of UMI. This threshold is determined by `STARsolo` using a hueristic approach.
