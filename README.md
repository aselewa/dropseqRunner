## Getting started

To use this pipeline, you will need `Anaconda` or `miniconda` installed. Python 3.x is also required but that should have been installed with `conda`.

### 0. Set up environment

Use the provided `environmeny.yaml` file to set up the conda environment.

```
conda env create -f environment.yaml
```

Once it is finished, confirm it is installed:

``
source activate dropRunner
``

### 1. Make reference genome indeces

Use `makeref.py` to create indeces for your reference genome of interest. You will need two things:

* fasta file of reference genome
* reference genome GTF annotations (you can get this from ENCODE)

```
python makeref.py --fasta refgenome.fa --gtf annots.gtf --outDir myref_indeces --cluster yes
```

This will run on SLURM using the `broadwl` partition. If you are not at UChicago, pass no to the cluster flag: `--cluster no`
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


**NOTE 1**: Paths to fastq files must be absolute paths. You may give multiple, comma-delimited fastq files for parallel processing. 

**NOTE 2**: Make sure your fastq files match the following pattern: **{project_name}_001.R1.fastq.gz** where {project_name} is a unique identifider.

