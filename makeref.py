
import os
import shutil
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=str, help='Genome fasta file in gzip format.')
    parser.add_argument('--gtf', type=str, help='Genome annotations in GTF format.')
    parser.add_argument('--outDir', type=str, help='Output directory of of genome index')
    parser.add_argument('--cluster', action='store_true', help='Add this flag if the job should run on the cluster.')
    args = parser.parse_args()
    
    assert os.file.exists(args.fasta), 'Please provide a valid fasta file.' 
    assert os.file.exists(args.gtf), 'Please provide a valid gtf file.' 
    assert not os.dir.exists(args.outDir), \
    'Cannot create output directory because it already exists. Please provide the location and name of a non-existing directory.')
    assert shutil.which('gtfToGenePred') is not None, \
    'gtfToGenePred not found. Did you forget to activate the conda environment? Use the conda environment in environment.yaml to quickly install all the required software.' 
    
    print('Setting up directory and creating auxililary files..')
    
    os.system(f'mkdir {args.outDir}')
    os.system(f'gtfToGenePred {args.gtf} tmp -genePredExt')
    cmd = f"""awk '{{print $12"\t"$0}}' tmp | cut -f1-11 > {args.outDir}/refFlat_for_picard.refFlat; rm tmp"""
    os.system(cmd)
     
    if os.stat(f'{args.outDir}/refFlat_for_picard.refFlat').st_size == 0:
        msg = 'FAILED to create auxililary files. Check above for error messages. Make sure that you supplied a valid GTF file.'
        raise Exception(msg)

    if args.cluster:
      
        cmd =f"""#!/bin/bash

#SBATCH --job-name=genome_index
#SBATCH --output=genome_index.log
#SBATCH --error=genome_index.error
#SBATCH --time=10:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=50G
#SBATCH --tasks-per-node=4

source activate dropRunner
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {args.outDir}/ --genomeFastaFiles {args.fasta} --sjdbGTFfile {args.gtf} --sjdbOverhang 59
"""

        with open('generate_indices.sbatch', 'w') as f:
            f.write(cmd)
        
        os.system('sbatch generate_indices.sbatch')
        print('Index generation has been submitted to the cluster. Type qstat -u CNETID to check the status.')
    
    else:
      
      print('Genome index generation will run locally on this machine. This may not complete due to STARs large memory requirement.')
      os.system(f'STAR --runThreadN 1 --runMode genomeGenerate --genomeDir {args.outDir}/ --genomeFastaFiles {args.fasta} --sjdbGTFfile {args.gtf} --sjdbOverhang 59')
