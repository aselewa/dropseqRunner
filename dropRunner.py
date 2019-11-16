import os
import argparse

install_dir = os.path.dirname(os.path.realpath(__file__))
work_dir = os.getcwd()

def check_gzip(files):
    '''
    Checks if all given files are in gzip format
    '''
    for f in files:
        if f.split('.')[-1] != 'gz':
            return False
    return True

def make_config(args, install_dir, work_dir):

    config=f"""proj_dir: {work_dir}/
genome_index: {args.indices}/
refFlat: {args.indices}/refFlat_for_picard.refFlat
scripts: {install_dir}/Scripts/
cell_num: 10000
barcode: "CCCCCCCCCCCCNNNNNNNN"
dir_log: log/
Sample: {args.sample}
Protocol: {args.protocol}
    """
    with open('config.yaml', 'w') as f:
         f.write(config)

def make_submit_snakemake(install_dir, work_dir):
    '''
    Creates the SLURM that submits Snakefile on the cluster
    '''
    cmd =f"""#!/bin/bash

#SBATCH --job-name=snakemake
#SBATCH --output=snakelog.out
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

source activate dropRunner

snakemake \\
    -kp \\
    --ri \\
    -j 150 \\
    --latency-wait 150 \\
    --cluster-config {install_dir}/cluster_solo.json \\
    -c "sbatch \\
        --mem={{cluster.mem}} \\
        --nodes={{cluster.n}} \\
        --tasks-per-node={{cluster.tasks}} \\
        --partition=broadwl \\
        --job-name={{cluster.name}} \\
    --output={{cluster.logfile}}" \\
    -s {install_dir}/Snakefile_solo.smk \\
    --configfile {work_dir}/config.yaml
"""

    with open('submit_snakemake.sbatch', 'w') as f:
        f.write(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--R1', type=str, help='Absolute path to gzipped read 1 fastq file (comma-delimited list of files if multiple.) REQUIRED')
    parser.add_argument('--R2', type=str, help='Absolute path to gzipped read 2 fastq file (comma-delimited list of files if multiple.) REQUIRED')
    parser.add_argument('--indices', type=str, help='Indeces folder made by makeref.py')
    parser.add_argument('--rurun', action='store_true', help='This flag re-runs a previously failed attempt.')
    parser.add_argument('--protocol', type=str, help='Protocol for producing data. Currently only drop-seq is available. Default: drop')
    parser.add_argument('--cluster', action='store_true', help='Provide this flag if this job should be run on the cluster.')
    parser.add_argument('--sample', type=str, help='sample name. Optional.')
    
    args = parser.parse_args()
    if args.rerun:
         if os.file.exists('submit_snakemake.sbatch'):
             os.system('sbatch submit_snakemake.sbatch')   
         else:
             raise Exception('sbatch file not found. Are you sure you ran this pipeline before?')

    if args.protocol == None:
        args.protocol = 'drop'
    if args.cluster == None:
        args.cluster = 'no'
    if args.sample == None:
        args.sample = 'NAME_NOT_PROVIDED'
        
    if args.R1 == None or args.R2 == None:
      raise Exception('Please provide gzipped fastq files for read 1 and read 2.')
    
    if args.indices == None:
        raise Exception('Please provide indices made by the makeref.py function!')
      
    os.system('mkdir .fastq')    

    r1, r2 = args.R1.split(','), args.R2.split(',')
    if not os.path.isabs(r1[0]) or not os.path.isabs(r2[0]):
      msg = 'Please give absolute path to the fastq files. No relative paths!'
      raise Exception(msg)
      
    if len(r1) != len(r2):
        msg='Number of files in read 1 and read 2 are not the same. Please provide a read 1 and read 2 file for each experiment.'
        raise Exception(msg)

    if check_gzip(r1) and check_gzip(r2):
        
        for R1,R2 in zip(r1, r2):
            f1_name = R1.split('/')[-1]
            f2_name = R2.split('/')[-1]
            os.system(f'ln -s {R1} .fastq/{f1_name}')
            os.system(f'ln -s {R2} .fastq/{f2_name}')
    else:
        msg = 'File format not recognized. Please make sure you provide gzipped fastq files (files should end with .fastq.gz)'
        raise TypeError(msg)

    make_config(args, install_dir, work_dir)
    
    if args.cluster:
        make_submit_snakemake(install_dir, work_dir)
        print('Submitting snakemake job..')
        os.system('sbatch submit_snakemake.sbatch')
        print('Snakemake job has been submitted to the cluster.\nType "qstat -u CNETID" to see the progress of the snakefile.')
        
    else:
        print('Running snakemake directly on this node. This may not finish because alignment requires >30GB of RAM.')
        os.system(f'source activate dropRunner; snakemake -kp --ri -s {install_dir}/Snakefile_solo.smk --configfile {work_dir}/config.yaml')



