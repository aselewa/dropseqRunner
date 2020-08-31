#!/usr/bin/env python

import os
import shutil
import argparse

src_dir = os.path.dirname(os.path.realpath(__file__))
install_dir = os.path.dirname(src_dir)

protocol_map = {'drop': 'Snakefile_drop.smk', 
                '10x-v2': 'Snakefile_10x.smk', 
                '10x-v3': 'Snakefile_10x.smk'}

def check_gzip(files):
    '''
    Checks if all given files are in gzip format
    '''
    for f in files:
        if f.split('.')[-1] != 'gz':
            return False
    return True

def make_config(indices, protocol, sample, work_dir):
     
    # make sure no trailing slash
    if indices[-1] == '/':
        indices = indices[:-1]
    if work_dir[-1] == '/':
        work_dir = work_dir[:-1]

    if protocol == '10x-v2':
        barcodes = f'{install_dir}/barcodes/whitelist.v2.txt.gz'
    elif protocol == '10x-v3':
        barcodes = f'{install_dir}/barcodes/whitelist.v3.txt.gz'
    else:
        barcodes = None

    config=f"""proj_dir: {work_dir}/
genome_index: {indices}/
scripts: {install_dir}/scripts/
cell_num: 10000
barcode: "CCCCCCCCCCCCNNNNNNNN"
dir_log: log/
Sample: {sample}
Protocol: {protocol}
10x_barcodes: {barcodes}
    """
    with open(os.path.join(work_dir, 'config.yaml'), 'w') as f:
         f.write(config)

def make_submit_snakemake(snakemake_path, config_path, cluster_json_path):
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

snakemake \\
    -kp \\
    --ri \\
    -j 150 \\
    --latency-wait 150 \\
    --cluster-config {cluster_json_path} \\
    -c "sbatch \\
        --mem={{cluster.mem}} \\
        --nodes={{cluster.n}} \\
        --tasks-per-node={{cluster.tasks}} \\
        --partition={{cluster.partition}} \\
        --account={{cluster.account}} \\
        --job-name={{cluster.name}} \\
        --output={{cluster.logfile}}" \\
    -s {snakemake_path} \\
    --configfile {config_path}
"""
    outpath = os.path.join(os.path.dirname(config_path),'submit_snakemake.sbatch')
    with open(outpath, 'w') as f:
        f.write(cmd)

def main(R1, R2, indices, protocol, cluster, sample, work_dir):
    
    try:
        smk = protocol_map[protocol]
    except:
        smk = protocol_map['hg38']

    SNAKEMAKE_DIR = os.path.join(install_dir, 'droprunner/')
    SMK_PATH = os.path.join(os.path.join(SNAKEMAKE_DIR, smk))
    CONFIG_PATH = os.path.join(work_dir, 'config.yaml')
    CLUSTER_CONFIG_PATH = os.path.join(install_dir, 'cluster_solo.json')

    make_config(indices, protocol, sample, work_dir)
    make_submit_snakemake(SMK_PATH, CONFIG_PATH, CLUSTER_CONFIG_PATH)

    FASTQ_DIR = os.path.join(work_dir, 'fastq')
    os.system(f'mkdir -p {FASTQ_DIR}')
    
    r1, r2 = R1.split(','), R2.split(',') 
    assert len(r1) == len(r2), \
    'Number of files in read 1 and read 2 are not the same. Please provide a read 1 and read 2 file for each experiment.'
    if check_gzip(r1) and check_gzip(r2):

        for R1,R2 in zip(r1, r2):
            f1_name = R1.split('/')[-1]
            f2_name = R2.split('/')[-1]
            
            assert 'R1' in f1_name, 'R1 filename does not contain "R1". Did you give R2 twice?'
            assert 'R2' in f2_name, 'R2 filename does not contain "R2". Did you give R1 twice?'

            os.system(f'ln -s {os.path.abspath(R1)} {FASTQ_DIR}/{f1_name}')
            os.system(f'ln -s {os.path.abspath(R2)} {FASTQ_DIR}/{f2_name}')
    else:
        msg = 'File format not recognized. Please make sure you provide gzipped fastq files (files should end with .fastq.gz)'
        raise TypeError(msg)

    if cluster:
        os.system('sbatch submit_snakemake.sbatch')
        print('Snakemake job has been submitted to the cluster.\nType "qstat -u CNETID" to see the progress of the snakefile.')
        
    else: 
        print('Running snakemake directly on this node. This may not finish because alignment requires >30GB of RAM.')
        os.system(f'snakemake -kp --ri -s {SMK_PATH} --configfile {CONFIG_PATH}')
