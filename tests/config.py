#!/usr/bin/env python

import os

src_dir = os.path.dirname(os.path.realpath(__file__))

def resources_dir():
    return os.path.join(src_dir, 'resources')

def expected_dir():
    return os.path.join(src_dir, 'results_expected')

def work_dir():
    return os.path.join(src_dir, 'results_to_compare')

def mtx_dir():
    return os.path.join(work_dir(), 'output/pbmc_test_v3_Solo.out/Gene/filtered')

def fasta_file():
    return os.path.join(resources_dir(), 'chr22.fa')

def gtf_file():
    return os.path.join(resources_dir(), 'hg38.chr22.ensGene.gtf')

def fastq_files():
    return [os.path.join(resources_dir(), 'pbmc_test_v3_R1.fastq.gz'), os.path.join(resources_dir(), 'pbmc_test_v3_R2.fastq.gz')]

def indeces():
    return os.path.join(work_dir(), 'test_STAR_indeces/')
    
def overhang():
    return 99

def Nbases():
    return 10

def protocol():
    return '10x-v3'

def sample():
    return 'test'