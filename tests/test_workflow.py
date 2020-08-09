
import os
import shutil

from tests import config
from droprunner import configurator

def test_conda():
    os.system('source activate dropRunner')
    assert shutil.which('snakemake') is not None, 'test conda failed - snakemake not found'
    assert shutil.which('STAR') is not None, 'test conda failed - STAR not found'

def test_makeref():
    cmd=f"""
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {config.indeces()} --genomeFastaFiles  {config.fasta_file()} --sjdbGTFfile {config.gtf_file()} --sjdbOverhang {config.overhang()} --genomeSAindexNbases {config.Nbases()}
    """
    res = os.system(cmd)
    assert res == 0, 'makeref failed'

def test_configurator():
    try:
        configurator.main(R1=config.fastq_files()[0], 
                        R2=config.fastq_files()[1], 
                        indeces=config.indeces(),
                        protocol=config.protocol(), 
                        cluster=None, 
                        sample=config.sample(), 
                        work_dir=config.work_dir())  
    except:
        assert False, 'main failed'

if __name__ == "__main__":
    os.system(f'mkdir {config.work_dir()}')
    os.system(f'mkdir {config.indeces()}')
    test_conda()
    test_makeref()
    test_configurator()
